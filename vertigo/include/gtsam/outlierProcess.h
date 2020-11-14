/*
 * outlierProcess.h
 *
 *  Created on: 08.06.2020
 *      Author: Milad Ramezani
 * Outlier process named Psi(z,alpha) in Barron's paper is a penalty term
 * added to the quatratic form of residuals when using Black-Ranagarajan
 *
 * Psi(alpha,z) by Jonathan T. Barron (https://arxiv.org/pdf/1701.03077.pdf)
 * Name                  Symbol               Formula
 * Outlier process       \Psi(z,alpha)        -log(z)+z-1 if alpha=0, zlog(z)-z+1     if alpha=-inf, |alpha-2|/alpha*[(1-0.5*alpha)z^(alpha/(alpha-2))+0.5*z*alpha-1] otherwise
 * Derivative wrt alpha  \Psi_alpha(z,alpha)  (z^(alpha/(alpha-2))*ln(z))/alpha+(2*z^(alpha/(alpha-2))-2)/(alpha^2)-0.5*(z(alpha/(alpha-2))-z)
 * Derivative wrt z      \Psi_z(z,alpha)      -(alpha-2)*(z^(alpha/(alpha-2))-z)/2*z
 * We use 1st order of Taylor expansion to approximate Psi(z,alpha)
 * Psi(z,alpha) = Psi(z0,alpha0) + Psi_alpha(z0,alpha0)*(alpha-alpha0) + Psi_z(z0,alpha0)*(z-z0)
 */


#ifndef OUTLIERPROCESS_H_
#define OUTLIERPROCESS_H_

#pragma once

#include <gtsam/slam/PriorFactor.h>
#include <gtsam/base/Lie.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include "shapeParameter.h"
//#include "betweenFactorAdaptive.h"
#include <boost/optional/optional_io.hpp>

namespace vertigo {

template<class VALUE>
class OutlierProcess : public gtsam::NoiseModelFactor1<ShapeParameter>
{
  public:
    OutlierProcess() {}
    OutlierProcess(gtsam::Key key, const VALUE& prior, const gtsam::SharedNoiseModel& modelPsi)
    : gtsam::NoiseModelFactor1<ShapeParameter>(modelPsi,key),
      priorFactor(key, prior){weight_=1.0;}


    gtsam::Vector evaluateError(const ShapeParameter& alpha,
                                boost::optional<gtsam::Matrix&> H =  boost::none) const
    {

      double sqrtPsi = sqrt(Psi(alpha.value(), weight_));
      double dPsi_dalpha = PsiDerivativeAlpha(alpha.value(), weight_);

      // Calculate error ||sqrt(Psi)||^2 using Rosen et al. to use non-Gaussian factor
      gtsam::Vector error = gtsam::Vector1(sqrtPsi);

      // Handle derivatives: 0.5 * sqrt(Psi(alpha, weight))^-0.5 * Psi_alpha
      double h = sqrtPsi > 0? 0.5 * (1.0 / sqrtPsi * dPsi_dalpha) : 1E-5;
      if (H) *H = gtsam::Vector1(h);

      // Sanity check if there are nan values
      if(!(*H).allFinite() /*|| (*H).isZero()*/){
        std::cout << "[OutlierProcess] H: " << H << std::endl;
        std::cout << "Psi(alpha=" << alpha.value() << ", z=" << weight_ << ") : " << sqrtPsi << std::endl;
        std::cout << "PsiDerivativeAlpha(alpha=" << alpha.value() << ", z=" << weight_ << ") : " << dPsi_dalpha << std::endl;
        std::cout << std::endl;
//        exit(-1);
      }

      if(!error.allFinite()){
        std::cout << "[OutlierProcess] error: " << error << std::endl;
        std::cout << "Psi(alpha=" << alpha.value() << ", z=" << weight_ << ") : " << sqrtPsi << std::endl;
        std::cout << "PsiDerivativeAlpha(alpha=" << alpha.value() << ", z=" << weight_ << ") : " << dPsi_dalpha << std::endl;
        std::cout << std::endl;
//        exit(-1);
      }

      return error;
    }

    // Interface to set the weight externally
    void setWeight(double weight) { weight_ = weight; }

  private:
    // Tolerance used for singularities
    const double epsilon_ = 1E-5;

    // Use practical implementation flag
    bool usePractical_ = false;

    // Use numerical derivatives flag
    bool useNumericalDerivatives_ = false;

    // Internal instance of prior factor
    gtsam::PriorFactor<VALUE> priorFactor;

    // Internal copy of weight
    double weight_;

    // Shape parameter
    ShapeParameter alpha_;

    double Psi(double alpha, double z) const{
      double psi;
      
      if(usePractical_){
        psi = PsiPractical(alpha, z);
      } else{
        psi = PsiAnalytical(alpha, z);
      }

      if(psi < 0.0){
        //std::cout << "Psi(alpha=" << alpha << ", z=" << z << ") : " << psi << std::endl;
        psi = 0.0;
      }

      return psi;
    }
    
    double PsiDerivativeAlpha(double alpha, double z) const{
      double dpsi;

      if(usePractical_){
        if(useNumericalDerivatives_){
          // Practical Jacobian using numerical differentiation
          double psi_m = PsiPractical(alpha - epsilon_, z); // psi(i-1)
          double psi = PsiPractical(alpha, z); // psi(i)
          dpsi = (psi - psi_m) / epsilon_;
        
        } else {
          // Analytical Jacobian using practical implementation
          dpsi = PsiPracticalDerivativeAlpha(alpha, z);
        }
        
      } else{
        if(useNumericalDerivatives_){
          // Analytical Jacobian using numerical differentiation
          double psi_m = PsiAnalytical(alpha - epsilon_, z); // psi(i-1)
          double psi = PsiAnalytical(alpha, z); // psi(i)
          dpsi = (psi - psi_m) / epsilon_;
        
        } else{
          // Analytical Jacobian
          dpsi = PsiAnalyticalDerivativeAlpha(alpha, z);
        }
      }

//      if(dpsi < epsilon_)
//        dpsi = epsilon_;


      // std::cout << "PsiDerivativeAlpha(alpha=" << alpha << ", z=" << z << ") : " << dpsi << std::endl;

      return dpsi;
    }
    
    // Psi function from Black-Rangarajan formulation, original proposal by Barron (Eq. 25, Appendix A)
    double PsiAnalytical(double alpha, double z) const {
      double psi;

      // Deal with singularities due to log(0)
      // z should always be positive since is a weighted quadratic error
      if(abs(z) < epsilon_){
        z = epsilon_;
      }

      if( abs(alpha) == 0.0){
        psi = -log(z) + z - 1.0;

      } else if(alpha <= ShapeParameter::MIN){
        psi = z * log(z) - z + 1;

       } else if(abs(alpha-2) <= 0.01){
         psi = 0.0; // this shouldn't be defined

      } else{
        psi = (abs(alpha - 2)/alpha) * ( (1 - 0.5 * alpha) * pow(z, alpha/(alpha-2)) + 0.5*alpha*z - 1 );
      }

      return psi;
    }

    double PsiAnalyticalDerivativeAlpha(double alpha, double z) const{
      double dpsi;

      // Deal with singularities due to log(0)
      // z should always be positive since is a weighted quadratic error
      if(abs(z) < epsilon_){
        z = epsilon_;
      }

      if( abs(alpha) == 0.0){
        dpsi = 0.0;

      } else if(alpha <= ShapeParameter::MIN){
        dpsi = 0.0;

      } else if(abs(alpha-2) <= 0.01){
        dpsi = 0.0;

      } else{
//        dpsi = (0.5*z + pow(z, alpha/(alpha - 2))*(-0.5*alpha + 1)*(-alpha/pow(alpha - 2, 2) + 1.0/(alpha - 2))*log(z) - 0.5*pow(z, alpha/(alpha - 2)))*abs(alpha - 2)/alpha + (0.5*alpha*z + pow(z, alpha/(alpha - 2))*(-0.5*alpha + 1) - 1)*(((alpha - 2) > 0) - ((alpha - 2) < 0))/alpha - (0.5*alpha*z + pow(z, alpha/(alpha - 2))*(-0.5*alpha + 1) - 1)*abs(alpha - 2)/pow(alpha, 2);
        dpsi = -(pow(z,alpha/(alpha-2))*log(z))/alpha + (2-2*pow(z,alpha/(alpha-2)))/(alpha*alpha) + (pow(z,alpha/(alpha-2))-z)/2;
      }

      return dpsi;
    }

    // Practical implementation of Psi function recommended by Barron (Appendix B)
    double PsiPractical(double alpha, double z) const {
      double psi;
      double b = abs(alpha - 2.0) + epsilon_;
      double d = alpha >= 0? alpha + epsilon_ : alpha - epsilon_;

      if ( alpha >= 2.0){
        psi = 0.0;

      } else if (alpha < 2.0){
        psi = (b/d) * ( (1-0.5*d) * pow(z, d/(d-2.0) ) + 0.5*d*z - 1.0);
      }

      return psi;
    }

    double PsiPracticalDerivativeAlpha(double alpha, double z) const {
      double dpsi;
      double b = abs(alpha - 2.0) + epsilon_;
      double d = alpha>=0? alpha + epsilon_ : alpha - epsilon_;

      double dPsi_db = (0.5*d*z + pow(z, d/(d - 2.0))*(-0.5*d + 1) - 1.0)/d;
      double dPsi_dd = b*(0.5*z + pow(z, d/(d - 2.0))*(-0.5*d + 1)*(-d/pow(d - 2.0, 2) + 1.0/(d - 2.0))*log(z) - 0.5*pow(z, d/(d - 2.0)))/d - b*(0.5*d*z + pow(z, d/(d - 2.0))*(-0.5*d + 1) - 1.0)/pow(d, 2);
      double db_dalpha = (((alpha - 2) > 0) - ((alpha - 2) < 0));
      double dd_dalpha = 1.0;

      dpsi = dPsi_db * db_dalpha + dPsi_dd * dd_dalpha;

      return dpsi;
    }

    double Psi_z(double alpha, double z) const {

      double b;
      b = abs(alpha-2)+epsilon_;
      double d = alpha>=0? alpha+epsilon_ : alpha-epsilon_;

      return -(b*(pow(z,d/(d-2.0))-z))/(2.0*z);
    }

};
} // vertigo namespace



#endif /* OUTLIERPROCESS_H_ */

