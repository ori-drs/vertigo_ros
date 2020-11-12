/*
 * betweenFactorAdaptive.h
 *
 *  Created on: 26.05.2020
 *      Author: Milad Ramezani
 */

#ifndef BETWEENFACTORADAPTIVE_H_
#define BETWEENFACTORADAPTIVE_H_

#include "planarSLAM.h"
#include <gtsam/nonlinear/NonlinearFactor.h>
#include "gtsam/base/LieScalar.h"
#include "gtsam/linear/NoiseModel.h"
#include "outlierProcess.h"

#include <iostream>
using std::cout;
using std::endl;

#include "shapeParameter.h"
#include "AdaptiveLossFunction.h"

namespace vertigo {

  template<class VALUE>
  class BetweenFactorAdaptive : public gtsam::NoiseModelFactor3<VALUE, VALUE, ShapeParameter>
  {
    public:
      BetweenFactorAdaptive() {}
      BetweenFactorAdaptive(gtsam::Key key1, gtsam::Key key2, gtsam::Key key3, const VALUE& measured, const gtsam::SharedNoiseModel& model, OutlierProcess<ShapeParameter>* outlierProcess)
      : gtsam::NoiseModelFactor3<VALUE, VALUE, ShapeParameter>(model, key1, key2, key3),
        betweenFactor(key1, key2, measured, model),outlierProcess_(outlierProcess) {}

      gtsam::Vector evaluateError(const VALUE& p1, const VALUE& p2, const ShapeParameter& alpha,
                                  boost::optional<gtsam::Matrix&> H1 = boost::none,
                                  boost::optional<gtsam::Matrix&> H2 =  boost::none,
                                  boost::optional<gtsam::Matrix&> H3 =  boost::none) const
        {

          // Calculate error
          gtsam::Vector error = betweenFactor.evaluateError(p1, p2, H1, H2);

          // Calculate weighted error
          double weightedError = this->noiseModel_->distance(error);
//          double weightedError = error.norm();

          // Recover weight and set it up in outlierProcess factor
          double c = 1.0; // c is scaling param set before optimisation
          double w = weightAdaptive (weightedError, alpha.value(), c);
          weight_ = w;
          outlierProcess_->setWeight(w);

          // Compute derivative of w wrt alpha
          double dw = weightAdaptiveDerivativeAlpha(weightedError, alpha.value(), c);

          // Scale factor error by the weight from the adaptive kernel
          // why this line was moved to below the Jacobians?(MR)


          // Jacobians
          // The Jacobians wrt the poses should be scaled
          if (H1) *H1 = *H1 * w;
          if (H2) *H2 = *H2 * w;
          // The Jacobian wrt alpha should be handled accordingly
          if (H3) *H3 =  error * dw;

          error *= w;
          
          // Sanity check if there are nan values
          // if(!(*H1).allFinite() || (*H1).isZero()){
          //   std::cout << "[BetweenFactorAdaptive] H1: " << H1 << std::endl;
          //   std::cout << "weightAdaptive(x=" << weightedError << ", alpha=" << alpha.value() << ") : " << w << std::endl;
          //   std::cout << "weightAdaptiveDerivativeAlpha(x=" << weightedError << ", alpha=" << alpha.value() << ") : " << dw << std::endl;
          //   std::cout << std::endl;
          //   exit(-1);
          // }
          // if(!(*H2).allFinite() || (*H2).isZero()){
          //   std::cout << "[BetweenFactorAdaptive] H2: " << H2 << std::endl;
          //   std::cout << "weightAdaptive(x=" << weightedError << ", alpha=" << alpha.value() << ") : " << w << std::endl;
          //   std::cout << "weightAdaptiveDerivativeAlpha(x=" << weightedError << ", alpha=" << alpha.value() << ") : " << dw << std::endl;
          //   std::cout << std::endl;
          //   exit(-1);
          // }

          if(!(*H3).allFinite() /*|| (*H3).isZero()*/){
            std::cout << "[BetweenFactorAdaptive] H3:\n " << H3 << std::endl;
            std::cout << "weightAdaptive(x=" << weightedError << ", alpha=" << alpha.value() << ") : " << w << std::endl;
            std::cout << "weightAdaptiveDerivativeAlpha(x=" << weightedError << ", alpha=" << alpha.value() << ") : " << dw << std::endl;
            std::cout << std::endl;
//            exit(-1);
          }
          if(!error.allFinite()){
            std::cout << "[BetweenFactorAdaptive] error: " << error << std::endl;
            std::cout << "weightAdaptive(x=" << weightedError << ", alpha=" << alpha.value() << ") : " << w << std::endl;
            std::cout << "weightAdaptiveDerivativeAlpha(x=" << weightedError << ", alpha=" << alpha.value() << ") : " << dw << std::endl;
            std::cout << std::endl;
            exit(-1);
          }

          return error;
        }

      double getWeight() const {return weight_;}

    private:
      // Tolerance used for singularities
      double epsilon_ = 1E-2;
      
      // Use practical implementation flag
      bool usePractical_ = false;

      // Use numerical derivatives flag
      bool useNumericalDerivatives_ = false;
      
      // Internal instance of non-adaptive BetweenFactor
      gtsam::BetweenFactor<VALUE> betweenFactor;
      
      // Adaptive weight
      mutable double weight_;
      
      // Hack: Pointer to outlierProcess factor
      OutlierProcess<ShapeParameter>* outlierProcess_;

    /**
    * A General and Adaptive Robust Loss Function by Jonathan T. Barron (https://arxiv.org/pdf/1701.03077.pdf)
    * which  present a generalization of the Cauchy/Lorentzian, Geman-McClure, Welsch/Leclerc,
    * generalized Charbonnier, Charbonnier/pseudo-Huber/L1-L2, and L2 loss functions
    * Name        Symbol            Adaptive Kernel
    * Residual    \rho(r,alpha,c)   1/2*(r/c)^2 if alpha=2, log(1/2*(r/c)^2+1) if alpha=0, 1-exp(-1/2*(r/c)^2)     if alpha=-inf, |alpha-2|/alpha*[((r/c)^2/|alpha-2| +1)^(alpha/2)-1] otherwise
    * Derivative  \phi(x)           r/c^2       if alpha=2, 2r/(r^2+2c^2)      if alpha=0, r/c^2*exp(-1/2*(r/c)^2) if alpha=-inf, r/c^2*((r/c)^2/|alpha-2|+1)^(alpha/2-1)              otherwise
    * Weight      w(x)=\phi(r)/r    1/c^2       if alpha=2, 2/(r^2+2c^2)       if alpha=0, 1/c^2*exp(-1/2*(r/c)^2) if alpha=-inf, 1/c^2*((r/c)^2/|alpha-2|+1)^(alpha/2-1)              otherwise
    * where r is the residual, alpha is the shape parameter and c>0 is scale param which controls the size of the loss's quadratic bowl near r=0 (default c=1)
    */
    double weightAdaptive(double x, double alpha, double c) const{
      double w;
      
      if(usePractical_){
        w = weightAdaptivePractical(x, alpha, c);
      } else{
        w = weightAdaptiveAnalytical(x, alpha, c);
      }

      // std::cout << "weightAdaptive(x=" << x << ", alpha=" << alpha << ", c=" << c << ") : " << w << std::endl;

      return w;
    }
    
    double weightAdaptiveDerivativeAlpha(double x, double alpha, double c) const{
      double dw;

      if(usePractical_){
        if(useNumericalDerivatives_){
          // Practical Jacobian using numerical differentiation
          double w_m = weightAdaptivePractical(x, alpha - epsilon_, c); // psi(i-1)
          double w = weightAdaptivePractical(x, alpha, c); // psi(i)
          dw = (w - w_m) / epsilon_;
        
        } else{
          // Analytical Jacobian using practical implementation
          dw = weightAdaptivePracticalDerivativeAlpha(x, alpha, c);
        }
        
      } else{
        if(useNumericalDerivatives_){
          // Analytical Jacobian using numerical differentiation
          double w_m = weightAdaptiveAnalytical(x, alpha - epsilon_, c); // psi(i-1)
          double w = weightAdaptiveAnalytical(x, alpha, c); // psi(i)
          dw = (w - w_m) / epsilon_;
        
        } else{
          // Analytical Jacobian
          dw = weightAdaptiveAnalyticalDerivativeAlpha(x, alpha, c);
        }
      }
      // do we need this here? (MR)
//      if(dw < epsilon_)
//        dw = 1E-2;

      // std::cout << "weightAdaptiveDerivativeAlpha(x=" << x << ", alpha=" << alpha << ", c=" << c << ") : " << dw << std::endl;

      return dw;
    }

    // Analytical expressions
    double weightAdaptiveAnalytical(double x, double alpha, double c) const {
      double w;

      if (abs(alpha-2) <= 0.01){
        w = 1.0 / pow(c, 2);

      } else if (alpha == 0.0){
        w = 2.0 / (pow(x, 2) + 2.0 * pow(c, 2));

      }else if (alpha <= ShapeParameter::MIN){
        w = 1.0 / pow(c, 2) * exp(-0.5 * pow(x/c, 2));

      } else{
        w = (1.0 / pow(c, 2)) * pow( pow(x/c, 2) / abs(alpha-2.0) + 1.0, (0.5*alpha-1));
      }

      return w;
    }

    double weightAdaptiveAnalyticalDerivativeAlpha(double x, double alpha, double c) const {
      double dw;

      if (abs(alpha-2) <= 0.01 || abs(alpha) == 0.0 || alpha <= ShapeParameter::MIN){
        dw = 1E-5;
      }
      else{
//        dw = pow(1 + pow(x, 2)/(pow(c, 2)*abs(alpha - 2)), 0.5*alpha - 1)*(0.5*log(1 + pow(x, 2)/(pow(c, 2)*abs(alpha - 2))) - pow(x, 2)*(0.5*alpha - 1)*(((alpha - 2) > 0) - ((alpha - 2) < 0))/(pow(c, 2)*(1 + pow(x, 2)/(pow(c, 2)*abs(alpha - 2)))*pow(alpha - 2, 2)))/pow(c, 2);
        dw = 0.5*weightAdaptiveAnalytical(x,alpha,c)*(pow(x/c,2)*pow(1-pow(x/c,2)/(alpha-2),-1))/(alpha-2);
      }

      return dw;
    }

    double weightAdaptivePractical(double x, double alpha, double c) const {

      double b = abs(alpha - 2.0) + epsilon_;
      double d = alpha >= 0? alpha + epsilon_ : alpha - epsilon_;
      
      double w = (1/pow(c,2)) * pow( pow(x/c, 2)/b + 1, (0.5*d - 1) );

      return w;
    }

    double weightAdaptivePracticalDerivativeAlpha(double x, double alpha, double c) const {
      
      double b = abs(alpha - 2.0) + epsilon_;
      double d = alpha >= 0? alpha + epsilon_ : alpha - epsilon_;

      double dw_db = -1.0*pow(x, 2)*pow(1 + pow(x, 2)/(b*pow(c, 2)), 0.5*d - 1)*(0.5*d - 1)/(pow(b, 2)*pow(c, 4)*(1 + pow(x, 2)/(b*pow(c, 2))));
      double dw_dd = 0.5*pow(1 + pow(x, 2)/(b*pow(c, 2)), 0.5*d - 1)*log(1 + pow(x, 2)/(b*pow(c, 2)))/pow(c, 2);
      double db_dalpha = (((alpha - 2) > 0) - ((alpha - 2) < 0));
      double dd_dalpha = 1.0;

      double dw = dw_db * db_dalpha + dw_dd * dd_dalpha;

      return dw;
    }

  };

} // vertigo namespace

#endif /* BETWEENFACTORADAPTIVE_H_ */
