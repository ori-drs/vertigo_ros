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
#include "betweenFactorAdaptive.h"
#include <boost/optional/optional_io.hpp>

namespace vertigo {

template<class VALUE>
class PriorFactorOutlierProcess : public gtsam::NoiseModelFactor1<ShapeParameter>
{
  public:
    PriorFactorOutlierProcess() {}
    PriorFactorOutlierProcess(gtsam::Key key, const VALUE& prior, const gtsam::SharedNoiseModel& modelPsi)
    : gtsam::NoiseModelFactor1<ShapeParameter>(modelPsi,key),
      priorFactor(key, prior){}


    gtsam::Vector evaluateError(const ShapeParameter& alpha,
                                boost::optional<gtsam::Matrix&> H =  boost::none) const
      {

      // Matias's method
      cout << "[priorFactor]alpha.weight_z is: " << alpha.weight_z() << endl;
      cout << "[priorFactor]alpha.value is: " << alpha.value() << endl;

      double sqrtPsi = sqrt(Psi(alpha.value(), alpha.weight_z()));
      double dPsi_dalpha = Psi_alpha(alpha.value(), alpha.weight_z());


//      cout << "[prior factor]sqrtPsi is: " << sqrtPsi << endl;
//      cout << "[prior factor]dPsi_dalpha is: " << dPsi_dalpha << endl;

      // calculate error ||sqrt(Psi)||^2 using Rosen et al. to use non Gaussian factor
      gtsam::Vector error = gtsam::Vector1(sqrtPsi);

      // handle derivatives: 0.5 * sqrt(Psi(alpha, weight))^-0.5 * Psi_alpha
      double h = sqrtPsi > 0? 0.5 * (1.0/sqrtPsi * dPsi_dalpha) : 1E-2;
      if (H) *H = gtsam::Vector1(h);

//      std::cout << "[prior factor]H is: \n" << H << std::endl;
//      std::cout << "[prior factor]error is: \n" << error << std::endl;

      return error;
      

      /*
        // calculate error
        gtsam::Vector error = priorFactor.evaluateError(alpha, H); // how to compute error?

        std::cout << "PriorFactorOutlierProcess: weight from shape parameter: " << std::to_string(alpha.weight_z()) << std::endl;

        gtsam::Vector a = gtsam::ones(3);
        a[0] = Psi(2.0,1.0);
        a[1] = Psi_z(2.0,1.0);
        a[2] = Psi_alpha(2.0,1.0);

        gtsam::Vector b = gtsam::ones(3);
        b[1] = alpha.weight_z() - 1.0;
        b[2] = error[0];


        error = a.transpose()*b;

        // handle derivatives
        if (H) *H = *H * Psi_alpha(2.0,1.0); //what should I put here?

//          std::cout << "error is: \n" << error << std::endl;
        return error;
        */
      }

  private:
    gtsam::PriorFactor<VALUE> priorFactor;
//    vertigo::BetweenFactorAdaptive<VALUE> betweenFactorAdaptive;
    ShapeParameter alpha_;
    double epsilon_ = 1E-5;
    double Psi(double alpha, double z) const {

      double b;
      b = abs(alpha-2)+epsilon_;
      double d = alpha>=0? alpha+epsilon_ : alpha-epsilon_;

      if (alpha=2.0) return 0;
      else if (alpha<2.0) return (b/d)*((1-0.5*d)*pow(z,d/(d-2.0))+0.5*d*z-1);
    }

    double Psi_alpha(double alpha, double z) const {

      double b;
      b = abs(alpha-2)+epsilon_;
      double d = alpha>=0? alpha+epsilon_ : alpha-epsilon_;

      double Psi_b, Psi_d;

      Psi_b = ((2.0-d)*pow(z,d/(d-2.0))+d*z-2.0)/(2.0*d);

      Psi_d = (b*((pow(z,d/(d-2.0))*(log(z)-1.0)+1.0)*d+2*pow(z,d/(d-2.0))-2.0))/((d-2.0)*d*d);

      return Psi_b+Psi_d;
//      if (alpha==0) return 0;
//      else if (alpha<=-10.0) return 0;
//      else{
//        return (pow(z,alpha/(alpha-2.0))*log(z))/alpha + (2.0*pow(z,alpha/(alpha-2.0))-2.0)*(alpha*alpha) - 0.5*(pow(z,alpha/(alpha-2.0))-z);
//      }
    }

    double Psi_z(double alpha, double z) const {

      double b;
      b = abs(alpha-2)+epsilon_;
      double d = alpha>=0? alpha+epsilon_ : alpha-epsilon_;

      return -(b*(pow(z,d/(d-2.0))-z))/(2.0*z);

//      if (alpha==0) return 1.0-1.0/z;
//      else if (alpha<=-10.0) return log(z);
//      else{
//        return -((alpha-2.0)*(pow(z,alpha/(alpha-2.0))-z)/2.0*z);
//      }
    }

};
} // vertigo namespace



#endif /* OUTLIERPROCESS_H_ */

