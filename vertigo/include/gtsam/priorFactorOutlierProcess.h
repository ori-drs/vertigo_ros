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
class PriorFactorOutlierProcess : public gtsam::NoiseModelFactor3<VALUE, VALUE, ShapeParameter>
{
  public:
    PriorFactorOutlierProcess() {}
    PriorFactorOutlierProcess(gtsam::Key key1, gtsam::Key key2, gtsam::Key key3, const VALUE& measured, const gtsam::SharedNoiseModel& modelBetween, const gtsam::SharedNoiseModel& modelPsi)
    : gtsam::NoiseModelFactor3<VALUE, VALUE, ShapeParameter>(modelPsi, key1, key2, key3),
      betweenFactorAdaptive(key1, key2, key3, measured, modelBetween) {}


    gtsam::Vector evaluateError(const VALUE& p1, const VALUE& p2, const ShapeParameter& alpha,
                                  boost::optional<gtsam::Matrix&> H1 = boost::none,
                                  boost::optional<gtsam::Matrix&> H2 =  boost::none,
                                  boost::optional<gtsam::Matrix&> H3 =  boost::none) const
      {

      // Matias's method
      gtsam::Vector errorBetween = betweenFactorAdaptive.evaluateError(p1, p2, alpha, H1, H2);
      if (H1) *H1 = gtsam::Vector1(0.0);
      if (H2) *H2 = gtsam::Vector1(0.0);

      double weight = betweenFactorAdaptive.getWeight();
      std::cout << "PriorFactorOutlierProcess: weight: " << std::to_string(weight) << std::endl;

      std::cout << "alpha is:       " << alpha.value() << std::endl;

      double psi = Psi(alpha.value(), weight);
      double sqrtPsi = sqrt(psi);
      double dPsi_dalpha = Psi_alpha(alpha.value(), weight);

      std::cout << "Psi is:         " << psi << std::endl;
      std::cout << "sqrtPsi is:     " << sqrtPsi << std::endl;
      std::cout << "dPsi_dalpha is: " << dPsi_dalpha << std::endl;

      // calculate error ||sqrt(Psi)||^2 using Rosen et al. to use non Gaussian factor
      gtsam::Vector error = gtsam::Vector1(sqrtPsi);
      // handle derivatives: 0.5 * sqrt(Psi(alpha, weight))^-0.5 * Psi_alpha
      double h3 = sqrtPsi > 0? 0.5 * (1.0/sqrtPsi * dPsi_dalpha) : 0.0;
      if (H3) *H3 = gtsam::Vector1(h3);
      //          std::cout << "error is: \n" << error << std::endl;

      std::cout << "error is: " << error << std::endl;
      std::cout << "H1 is: " << H1 << std::endl;
      std::cout << "H2 is: " << H2 << std::endl;
      std::cout << "H3 is: " << H3 << std::endl;
      

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
    //gtsam::PriorFactor<VALUE> priorFactor;
    vertigo::BetweenFactorAdaptive<VALUE> betweenFactorAdaptive;

    double Psi(double alpha, double z) const {
      if (alpha==0) return -log(z)+z-1.0;
      else if (alpha<=-10.0) return z*log(z)+1.0;
      else{
        double A = (abs(alpha - 2.0)/alpha);
        double B = (1.0 - 0.5*alpha);
        double C = pow(z, alpha/(alpha - 2.0));
        double D = 0.5*alpha*z - 1.0;

        std::cout << "A is: " << A << std::endl;
        std::cout << "B is: " << B << std::endl;
        std::cout << "C is: " << C << std::endl;
        std::cout << "D is: " << D << std::endl;

        double result = A * (B*C + D);
        std::cout << "result is: " << result << std::endl;

        return result;
      }
    }

    double Psi_alpha(double alpha, double z) const {
      if (alpha==0) return 0;
      else if (alpha<=-10.0) return 0;
      else{
        return (pow(z,alpha/(alpha-2.0))*log(z))/alpha + (2.0*pow(z,alpha/(alpha-2.0))-2.0)*(alpha*alpha) - 0.5*(pow(z,alpha/(alpha-2.0))-z);
      }
    }

    double Psi_z(double alpha, double z) const {
      if (alpha==0) return 1.0-1.0/z;
      else if (alpha<=-10.0) return log(z);
      else{
        return -((alpha-2.0)*(pow(z,alpha/(alpha-2.0))-z)/2.0*z);
      }
    }

};
} // vertigo namespace



#endif /* OUTLIERPROCESS_H_ */

