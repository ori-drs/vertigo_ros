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

namespace vertigo {

template<class VALUE>
class PriorFactorOutlierProcess : public gtsam::NoiseModelFactor1<ShapeParameter>
{
  public:
    PriorFactorOutlierProcess() {}
    PriorFactorOutlierProcess(gtsam::Key key1, const VALUE& prior, const gtsam::SharedNoiseModel& model)
    : gtsam::NoiseModelFactor1<ShapeParameter>(model, key1),
      priorFactor(key1, prior, model) {}

    gtsam::Vector evaluateError(const ShapeParameter& alpha,
        boost::optional<gtsam::Matrix&> H = boost::none) const
      {

        // calculate error
        gtsam::Vector error = priorFactor.evaluateError(alpha, H); // how to compute error?
        error = Psi(2.0,1.0)+Psi_z(2.0,1.0)*(alpha.weight_z()-1.0)+Psi_alpha(2.0,1.0)*error;

        // handle derivatives
        if (H) *H = *H * Psi_alpha(2.0,1.0); //what should I put here?

//          std::cout << "error is: \n" << error << std::endl;
        return error;
      }

  private:
    gtsam::PriorFactor<VALUE> priorFactor;

    double Psi(double alpha, double z) const {
      if (alpha==0) return -log(z)+z-1.0;
      else if (alpha<=-10.0) return z*log(z)+1.0;
      else{
        return (abs(alpha-2.0)/alpha)*((1.0-0.5*alpha)*pow(z,alpha/(alpha-2.0))+0.5*alpha*z-1.0);
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

