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
      BetweenFactorAdaptive(gtsam::Key key1, gtsam::Key key2, gtsam::Key key3, const VALUE& measured, const gtsam::SharedNoiseModel& model)
      : gtsam::NoiseModelFactor3<VALUE, VALUE, ShapeParameter>(model, key1, key2, key3),
        betweenFactor(key1, key2, measured, model) {}

      gtsam::Vector evaluateError(const VALUE& p1, const VALUE& p2, const ShapeParameter& alpha,
                                  boost::optional<gtsam::Matrix&> H1 = boost::none,
                                  boost::optional<gtsam::Matrix&> H2 =  boost::none,
                                  boost::optional<gtsam::Matrix&> H3 =  boost::none) const
        {

          // calculate error
          gtsam::Vector error = betweenFactor.evaluateError(p1, p2, H1, H2);

          double error_dis = this->noiseModel_->distance(error);

          double c = 1; // c is scalling param set before optimisation
          double w = weight_adaptive (error_dis, alpha.value(), c);
          weight_ = w;
          alpha.setWeight_z(w);


          error *= w;

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


          // handle derivatives MR: Need to be modified according to adaptive kernel
          if (H1) *H1 = *H1 * w;
          if (H2) *H2 = *H2 * w;
//          if (H3) *H3 = error;

          if (alpha.value()==2 || alpha.value()==0 || alpha.value()<= -10){
            if (H3) *H3 = error;
          } else {
            if (H3) *H3 = error * (-w*(alpha.value()-1)*pow(pow(error_dis/c,2)/(alpha.value()-2)+1,-1)*
                                   pow(error_dis/c,2)/pow(alpha.value()-2,2));
          }

//          std::cout << "error is: \n" << error << std::endl;
          return error;
        }

      double getWeight() const {return weight_;}

    private:
      gtsam::BetweenFactor<VALUE> betweenFactor;
      mutable double weight_;

      double weight_adaptive(double x, double alpha, double c) const {
        if (alpha == 2) return 1/pow(c,2);
        else if (alpha == 0) return 2/(pow(x,2) + 2*pow(c,2));
        else if (alpha <= -10) return (1/pow(c,2)) * exp(-0.5*pow(x/c,2));
        else return (1/pow(c,2))*pow(pow(x/c,2)/abs(alpha-2)+1,0.5*alpha-1);
      }

  };

} // vertigo namespace

#endif /* BETWEENFACTORADAPTIVE_H_ */
