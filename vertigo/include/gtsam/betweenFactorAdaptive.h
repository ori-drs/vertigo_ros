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

          // calculate error
          gtsam::Vector error = betweenFactor.evaluateError(p1, p2, H1, H2);

//          double error_dis = this->noiseModel_->distance(error);
          double error_dis = this->noiseModel_->distance(error);

          double c = 1.0; // c is scalling param set before optimisation
          double w = weight_adaptive (error_dis, alpha.value(), c);
          weight_ = w;
          outlierProcess_->setWeight(w);
//          alpha.setWeight_z(w);
//          std::cout << "[BetweenFactorAdaptive] set weight to shape parameter: " << std::to_string(weight_) << std::endl;


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
          if (abs(alpha.value()-2.0) <= 0.01 || abs(alpha.value()) <= 0.01 || alpha.value() <= -10.0){

            if (H3) *H3 = error;

//              cout << "H3 is:\n " << H3 << endl;



//             std::cout << "alpha.value() is:" <<alpha.value() << std::endl;
            // std::cout << "H3 is:\n" <<H3 << std::endl;
          } else {
//             cout << "alpha is: " << alpha.value() << endl;
//             if (H3) *H3 = error * weight_adaptive_alpha(error_dis, alpha.value(), c);

            if (H3) *H3 = error * (-w*pow(1-pow(error_dis/c,2)/(alpha.value()-2),-1)*
                                   error_dis/(pow(c,2)*abs(alpha.value()-2)));


//             if (H3) *H3 = error * w*(0.5*(log(pow(error_dis/c,2)/abs(alpha.value()-2)+1))-(0.5*alpha.value()-1)*pow(pow(error_dis/c,2)/abs(alpha.value()-2)+1,-1)*
//             pow(error_dis/c,2)/((alpha.value()-2)*abs(alpha.value()-2)));
//            cout << "weight_adaptive_alpha is: " << weight_adaptive_alpha(error_dis, alpha.value(), c) << endl;
          }


//          std::cout << "error is: \n" << error << std::endl;
          return error;
        }

      double getWeight() const {return weight_;}

    private:
      gtsam::BetweenFactor<VALUE> betweenFactor;
      mutable double weight_;
      OutlierProcess<ShapeParameter>* outlierProcess_;

      double epsilon = 1E-5;

      double weight_adaptive(double x, double alpha, double c) const {
        // practical
//        double b, d;
//        b = abs(alpha-2)+epsilon;
//        if (alpha>=0) d = alpha+epsilon;
//        if (alpha<0)  d = alpha-epsilon;

//        return (1/pow(c,2))*pow(pow(x/c,2)/b+1,(0.5*d-1));

        // analytical
        if (abs(alpha-2.0) <= 0.01) return 1/pow(c,2);
//        else if (alpha == 0.0) return 2/(pow(x,2)+2*pow(c,2));
//        else if (alpha <= -10) return 1/pow(c,2) * exp(-0.5*pow(x/c,2));
        else return (1/pow(c,2))*pow(pow(x/c,2)/abs(alpha-2)+1,(alpha/2-1));

      }

      double weight_adaptive_alpha(double x, double alpha, double c) const {
        // practical
//        double b, d;
//        b = abs(alpha-2)+epsilon;
//        if (alpha>=0) d = alpha+epsilon;
//        if (alpha<0)  d = alpha-epsilon;

//        double w_b, w_d;
//        w_b = -((d/2-1)*pow(x,2)*pow(pow(x,2)/(pow(c,2)*b)+1,(d/2-2)))/(pow(c,4)*pow(b,2))*(alpha-2)/abs(alpha-2);
//        w_d = (pow(pow(x,2)/(b*pow(c,2))+1,(d/2-1))*log(pow(x,2)/(b*pow(c,2))+1))/(2*pow(c,2));
//        return (w_b + w_d);

        // analytical
        if (abs(alpha-2.0) <= 0.01) return 0;
        else return (1/pow(c,2))*pow(pow(x,2)/(pow(c,2)*abs(alpha-2))+1,alpha/2-1)*(log(pow(x,2)/(pow(c,2)*abs(alpha-2))+1))/2-
                    (pow(x,2)*(alpha/2-1))/(pow(c,2)*(pow(x,2)/(pow(c,2)*abs(alpha-2))+1)*abs(alpha-2)*(alpha-2));

      }

  };

} // vertigo namespace

#endif /* BETWEENFACTORADAPTIVE_H_ */
