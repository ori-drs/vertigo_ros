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

          // Recover weight and set it up in outlierProcess factor
          double c = 1.0; // c is scaling param set before optimisation
          double w = weightAdaptive (weightedError, alpha.value(), c);
          weight_ = w;
          outlierProcess_->setWeight(w);

          // Jacobians
          // The Jacobians wrt the poses should be scaled
          if (H1) *H1 = *H1 * w;
          if (H2) *H2 = *H2 * w;
          // The Jacobian wrt alpha should be handled accordingly
          if (H3) *H3 = error * weightAdaptiveDerivativeAlpha(weightedError, alpha.value(), c);

          // Scale factor error by the weight from the adaptive kernel
          error *= w;

          // std::cout << "H1: " << H1 << std::endl;
          // std::cout << "H2: " << H2 << std::endl;
          // std::cout << "H3: " << H3 << std::endl;
          // std::cout << "error: " << error << std::endl;

          return error;
        }

      double getWeight() const {return weight_;}

    private:
      // Tolerance used for singularities
      double epsilon_ = 1E-5;
      
      // Use practical implementation flag
      bool usePractical_ = true;
      
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
        dw = weightAdaptivePracticalDerivativeAlpha(x, alpha, c);

        //  // Check derivative by simple finite differences
        // double w_m = weightAdaptivePractical(x, alpha - epsilon_, c); // psi(i-1)
        // double w = weightAdaptivePractical(x, alpha, c); // psi(i+1)
        // double dw_numeric = (w - w_m) / epsilon_;

        // std::cout << "dw: "<< dw << ", dw_numeric: " << dw_numeric << std::endl;

      } else{
        dw = weightAdaptiveAnalyticalDerivativeAlpha(x, alpha, c);
      }

      // std::cout << "weightAdaptiveDerivativeAlpha(x=" << x << ", alpha=" << alpha << ", c=" << c << ") : " << dw << std::endl;

      return dw;
    }

    // Analytical expressions
    double weightAdaptiveAnalytical(double x, double alpha, double c) const {
      double w;

      if (alpha >= 2){
        w = 1.0 / pow(c,2);

      } else if (abs(alpha) <= epsilon_){
        w = 2.0 / (pow(x, 2) + 2 * pow(c, 2));

      } else if (alpha <= ShapeParameter::MIN){
        w = 1 / pow(c, 2) * exp(-0.5 * pow(x/c, 2));

      } else{
        w = (1/pow(c,2)) * pow( pow(x/c, 2) / abs(alpha-2) + 1, (0.5*alpha-1));
      }

      return w;
    }

    double weightAdaptiveAnalyticalDerivativeAlpha(double x, double alpha, double c) const {
      double dw;

      if (alpha >= 2 || abs(alpha)<= epsilon_ || alpha <= ShapeParameter::MIN){
        dw = 0;
      }
      else{
        dw = pow(1 + pow(x, 2)/(pow(c, 2)*fabs(alpha - 2)), 0.5*alpha - 1)*(0.5*log(1 + pow(x, 2)/(pow(c, 2)*fabs(alpha - 2))) - pow(x, 2)*(0.5*alpha - 1)*(((alpha - 2) > 0) - ((alpha - 2) < 0))/(pow(c, 2)*(1 + pow(x, 2)/(pow(c, 2)*fabs(alpha - 2)))*pow(alpha - 2, 2)))/pow(c, 2);
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
