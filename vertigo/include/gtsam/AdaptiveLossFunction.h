#ifndef ADAPTIVELOSSFUNCTION_H_
#define ADAPTIVELOSSFUNCTION_H_

#pragma once

//#include <Eigen/Dense>
#include <iostream>
#include <math.h>
#include "gtsam/linear/NoiseModel.h"

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

namespace vertigo {

/// Adaptive implements the "Adaptive" robust error model (Jonathan T. Barron)
class GTSAM_EXPORT Adaptive : public gtsam::noiseModel::mEstimator::Base {
protected:
  double c_, csquared_;
  double alpha_;

public:

  typedef boost::shared_ptr<Adaptive> shared_ptr;

  Adaptive (double c = 1.0, double alpha = 2.0, const ReweightScheme reweight = Block);
  double weight(double error) const {
    if (alpha_ == 2.0) return (error/csquared_);
    else if (alpha_ == 0.0) return (2*error/(error*error+2*csquared_));
    else if (alpha_ < -20)  return (error/csquared_)*exp(-0.5*(error*error/csquared_));
    else (error/csquared_)*pow(((error*error/csquared_)/abs(alpha_-2)+1),0.5*alpha_-1);
  }
  void print(const std::string &s) const;
  bool equals(const Base& expected, double tol=1e-8) const;
  static shared_ptr Create(double k, double alpha, const ReweightScheme reweight = Block) ;

private:
  /** Serialization function */
  friend class boost::serialization::access;
  template<class ARCHIVE>
  void serialize(ARCHIVE & ar, const unsigned int /*version*/) {
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Base);
    ar & BOOST_SERIALIZATION_NVP(c_);
  }
};

} // vertigo namespace


#endif /* ADAPTIVELOSSFUNCTION_H_ */
