/*
 * shapeParamater.h
 *
 *  Created on: 26.05.2020
 *      Author: Milad Ramezani
 */

#ifndef SHAPEPARAMETER_H_
#define SHAPEPARAMETER_H_

#pragma once

#include <gtsam/base/Lie.h>

namespace vertigo {

  /**
   * ShapeParameter is a wrapper around double to allow it to be a Lie type
   */
  struct ShapeParameter {

    /** default constructor */
    ShapeParameter() : d_(2.0) {}

    /** wrap a double */
    ShapeParameter(double d) : d_(d) {
      if (d_ > 2.0) d_= 2.0; // maximum shape param is 2.0
      if (d_ < -10) d_= -10; // maximum shape param is 2.0
    }

    /** access the underlying value */
    double value() const { return d_; }

    /** print @param s optional string naming the object */
    inline void print(const std::string& name="") const {
      std::cout << name << ": " << d_ << std::endl;
    }

    /** equality up to tolerance */
    inline bool equals(const ShapeParameter& expected, double tol=1e-5) const {
      return fabs(expected.d_ - d_) <= tol;
    }

    // Manifold requirements

    /** Returns dimensionality of the tangent space */
    inline size_t dim() const { return 1; }
    inline static size_t Dim() { return 1; }

    /** Update the ShapeParameter with a tangent space update */
    inline ShapeParameter retract(const gtsam::Vector& v) const {
      double x = value() + v(0);

      if (x>2.0) x=2.0;
      else if (x<-10.0) x=-10.0;

      return ShapeParameter(x);
    }

    /** @return the local coordinates of another object */
    inline gtsam::Vector localCoordinates(const ShapeParameter& t2) const { return gtsam::Vector1(t2.value() - value()); }

    // Group requirements

    /** identity */
    inline static ShapeParameter identity() {
      return ShapeParameter();
    }

    /** compose with another object */
    inline ShapeParameter compose(const ShapeParameter& p) const {
      return ShapeParameter(d_ + p.d_);
    }

    /** between operation */
    inline ShapeParameter between(const ShapeParameter& l2,
        boost::optional<Matrix&> H1=boost::none,
        boost::optional<Matrix&> H2=boost::none) const {
      if(H1) *H1 = -eye(1);
      if(H2) *H2 = eye(1);
      return ShapeParameter(l2.value() - value());
    }

    /** invert the object and yield a new one */
    inline ShapeParameter inverse() const {
      return ShapeParameter(-1.0 * value());
    }

    // Lie functions

    /** Expmap around identity */
    static inline ShapeParameter Expmap(const gtsam::Vector& v) { return ShapeParameter(v(0)); }

    /** Logmap around identity - just returns with default cast back */
    static inline gtsam::Vector Logmap(const ShapeParameter& p) { return gtsam::Vector1(p.value()); }

  private:
      double d_;
  };
}

namespace gtsam {
// Define Key to be Testable by specializing gtsam::traits
template<typename T> struct traits;
template<> struct traits<vertigo::ShapeParameter> {
  static void Print(const vertigo::ShapeParameter& key, const std::string& str = "") {
    key.print(str);
  }
  static bool Equals(const vertigo::ShapeParameter& key1, const vertigo::ShapeParameter& key2, double tol = 1e-8) {
    return key1.equals(key2, tol);
  }
  static int GetDimension(const vertigo::ShapeParameter & key) {return key.Dim();}

  typedef OptionalJacobian<3, 3> ChartJacobian;
  typedef gtsam::Vector TangentVector;
  static TangentVector Local(const vertigo::ShapeParameter& origin, const vertigo::ShapeParameter& other,
  ChartJacobian Horigin = boost::none, ChartJacobian Hother = boost::none) {
    return origin.localCoordinates(other);
  }
  static vertigo::ShapeParameter Retract(const vertigo::ShapeParameter& g, const TangentVector& v,
        ChartJacobian H1 = boost::none, ChartJacobian H2 = boost::none) {
      return g.retract(v);
    }
};
}



#endif /* SHAPEPARAMETER_H_ */

