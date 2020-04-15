/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation, 
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 *  @file  fullSLAM.h
 *  @brief: measurements in 3D space
 *  @author Milad Ramezani
 */

#pragma once

#include <gtsam/slam/PriorFactor.h>
#include <gtsam/slam/BetweenFactor.h>
#include <gtsam/slam/RangeFactor.h>
#include <gtsam/slam/BearingFactor.h>
#include <gtsam/slam/BearingRangeFactor.h>
#include <gtsam/nonlinear/NonlinearEquality.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#include <gtsam/geometry/Pose2.h> // TODO: must be removed
#include <gtsam/geometry/Pose3.h>

// Use planarSLAM namespace for specific SLAM instance
namespace fullSLAM {

  using namespace gtsam;

  /// Convenience function for constructing a pose key
  inline Symbol PoseKey(Eigen::Index j) { return Symbol('x', j); }

  /// Convenience function for constructing a pose key
  inline Symbol PointKey(Eigen::Index j) { return Symbol('l', j); }

  /*
   * List of typedefs for factors
   */
  /// A hard constraint for PoseKeys to enforce particular values
  typedef NonlinearEquality<Pose3> Constraint;
  /// A prior factor to bias the value of a PoseKey
  typedef PriorFactor<Pose3> Prior;
  /// A factor between two PoseKeys set with a Pose2
  typedef BetweenFactor<Pose3> Odometry;
  /// A factor between a PoseKey and a PointKey to express difference in rotation (set with a Rot3)
  typedef BearingFactor<Pose3, Point3> Bearing;
  /// A factor between a PoseKey and a PointKey to express distance between them (set with a double)
  typedef RangeFactor<Pose3, Point3> Range;
  /// A factor between a PoseKey and a PointKey to express difference in rotation and location
  typedef BearingRangeFactor<Pose3, Point3> BearingRange;

  /**  Values class, using specific PoseKeys and PointKeys
   * Mainly as a convenience for MATLAB wrapper, which does not allow for identically named methods
   */
  struct Values: public gtsam::Values {

    /// Default constructor
    Values() {}

    /// Copy constructor
    Values(const gtsam::Values& values) :
      gtsam::Values(values) {
    }

    /// get a pose
    Pose3 pose(Eigen::Index key) const { return at<Pose3>(PoseKey(key)); }

    /// get a point
    Point3 point(Eigen::Index key) const { return at<Point3>(PointKey(key)); }

    /// insert a pose
    void insertPose(Eigen::Index key, const Pose3& pose) { insert(PoseKey(key), pose); }

    /// insert a point
    void insertPoint(Eigen::Index key, const Point3& point) { insert(PointKey(key), point); }
  };

  /// Creates a NonlinearFactorGraph with the Values type
  struct Graph: public NonlinearFactorGraph {

    /// Default constructor for a NonlinearFactorGraph
    Graph(){}

    /// Creates a NonlinearFactorGraph based on another NonlinearFactorGraph
    Graph(const NonlinearFactorGraph& graph);

    /// Biases the value of PoseKey key with Pose3 p given a noise model
    void addPrior(Eigen::Index poseKey, const Pose3& pose, const SharedNoiseModel& noiseModel);

    /// Creates a hard constraint to enforce Pose3 p for PoseKey poseKey's value
    void addPoseConstraint(Eigen::Index poseKey, const Pose3& pose);

    /// Creates a factor with a Pose3 between PoseKeys poseKey and pointKey (poseKey.e. an odometry measurement)
    void addOdometry(Eigen::Index poseKey1, Eigen::Index poseKey2, const Pose3& odometry, const SharedNoiseModel& model);

    /// Creates a factor with a Rot3 between a PoseKey poseKey and PointKey pointKey for difference in rotation
//    void addBearing(Eigen::Index poseKey, Eigen::Index pointKey, const Rot3& bearing, const SharedNoiseModel& model);

    /// Creates a factor with a Rot3 between a PoseKey poseKey and PointKey pointKey for difference in location
    void addRange(Eigen::Index poseKey, Eigen::Index pointKey, double range, const SharedNoiseModel& model);

    /// Creates a factor with a Rot3 between a PoseKey poseKey and PointKey pointKey for difference in rotation and location
//    void addBearingRange(Eigen::Index poseKey, Eigen::Index pointKey, const Rot3& bearing, double range, const SharedNoiseModel& model);

    /// Optimize
    Values optimize(const Values& initialEstimate) {
      return LevenbergMarquardtOptimizer(*this, initialEstimate).optimize();
    }
  };

} // fullSLAM


