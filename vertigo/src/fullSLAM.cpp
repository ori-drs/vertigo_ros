/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation, 
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 *  @file  fullSLAM.cpp
 *  @brief: measurements in 3D space
 *  @author Milad Ramezani
 **/

#include "fullSLAM.h"
#include <gtsam/nonlinear/NonlinearFactorGraph.h>

// Use planarSLAM namespace for specific SLAM instance
namespace fullSLAM {

  /* ************************************************************************* */
  Graph::Graph(const NonlinearFactorGraph& graph) :
      NonlinearFactorGraph(graph) {}

  /* ************************************************************************* */
  void Graph::addPrior(Eigen::Index i, const Pose3& p, const SharedNoiseModel& model) {
    sharedFactor factor(new Prior(PoseKey(i), p, model));
    push_back(factor);
  }

  /* ************************************************************************* */
  void Graph::addPoseConstraint(Eigen::Index i, const Pose3& p) {
    sharedFactor factor(new Constraint(PoseKey(i), p));
    push_back(factor);
  }

  /* ************************************************************************* */
  void Graph::addOdometry(Eigen::Index i, Eigen::Index j, const Pose3& odometry, const SharedNoiseModel& model) {
    sharedFactor factor(new Odometry(PoseKey(i), PoseKey(j), odometry, model));
    push_back(factor);
  }

  /* ************************************************************************* */
//  void Graph::addBearing(Eigen::Index i, Eigen::Index j, const Rot2& z, const SharedNoiseModel& model) {
//    sharedFactor factor(new Bearing(PoseKey(i), PointKey(j), z, model));
//    push_back(factor);
//  }

  /* ************************************************************************* */
  void Graph::addRange(Eigen::Index i, Eigen::Index j, double z, const SharedNoiseModel& model) {
    sharedFactor factor(new Range(PoseKey(i), PointKey(j), z, model));
    push_back(factor);
  }

  /* ************************************************************************* */
//  void Graph::addBearingRange(Eigen::Index i, Eigen::Index j, const Rot2& z1,
//      double z2, const SharedNoiseModel& model) {
//    sharedFactor factor(new BearingRange(PoseKey(i), PointKey(j), z1, z2, model));
//    push_back(factor);
//  }

  /* ************************************************************************* */

} // fullSLAM

