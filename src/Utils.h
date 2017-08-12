//
// Created by x on 17-8-12.
//

#ifndef MPC_COORDUTILS_H
#define MPC_COORDUTILS_H

#include "Eigen/Core"

/**
 * "Normalize" the given angle into (-pi, pi]
 */
inline double normalize_angle(double angle) {
  while (angle <= -M_PI) angle += 2 * M_PI;
  while (angle >   M_PI) angle -= 2 * M_PI;
  return angle;
}

/**
 * Convert map coordinates to vehicle coordinates.
 *
 * Transformation from map coordinates to vehicle coordinates can be carried out by
 * applying affine transform. (In details, first a translation then a rotation)
 *
 * Supports batch transformation using Eigen to do the matrix stuff.
 *
 * All vectors shall be represented using homogeneous coordinates, with the last
 * components set as 1.
 *
 * To do batch transformation, vectors can be stored as a matrix.
 *
 * @param vehicle_x: the x position of the vehicle, in map coordinate
 * @param vehicle_y: the y position of the vehicle, in map coordinate
 * @param vehicle_orientation: the orientation of the vehicle, in map coordinate
 * @param points_in_map: points to be transformed, in map coordinate
 */
Eigen::MatrixXd mapCoordinates2VechileCoordinates(
    double vehicle_x, double vehicle_y, double vehicle_orientation,
    const Eigen::MatrixXd& points_in_map);

#endif //MPC_COORDUTILS_H
