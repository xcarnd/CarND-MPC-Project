//
// Created by x on 17-8-12.
//

#include "Utils.h"

Eigen::MatrixXd mapCoordinates2VehicleCoordinates(
    double vehicle_x, double vehicle_y, double vehicle_orientation,
    const Eigen::MatrixXd& points_in_map) {
  // denotes vx = vehicle_x, vy = vehicle_y, vpsi = vehicle_orientation for short
  //
  // said, vehicle at (vx, vy, vpsi) in map coordinate.
  //
  // to transform points in map coordinates to vehicle coordinate,
  // image we can fix the points, then transform can be archived by:
  // 1. move the rotated origin to (vx, vy)
  // 2. rotate the coordinate origin with (vpsi), which will make the x axis
  // pointing to the same direction as the x axis in vehicle coordinate.
  //
  // the same effect can be archived by making the coordinate fixed and move the points
  // in the opposite manner. that is, rotate the points with -vpsi, then translate by (-vx, -vy)
  //
  // prepare the transform matrix.
  // if c = cos(-vpsi), s = sin(-vpsi), then
  // rotation matrix R will be:
  // | c -s 0 |
  // | s  c 0 |
  // | 0  0 1 |
  // translation matrix T will be:
  // | 1 0 -vx |
  // | 0 1 -vy |
  // | 0 0   1 |
  // the final affine transform matrix is R * T, which can be pre-multiplied:
  // | c, -s, -c * vx + s * vy |
  // | s,  c, -s * vx - c * vy |
  // | 0,  0,                1 |
  const double p = -vehicle_orientation;
  const double c = std::cos(p);
  const double s = std::sin(p);
  Eigen::Matrix3d t;
  t <<   c, -s, -c * vehicle_x + s * vehicle_y,
         s,  c, -s * vehicle_x - c * vehicle_y,
         0,  0,                              1;
  return t * points_in_map;
}
