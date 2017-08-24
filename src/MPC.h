#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"

using namespace std;

// Define N and dt
constexpr int N = 6;
// dt shall be at least 0.1 second since there are 100ms latency.
// 0.05 for the interval between command is issued and MPC received sensor updates.
constexpr double dt = 0.1 + 0.1;

class MPC {
 public:
  MPC();

  virtual ~MPC();

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuatotions.
  vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);
};

#endif /* MPC_H */
