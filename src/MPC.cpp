#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"
#include "Utils.h"

using CppAD::AD;

constexpr std::size_t x_start = 0;
constexpr std::size_t y_start = x_start + N;
constexpr std::size_t psi_start = y_start + N;
constexpr std::size_t v_start = psi_start + N;
constexpr std::size_t cte_start = v_start + N;
constexpr std::size_t epsi_start = cte_start + N;
constexpr std::size_t delta_start = epsi_start + N;
constexpr std::size_t a_start = delta_start + (N - 1);

namespace {
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
}

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
constexpr double Lf = 2.67;

double v_ref = mph_to_mps(15.0); 

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  Eigen::VectorXd coeffs_d;
  FG_eval(Eigen::VectorXd coeffs) {
    this->coeffs = coeffs;
    // if y = a_n * x^n + a_(n-1) * x^(n-1) + ... + a_1 * x + a_0
    // then y' = a_n * n * x^(n-1) + a_(n-1) * (n-1) * x^(n-2) + ... + a_1
    std::vector<double> d;
    std::size_t highest_order = coeffs.size() - 1;
    for (std::size_t i = 1; i <= highest_order; ++i) {
      d.push_back(coeffs(i) * i);
    }
    this->coeffs_d = Eigen::Map<Eigen::VectorXd>(d.data(), d.size());
    //    std::cout<<"coeffs"<<this->coeffs<<std::endl;
    //    std::cout<<"coeffs_d"<<this->coeffs_d<<std::endl;
  }

  
  AD<double> polyeval(const Eigen::VectorXd& coeffs, AD<double> v) {
    AD<double> result = 0;
    for (int i = 0; i < coeffs.size(); ++i) {
      result += coeffs(i) * CppAD::pow(v, i);
    }
    return result;
  }


  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // TODO: implement MPC
    // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.

    // fg[0] stores the cost
    fg[0] = 0;
    // minimize cte and epsi and difference with the target v
    for (std::size_t i = 0; i < N; ++i) {
      fg[0] += CppAD::pow(vars[cte_start + i], 2);
      fg[0] += CppAD::pow(vars[epsi_start + i], 2);
      fg[0] += CppAD::pow(vars[v_start + i] - v_ref, 2);
    }

    // minimize actuators used
    for (std::size_t i = 0; i < N - 1; ++i) {
      fg[0] += 100 * CppAD::pow(vars[delta_start + i], 2);
      fg[0] += 10 * CppAD::pow(vars[a_start + i], 2);
    }

    // minimize the gap between two consecutive inputs
    for (std::size_t i = 0; i < N - 2; ++i) {
      fg[0] += 500 * CppAD::pow(vars[delta_start + i + 1] - vars[delta_start + i], 2);
      fg[0] += 50 * CppAD::pow(vars[a_start + i + 1] - vars[a_start + i], 2);
    }

    fg[1 + x_start] = vars[x_start];
    fg[1 + y_start] = vars[y_start];
    fg[1 + psi_start] = vars[psi_start];
    fg[1 + v_start] = vars[v_start];
    fg[1 + cte_start] = vars[cte_start];
    fg[1 + epsi_start] = vars[epsi_start];

    // motion model
    for (std::size_t i = 1; i < N; ++i) {
      AD<double> x0 = vars[x_start + i - 1];
      AD<double> y0 = vars[y_start + i - 1];
      AD<double> psi0 = vars[psi_start + i - 1];
      AD<double> v0 = vars[v_start + i - 1];
      AD<double> cte0 = vars[cte_start + i - 1];
      AD<double> epsi0 = vars[epsi_start + i - 1];

      AD<double> x1 = vars[x_start + i];
      AD<double> y1 = vars[y_start + i];
      AD<double> psi1 = vars[psi_start + i];
      AD<double> v1 = vars[v_start + i];
      AD<double> cte1 = vars[cte_start + i];
      AD<double> epsi1 = vars[epsi_start + i];

      AD<double> delta0 = vars[delta_start + i - 1];
      AD<double> a0 = vars[a_start + i - 1];

      AD<double> y_ref = this->polyeval(this->coeffs, x0);
      AD<double> d_y_ref = this->polyeval(this->coeffs_d, x0);
      AD<double> psi_ref = CppAD::atan(d_y_ref);

      // x1 = x0 + v0 * dt * cos(psi0)
      fg[1 + x_start + i] = x0 + v0 * dt * CppAD::cos(psi0) - x1;
      // y1 = y0 + v0 * dt * cos(psi0)
      fg[1 + y_start + i] = y0 + v0 * dt * CppAD::sin(psi0) - y1;
      // psi1 = psi0 + v0 / Lf * delta0 * dt
      fg[1 + psi_start + i] = psi0 + v0 / Lf * delta0 * dt - psi1;
      // v1 = v0 + a0 * dt
      fg[1 + v_start + i] = v0 + a0 * dt - v1;
      // cte1 = cte0 + v0 * dt * sin(epsi0)
      //      = (y0 - y_ref) + v0 * dt * sin(epsi0)
      fg[1 + cte_start + i] = (y0 - y_ref) + v0 * dt * CppAD::sin(epsi0) - cte1;
      // epsi1 = epsi0 + v0 / Lf * delta0 * dt
      //       = (psi0 - psi_ref) + v0 / Lf * delta0 * dt
      fg[1 + epsi_start + i] = (psi0 - psi_ref) + v0 / Lf * delta0 * dt - epsi1;
    }
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  double x = state(0);
  double y = state(1);
  double psi = state(2);
  double v = state(3);
  double cte = state(4);
  double epsi = state(5);

  std::cout<<x<<" "<<y<<" "<<psi<<" "<<v<<" "<<cte<<" "<<epsi<<std::endl;

  // TODO: Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  //
  // 4 * 10 + 2 * 9
  size_t n_vars = 6 * N + 2 * (N - 1);
  // TODO: Set the number of constraints
  size_t n_constraints = 6 * N;
  

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (std::size_t i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }
  vars[x_start] = x;
  vars[y_start] = y;
  vars[psi_start] = psi;
  vars[v_start] = v;
  vars[cte_start] = cte;
  vars[epsi_start] = epsi;

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  // unlimited variables except the inputs
  for (std::size_t i = 0; i < delta_start; ++i) {
    vars_lowerbound[i] = -1e6;
    vars_upperbound[i] = 1e6;
  }
  for (std::size_t i = delta_start; i < a_start; ++i) {
    vars_lowerbound[i] = -1 * deg2rad(25);
    vars_upperbound[i] = 1 * deg2rad(25);
  }
  for (std::size_t i = a_start; i < n_vars; ++i) {
    vars_lowerbound[i] = -1;
    vars_upperbound[i] = 1;
  }

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (std::size_t i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }
  constraints_lowerbound[x_start] = x;
  constraints_lowerbound[y_start] = y;
  constraints_lowerbound[psi_start] = psi;
  constraints_lowerbound[v_start] = v;
  constraints_lowerbound[cte_start] = cte;
  constraints_lowerbound[epsi_start] = epsi;

  constraints_upperbound[x_start] = x;
  constraints_upperbound[y_start] = y;
  constraints_upperbound[psi_start] = psi;
  constraints_upperbound[v_start] = v;
  constraints_upperbound[cte_start] = cte;
  constraints_upperbound[epsi_start] = epsi;

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;
  std::cout << "Optimizing done? " << ok << std::endl;
  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  // TODO: Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.

  std::vector<double> x_vals;
  std::vector<double> y_vals;
  std::vector<double> psi_vals;
  std::vector<double> v_vals;
  std::vector<double> cte_vals;
  std::vector<double> epsi_vals;
  std::vector<double> delta_vals;
  std::vector<double> a_vals;
  
  std::cout<<"x = [";
  for (std::size_t t = 0; t < N; ++t) {
    std::cout<<solution.x[x_start + t]<<",";
    x_vals.push_back(solution.x[x_start + t]);
  }
  std::cout<<"]"<<std::endl;

  std::cout<<"y = [";
  for (std::size_t t = 0; t < N; ++t) {
    std::cout<<solution.x[y_start + t]<<",";
    y_vals.push_back(solution.x[y_start + t]);
  }
  std::cout<<"]"<<std::endl;

  std::cout<<"psi = [";
  for (std::size_t t = 0; t < N; ++t) {
    std::cout<<solution.x[psi_start + t]<<",";
    psi_vals.push_back(solution.x[psi_start + t]);
  }
  std::cout<<"]"<<std::endl;

  std::cout<<"v = [";
  for (std::size_t t = 0; t < N; ++t) {
    std::cout<<solution.x[v_start + t]<<",";
    v_vals.push_back(solution.x[v_start + t]);
  }
  std::cout<<"]"<<std::endl;
  
  std::cout<<"cte = [";
  for (std::size_t t = 0; t < N; ++t) {
    std::cout<<solution.x[cte_start + t]<<",";
    cte_vals.push_back(solution.x[cte_start + t]);
  }
  std::cout<<"]"<<std::endl;
  
  std::cout<<"epsi = [";
  for (std::size_t t = 0; t < N; ++t) {
    std::cout<<solution.x[epsi_start + t]<<",";
    epsi_vals.push_back(solution.x[epsi_start + t]);
  }
  std::cout<<"]"<<std::endl;
  
  std::cout<<"delta = [";
  for (std::size_t t = 0; t < N - 1; ++t) {
    std::cout<<solution.x[delta_start + t]<<",";
    delta_vals.push_back(solution.x[delta_start + t]);
  }
  std::cout<<"]"<<std::endl;

  std::cout<<"a = [";
  for (std::size_t t = 0; t < N - 1; ++t) {
    std::cout<<solution.x[a_start + t]<<",";
    a_vals.push_back(solution.x[a_start + t]);
  }
  std::cout<<"]"<<std::endl;

  std::cout<<std::endl;

  std::cout<<solution.x[delta_start]<<"  --  "<<solution.x[a_start]<<std::endl;
  
  std::vector<double> result = {solution.x[delta_start], solution.x[a_start]};
  for (std::size_t t = 0; t < N; ++t) {
    result.push_back(solution.x[x_start + t]);
  }
  for (std::size_t t = 0; t < N; ++t) {
    result.push_back(solution.x[y_start + t]);
  }
  return result;
}
