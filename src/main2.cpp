#include <math.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "MPC.h"
#include "json.hpp"
#include "Utils.h"
#include <fstream>

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.rfind("}]");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

// Evaluate a polynomial.
double polyeval(const Eigen::VectorXd& coeffs, double x) {
  double result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}

// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(const Eigen::VectorXd& xvals, const Eigen::VectorXd& yvals,
                        int order) {
  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);
  Eigen::MatrixXd A(xvals.size(), order + 1);

  for (int i = 0; i < xvals.size(); i++) {
    A(i, 0) = 1.0;
  }

  for (int j = 0; j < xvals.size(); j++) {
    for (int i = 0; i < order; i++) {
      A(j, i + 1) = A(j, i) * xvals(j);
    }
  }

  auto Q = A.householderQr();
  auto result = Q.solve(yvals);
  return result;
}

int main() {
  std::ifstream infile("data.txt");
  vector<string> alldata;
  std::string line;
  while (std::getline(infile, line)) {
    alldata.push_back(line);
  }
  infile.close();
  // MPC is initialized here!
  MPC mpc;

    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
  
  for (size_t i = 0; i < alldata.size(); ++i) {
  string sdata = alldata[i];
    cout << sdata << endl;
    if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
      string s = hasData(sdata);
      if (s != "") {
        auto j = json::parse(s);
        string event = j[0].get<string>();
        if (event == "telemetry") {
          // j[1] is the data JSON object
          vector<double> ptsx = j[1]["ptsx"];
          vector<double> ptsy = j[1]["ptsy"];
          double px = j[1]["x"];
          double py = j[1]["y"];
          double psi = j[1]["psi"];
          double v = j[1]["speed"];
          v = mph2mps(v);

          // fitting polyline for the given waypoints
          Eigen::Map<Eigen::VectorXd> xvals(ptsx.data(), ptsx.size());
          Eigen::Map<Eigen::VectorXd> yvals(ptsy.data(), ptsy.size());
          // fitting x = f(y) OR y = f(x)
          // for functions, domain values must be always increasing or decreasing,
          // which is not always true for fitting eithe one in this project
          int fit_type = IsMonotone(ptsx) ? TYPE_X_MAPS_TO_Y : TYPE_Y_MAPS_TO_X;

          Eigen::VectorXd coeffs;
          double cte;
          double epsi;
          double d_sign;
          if (fit_type == TYPE_X_MAPS_TO_Y) {
            coeffs = polyfit(xvals, yvals, 3);

            Eigen::VectorXd coeffs_diff(3);
            coeffs_diff << 3 * coeffs(3),
                           2 * coeffs(2),
                               coeffs(1);

            // cte is calculated as (actual y - reference y)
            cte = py - polyeval(coeffs, px);
          
            // cte_psi is calculated as (actual psi - reference psi)
            // why not using psi_ref = atan(f'(x))?
            // because the result of atan is in (-0.5*pi, 0.5pi],
            // without further information it can not be turned back to [0, 2*pi)
            double sign_x = signnum(xvals[1] - xvals[0]);
            d_sign = sign_x;
            double x_0 = px;
            double x_1 = px + 0.01 * sign_x;
            double y_0 = polyeval(coeffs, x_0);
            double y_1 = polyeval(coeffs, x_1);
            double psi_ref = std::atan2(y_1 - y_0, x_1 - x_0);
            psi_ref = normalize_angle2(psi_ref);
            
            epsi = psi - psi_ref;
          } else {
            coeffs = polyfit(yvals, xvals, 3);

            Eigen::VectorXd coeffs_diff(3);
            coeffs_diff << 3 * coeffs(3),
                           2 * coeffs(2),
                               coeffs(1);

            // cte is calculated as (actual x - reference x)
            cte = px - polyeval(coeffs, py);

            // cte_psi is calculated as (actual psi - reference psi)
            double sign_y = signnum(yvals[1] - yvals[0]);
            d_sign = sign_y;
            double y_0 = py;
            double y_1 = y_0 + 0.01 * sign_y;
            double x_0 = polyeval(coeffs, y_0);
            double x_1 = polyeval(coeffs, y_1);
            double psi_ref = std::atan2(y_1 - y_0, x_1 - x_0);
            psi_ref = normalize_angle2(psi_ref);

            epsi = psi - psi_ref;
          }

          /*
          * TODO: Calculate steering angle and throttle using MPC.
          *
          * Both are in between [-1, 1].
          *
          */
          cte = cte;
          epsi = epsi;
          Eigen::VectorXd state(6);
          state << px, py, psi, v, cte, epsi;
          const vector<double> result = mpc.Solve(state, coeffs, fit_type, d_sign);
          double steer_value = result[0];
          double throttle_value = result[1];

          json msgJson;
          // NOTE: Remember to divide by deg2rad(25) before you send the steering value back.
          // Otherwise the values will be in between [-deg2rad(25), deg2rad(25] instead of [-1, 1].
          msgJson["steering_angle"] = -steer_value / deg2rad(25);
          msgJson["throttle"] = throttle_value;

          //Display the MPC predicted trajectory 
          vector<double> mpc_x_vals;
          vector<double> mpc_y_vals;

          Eigen::MatrixXd mpc_wps(3, N-1);
          for (std::size_t t = 1; t < N; ++t) {
            mpc_wps.col(t-1) <<
                result[2 + 2 * t],
                result[2 + 2 * t + 1],
                1;
          }
          std::cout<<"trjx = [";
          for (std::size_t t = 1; t < N; ++t) {
            std::cout<<result[2 + 2 * t]<<",";
          }
          std::cout<<"]"<<std::endl;;
          std::cout<<"trjy = [";
          for (std::size_t t = 1; t < N; ++t) {
            std::cout<<result[2 + 2 * t + 1]<<",";
          }
          std::cout<<"]"<<std::endl;;
          Eigen::MatrixXd trans_mpc_wps = mapCoordinates2VechileCoordinates(px, py, psi, mpc_wps);
          for (std::size_t i = 0; i < N - 1; ++i) {
            mpc_x_vals.push_back(trans_mpc_wps(0, i));
            mpc_y_vals.push_back(trans_mpc_wps(1, i));
          }

          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Green line

          msgJson["mpc_x"] = mpc_x_vals;
          msgJson["mpc_y"] = mpc_y_vals;

          //Display the waypoints/reference line
          vector<double> next_x_vals;
          vector<double> next_y_vals;

          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Yellow line
          Eigen::MatrixXd waypoints(3, N);
            

          if (fit_type == TYPE_X_MAPS_TO_Y) {
            double wp_start = ptsx[0];
            double wp_length = ptsx[ptsx.size() - 1] - wp_start;
            double step = wp_length / (N - 1);
            for (std::size_t t = 0; t < N; ++t) {
              double x = wp_start + t * step;
              waypoints.col(t) << x, 
                                  polyeval(coeffs, x),
                                  1;
            }
          } else {
            double wp_start = ptsy[0];
            double wp_length = ptsy[ptsy.size() - 1] - wp_start;
            double step = wp_length / (N - 1);
            for (std::size_t t = 0; t < N; ++t) {
              double y = wp_start + t * step;
              waypoints.col(t) << polyeval(coeffs, y),
                                  y,
                                  1;
            }
          }
          std::cout<<"refx = [";
          for (std::size_t t = 0; t < N; ++t) {
            std::cout<<waypoints.col(t)(0)<<",";
          }
          std::cout<<"]"<<std::endl;;
          std::cout<<"refy = [";
          for (std::size_t t = 0; t < N; ++t) {
            std::cout<<waypoints.col(t)(1)<<",";
          }
          std::cout<<"]"<<std::endl;


          Eigen::MatrixXd trans_waypoints = mapCoordinates2VechileCoordinates(px, py, psi, waypoints);
          for (std::size_t i = 0; i < N; ++i) {
            next_x_vals.push_back(trans_waypoints(0, i));
            next_y_vals.push_back(trans_waypoints(1, i));
          }


          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;


          auto msg = "42[\"steer\"," + msgJson.dump() + "]";
          std::cout << msg << std::endl;
          // Latency
          // The purpose is to mimic real driving conditions where
          // the car does actuate the commands instantly.
          //
          // Feel free to play around with this value but should be to drive
          // around the track with 100ms latency.
          //
          // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE
          // SUBMITTING.
        }
      } else {
        // Manual driving
      }
    }
  }
}
