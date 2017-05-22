#include <iostream>
#include "tools.h"
#include <cmath>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  // Check the validity of the following inputs:
  // * the estimation vector size should not be zero
  // * the estimation vector size should equal ground truth vector size

  if (estimations.size() != ground_truth.size()
          || estimations.size() == 0) {
    cout << "Invalid estimation or ground_truth data" << endl;
    return rmse;
  }

  // accumulate squared residuals
  for (unsigned int i=0; i < estimations.size(); i+=1) {
    VectorXd residual = estimations[i] - ground_truth[i];

    // Coefficient-wise multiplication
    residual = pow(residual.array(), 2);
    rmse += residual;
  }

  // Calculate the mean
  rmse = rmse/estimations.size();

  // Calculate the squared root
  rmse = rmse.array().sqrt();

  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
    * Calculate a Jacobian here.
  */
  MatrixXd Hj(3, 4);

  float_t threshold = 1e-4;
//  Recover state parameters
  float_t px = x_state(0);
  float_t py = x_state(1);
  float_t vx = x_state(2);
  float_t vy = x_state(3);

//  Pre-compute a set of terms to avoid repeated calculation
  float_t c1 = pow(px, 2) + pow(py, 2);
  float_t c2 = sqrt(c1);
  float_t c3 = c1 * c2;

//  Check division by zero
  if(fabs(c1) < threshold) {
    cout << "CalculateJacobian() - Error - Division by Zero" << endl;
    return Hj;
  }

//  Computer the Jacobian Matrix
  Hj << px/2, py/2, 0, 0,
        -py/c1, px/c1, 0, 0,
        py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;
  return Hj;
}
