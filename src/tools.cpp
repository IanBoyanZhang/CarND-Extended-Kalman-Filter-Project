#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
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
  MatrixXd Hj = MatrixXd::Zero(3, 4);

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

/**
 * @param x_state
 * @return
 */
VectorXd Tools::Cart2Polar(const VectorXd &x_state) {
  VectorXd H_x(3);

  float_t px = x_state(0);
  float_t py = x_state(1);
  float_t vx = x_state(2);
  float_t vy = x_state(3);

  float_t ro = sqrt(pow(px, 2) + pow(py, 2));
  float_t theta = atan2f(py, px);

  float_t threshold = 1e-4;
  float_t ro_dot = 0;
  if (fabs(ro) >= threshold) {
    ro_dot = (px * vx + py * vy)/ro;
  } else {
    cout << "Cart2Polar() - Error - Division by Zero" << endl;
  }

  H_x << ro, theta, ro_dot;
  return H_x;
}

/**
 *
 * @param measurement_pack
 * @return
 */
VectorXd Tools::Polar2Cart(const MeasurementPackage &measurement_pack) {
  float_t ro = measurement_pack.raw_measurements_[0];
  float_t theta = measurement_pack.raw_measurements_[1];
  float_t ro_dot = measurement_pack.raw_measurements_[2];

  float_t px = ro*cos(theta);
  float_t py = ro*sin(theta);

  VectorXd radar_cart_x(4);
  radar_cart_x << px, py, 0, 0;
  return radar_cart_x;
}

/**
 *
 * @param theta_
 * @return
 */
VectorXd Tools::NormalizeAngle(VectorXd &z_diff) {
  float_t TWO_PI = 2. * M_PI;
  while(z_diff(1) > M_PI) { z_diff(1) -= TWO_PI; }
  while(z_diff(1) < -M_PI) { z_diff(1) += TWO_PI; }
  return z_diff;
}

/**
 * Alternative solution
 */
/*
VectorXd Tools::NormalizeAngle(VectorXd &z_diff) {
  return fmodf(z_diff, 2. * M_PI);
}
 */
