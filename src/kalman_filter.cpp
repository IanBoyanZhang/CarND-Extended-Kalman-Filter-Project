#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
    * predict the state
  */
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
    * update the state by using Kalman Filter equations
  */

  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;

  y = tools.NormalizeAngle(y);

  MatrixXd Ht = H_k_.transpose();
  MatrixXd S = H_k_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;

  MatrixXd K = PHt * Si;

  // New estimate
  x_ = x_ + (K * y);
  P_ = (I - K * H_k_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
    * update the state by using Extended Kalman Filter equations
  */
  VectorXd y = z - H_;

  y = tools.NormalizeAngle(y);

  MatrixXd Ht = H_k_.transpose();
  cout << "P" << P_ << endl;
  MatrixXd S = H_k_ * P_ * Ht + R_;
//  cout << "R_" << H_k_ * P_ << endl;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;

  MatrixXd K = PHt * Si;

  // New estimate
  x_ = x_ + (K * y);
  P_ = (I - K * H_k_) * P_;
}
