#include <iostream>
#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;
  is_setup = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}
/**
 *
 * long long or time_t?
 * Compute the time elapsed between the current and previous measurement
 * @param curr_time
 * @param prev_time
 * @return float_t dt
 */
float_t FusionEKF::GetTimeDiff(long long curr_time, long long prev_time) {
  return (float_t)(curr_time - prev_time)/1000000.0;
}

/**
 * Better estimation of ax and ay over time?
 * @param dt
 * @return MatrixXd * Q process noise
 */
MatrixXd FusionEKF::ConstructQ(float_t dt) {
  float_t dt_2 = pow(dt, 2);
  float_t dt_3 = pow(dt, 3);
  float_t dt_4 = pow(dt, 4);

  MatrixXd Q(4, 4);
  Q << dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
       0, dt_4/4 * noise_ay, 0, dt_3/2*noise_ay,
       dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
       0, dt_3/2*noise_ay, 0, dt_2*noise_ay;
  return Q;
}

bool isRadar(const MeasurementPackage &measurement_pack) {
  return measurement_pack.sensor_type_ == MeasurementPackage::RADAR;
}

bool isLaser(const MeasurementPackage &measurement_pack) {
  return measurement_pack.sensor_type_ == MeasurementPackage::LASER;
}

/**
 * TODO: Refactor code to separate LIDAR and RADAR data
 * @param measurement_pack
 */

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  VectorXd meas_cart;

  /**
   * TODO: Refactor this part into a function
   */
  if (!is_setup) {
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 0, 0, 0, 0;

    ekf_.F_ = MatrixXd(4, 4);
    ekf_.F_ << 1, 0, 1, 0,
               0, 1, 0, 1,
               0, 0, 1, 0,
               0, 0, 0, 1;
    // State covariance matrix
    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ << 1, 0, 0, 0,
              0, 1, 0, 0,
              0, 0, 1000, 0,
              0, 0, 0, 1000;

    is_setup = true;
    return;
  }

  if (!is_initialized_) {
    /**
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      meas_cart = tools.Polar2Cart(measurement_pack);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      meas_cart = measurement_pack.raw_measurements_;
    }
    ekf_.x_ << meas_cart[0], meas_cart[1], 0, 0;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }


 /*****************************************************************************
   *  Parameter preparation
   ****************************************************************************/
  float_t dt = GetTimeDiff(measurement_pack.timestamp_, previous_timestamp_);
  previous_timestamp_ = measurement_pack.timestamp_;

  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;
  ekf_.Q_ = ConstructQ(dt);

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
 if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.H_ = tools.Cart2Polar(ekf_.x_);
    ekf_.H_k_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;
  } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.H_k_ = H_laser_;
    ekf_.R_ = R_laser_;
  }

  float_t threshold = 1e-3;
  if (dt >= threshold) {
    ekf_.Predict();
  }
  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    meas_cart = measurement_pack.raw_measurements_;
    ekf_.UpdateEKF(meas_cart);
  } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
    // Laser updates
    meas_cart = measurement_pack.raw_measurements_;
    ekf_.Update(meas_cart);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
