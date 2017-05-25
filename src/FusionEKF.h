#ifndef FusionEKF_H_
#define FusionEKF_H_
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include "measurement_package.h"
#include "Eigen/Dense"
#include "kalman_filter.h"
#include "tools.h"

class FusionEKF {
public:
  /**
  * Constructor.
  */
  FusionEKF();

  /**
  * Destructor.
  */
  virtual ~FusionEKF();

  /**
  * Run the whole flow of the Kalman Filter from here.
  */
  void ProcessMeasurement(const MeasurementPackage &measurement_pack);

  /**
  * Kalman Filter update and prediction math lives in here.
  */
  KalmanFilter ekf_;


private:
  // check whether the tracking toolbox was initialized or not (first measurement)
  bool is_initialized_;

  float_t noise_ax_ = 9;
  float_t noise_ay_ = 9;

  // previous timestamp
  long long previous_timestamp_;

  // tool object used to compute Jacobian and RMSE
  Tools tools;
  Eigen::MatrixXd R_laser_;
  Eigen::MatrixXd R_radar_;
  Eigen::MatrixXd H_laser_;
  Eigen::MatrixXd H_laser_t;
//  Eigen::MatrixXd Hj_;

  float_t GetTimeDiff(long long curr_time, long long prev_time);
  Eigen::MatrixXd ConstructQ(float_t dt);

  MeasurementPackage::SensorType GetSensorType(const MeasurementPackage
                                               &measurement_pack);
};

#endif /* FusionEKF_H_ */
