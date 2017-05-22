#ifndef FusionEKF_H_
#define FusionEKF_H_

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
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
  bool pos_is_initialized_;
  bool velo_is_initialized_;
  bool is_setup;

  // previous timestamp
  long long previous_timestamp_;

  // tool object used to compute Jacobian and RMSE
  Tools tools;
  Eigen::MatrixXd R_laser_;
  Eigen::MatrixXd R_radar_;
  Eigen::MatrixXd H_laser_;
  Eigen::MatrixXd Hj_;

  Eigen::MatrixXd R_;
  Eigen::MatrixXd H_;

  float_t GetTimeDiff(long long curr_time, long long prev_time);
  Eigen::MatrixXd ConstructQ(float_t dt);
  Eigen::VectorXd Polar2Cart(const MeasurementPackage &measurement_pack);

  bool isRadar(const MeasurementPackage &measurement_pack);
  bool isLaser(const MeasurementPackage &measurement_pack);
};

#endif /* FusionEKF_H_ */
