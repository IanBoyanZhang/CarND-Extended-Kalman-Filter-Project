#ifndef TOOLS_H_
#define TOOLS_H_
#define _USE_MATH_DEFINES
#include <vector>
#include <cmath>
#include "measurement_package.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

class Tools {
public:
  /**
  * Constructor.
  */
  Tools();

  /**
  * Destructor.
  */
  virtual ~Tools();

  /**
  * A helper method to calculate RMSE.
  */
  VectorXd CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth);

  /**
  * A helper method to calculate Jacobians.
  */
  MatrixXd CalculateJacobian(const VectorXd& x_state);

  /**
   * Coordinate transformation
   * @param x_state
   * @return
   */
  VectorXd Cart2Polar(const VectorXd& x_state);
  VectorXd Polar2Cart(const MeasurementPackage &measurement_pack);

  VectorXd NormalizeAngle(VectorXd &z_diff);
};

#endif /* TOOLS_H_ */
