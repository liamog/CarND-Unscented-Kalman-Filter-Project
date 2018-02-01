#ifndef _USERS_LIAM_DEV_UDACITY_CARND_UNSCENTED_KALMAN_FILTER_PROJECT_SRC_UKF_H
#define _USERS_LIAM_DEV_UDACITY_CARND_UNSCENTED_KALMAN_FILTER_PROJECT_SRC_UKF_H

#include "Eigen/Dense"
#include "measurement_package.h"

#include <fstream>
#include <string>
#include <vector>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
 public:
  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix.
   * @param delta_t Time between k and k+1 in seconds.
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(const VectorXd& measurement);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(const VectorXd& measurement);

  /**
   * Generates the sigma points.
   */
  void GenerateAugmentedSigmaPoints(MatrixXd* Xsig_out);

  /**
   * Generates the sigam points and then predicts the values for Xsig_pred_
   * for the sigma points for delta_t.
   * @param delta_t
   */

  void PredictSigmaPoints(double delta_t);
  /**
   * Calculates the mean and covariance for the predicted sigma points.
   */
  void PredictMeanAndCovarianceSigmaPoints();

  /**
   * Update the predicted state in radar measurement space.
   */
  void PredictRadarMeasurement();

  /**
   * Accessor for the state x.
   * @return state x vector.
   */
  const VectorXd& x() const { return x_; }

 private:
  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  ///* Measurement noise radar
  MatrixXd R_rad_;

  ///* Measurement noise lidar
  MatrixXd R_lid_;

  ///* Process noise matrix.
  MatrixXd Q_;

  MatrixXd S_rad_;
  VectorXd z_rad_pred_;
  MatrixXd T_rad_;

  MatrixXd H_lid_;
  MatrixXd H_lid_t_;
  MatrixXd I_;

  ///* time when the state is true, in us
  long long time_us_;

  double current_timestamp_;

  ///* Weights of sigma points
  VectorXd weights_;

  // NIS samples for Radar
  std::fstream radar_nis_samples_;
  // NIS samples for Lidar
  std::fstream lidar_nis_samples_;

  int counter_ = 0;
};

#endif //_USERS_LIAM_DEV_UDACITY_CARND_UNSCENTED_KALMAN_FILTER_PROJECT_SRC_UKF_H
