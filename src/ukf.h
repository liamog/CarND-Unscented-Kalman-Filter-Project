#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

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
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(double delta_t, const VectorXd& measurement);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(double delta_t, const VectorXd& measurement);

  /**
   * Generates the sigma points.
   */
  void AugmentedSigmaPoints(MatrixXd* Xsig_out);

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
   * @param delta_t
   */
  void PredictRadarMeasurement(double delta_t);

  /**
   * Accessor for the state x.
   * @return state x vector.
   */
  const VectorXd &x() const {return x_;}

 private:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  // Measurement noise radar
  MatrixXd R_rad_;

  // Measurement noise lidar
  MatrixXd R_lid_;

  MatrixXd S_rad_;
  VectorXd z_rad_pred_;
  MatrixXd T_rad_;

  MatrixXd H_lid_;
  MatrixXd H_lid_t_;
  MatrixXd I_;

  ///* time when the state is true, in us
  long long time_us_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_;

  ///* Radar measurement space dimension
  int n_z_;

  ///* Augmented state dimension
  int n_aug_;

  ///* Sigma point spreading parameter
  double lambda_;

  // Number of sigma points
  int n_sigma_;

  // NIS samples for Radar
  std::vector<std::pair<double, double>> radar_nis_samples_;
  // NIS samples for Lidar
  std::vector<std::pair<double, double>> lidar_nis_samples_;

};

#endif /* UKF_H */
