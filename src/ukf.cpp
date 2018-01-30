#include "ukf.h"
#include "Eigen/Dense"
#include "tools.h"

#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

#define PRINT(x) Tools::Print(x, #x);

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  n_x_ = 5;
  n_aug_ = 7;
  n_z_ = 3;
  n_sigma_ =  (2 * n_aug_) + 1;
  lambda_ = 3 - n_aug_;

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
  // clang-format off
  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;
  // clang-format on

  // initial sigma point predicted matrix
  Xsig_pred_ = MatrixXd(n_x_, n_sigma_);

  // Radar Measurement covariance matrix
  S_rad_ = MatrixXd(n_z_, n_z_);
  T_rad_ = MatrixXd(n_x_, n_z_);

  // Discard all but px and py when taking LIDAR measurements.
  H_lid_= MatrixXd(n_x_, n_x_);
  // clang-format off
  H_lid_ << 1, 0, 0, 0, 0,
            0, 1, 0, 0, 0,
            0, 0, 0, 0, 0,
            0, 0, 0, 0, 0,
            0, 0, 0, 0, 0;
  // clang-format on
  H_lid_t_ = H_lid_.transpose();
  I_ = MatrixXd::Identity(n_x_, n_x_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  // TODO 30ms2 is way way too high.
  std_a_ = 3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  // TODO yaw rate is way too high
  std_yawdd_ = M_PI / 16;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  // Predicted state in radar measurement space.
  z_rad_pred_ = VectorXd(n_z_);


  R_rad_ = MatrixXd(3, 3);
  R_rad_ << std_radr_ * std_radr_, 0, 0,
            0, std_radphi_ * std_radphi_ ,0,
            0, 0, std_radrd_ * std_radrd_;

  R_lid_ = MatrixXd(5, 5);
  R_lid_ << std_laspx_ * std_laspx_, 0, 0, 0, 0,
              0, std_laspy_* std_laspy_,0, 0, 0,
              0 ,0 ,0 ,0 ,0,
              0 ,0 ,0 ,0 ,0,
              0 ,0 ,0 ,0 ,0;

  weights_ = VectorXd(n_sigma_);
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<n_sigma_; i++) {
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }
//  std::cout << "weights=" << weights_ << std::endl;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage measurement_pack) {
  if (!is_initialized_) {

    // first measurement
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR &&
        use_radar_) {
      Eigen::VectorXd x_in(4);
//      cout << "UKF::ProcessMeasurement First Measurement RADAR" << endl;

      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      double ro = measurement_pack.raw_measurements_[0];
      double theta =
          Tools::NormalizeAngle(measurement_pack.raw_measurements_[1]);
      double px = cos(theta) * ro;
      double py = sin(theta) * ro;
      x_ << px, py, 0, 0 ,0;

      is_initialized_ = true;
    } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER &&
        use_laser_) {
//      cout << "UKF::ProcessMeasurement First Measurement LIDAR" << endl;

      /**
      Initialize state.dd
      */
      x_ << measurement_pack.raw_measurements_[0],
            measurement_pack.raw_measurements_[1], 0, 0, 0;
      is_initialized_ = true;
    }

    time_us_ = measurement_pack.timestamp_;
    // done initializing, no need to predict or update
    return;
  }


  double delta_t = (measurement_pack.timestamp_ - time_us_) / 1000000.0;
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR &&
      use_radar_) {
    Tools::PrintIn("UKF::ProcessMeasurement RADAR Measurement");
    PRINT(x_);
    PRINT(P_);

    // Prediction step.
    Prediction(delta_t);
    PredictRadarMeasurement(delta_t);
    UpdateRadar(delta_t, measurement_pack.raw_measurements_);
    time_us_ = measurement_pack.timestamp_;
    Tools::PrintOut("UKF::ProcessMeasurement RADAR Measurement");
  } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER &&
      use_laser_) {
    Tools::PrintIn("UKF::ProcessMeasurement LIDAR Measurement");
    Prediction(delta_t);
    VectorXd measurements(n_x_);
    measurements.fill(0.0);
    measurements(0) = measurement_pack.raw_measurements_(0);
    measurements(1) = measurement_pack.raw_measurements_(1);
    UpdateLidar(delta_t, measurements);
    time_us_ = measurement_pack.timestamp_;
    Tools::PrintOut("UKF::ProcessMeasurement LIDAR Measurement");
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  PredictSigmaPoints(delta_t);
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(double delta_t, const VectorXd& measurement) {
  Tools::PrintIn("UpdateLidar");

  PRINT(measurement);
  PRINT(x_);

  VectorXd y = measurement - H_lid_ * x_;
  PRINT(y);

  //calculate innovation covariance matrix S
  MatrixXd S(n_x_, n_x_);
  S.fill(0.0);
  for (int ii = 0; ii < n_sigma_; ++ii) {  //iterate over sigma points
    // state difference
    VectorXd x_diff = Xsig_pred_.col(ii) - x_;
    //angle normalization
    x_diff(3) = Tools::NormalizeAngle(x_diff(3));
    S = S + (weights_(ii) * x_diff * x_diff.transpose());
  }
  S = S + R_lid_;
  PRINT(S);

  MatrixXd K = P_ * H_lid_t_ * S.inverse();
  x_ = x_ + (K * y);
  x_(3) = Tools::NormalizeAngle(x_(3));
  P_ = (I_ - K * H_lid_) * P_;

  PRINT(x_);
  PRINT(P_);

  double nis = y.transpose() * S.inverse() *  y;
  lidar_nis_samples_.push_back(std::make_pair(delta_t, nis));

  Tools::PrintOut("UpdateLidar");
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(double delta_t, const VectorXd& measurement) {
  Tools::PrintIn("UpdateRadar");

  PRINT(z_rad_pred_);
  PRINT(measurement);

  VectorXd normalized_measurement(measurement);
  normalized_measurement(0) = measurement(0);
  normalized_measurement(1) = Tools::NormalizeAngle(measurement(1));
  normalized_measurement(2) = measurement(2);

  VectorXd z_diff = measurement - z_rad_pred_;
  z_diff(1) = Tools::NormalizeAngle(z_diff(1));

  PRINT(z_diff);
  MatrixXd S_rad_inverse = S_rad_.inverse();
  MatrixXd K = T_rad_ * S_rad_inverse;
  x_ = x_ + (K * z_diff);
  x_(3) = Tools::NormalizeAngle(x_(3));
  P_ = P_ - (K * S_rad_ * K.transpose());
  PRINT(measurement);

  PRINT(x_);
  PRINT(P_);

  double nis = z_diff.transpose() * S_rad_inverse *  z_diff;
  radar_nis_samples_.push_back(std::make_pair(delta_t, nis));
  Tools::PrintOut("UpdateRadar");
}

void UKF::PredictSigmaPoints(double delta_t) {
  Tools::PrintIn("PredictSigmaPoints");

  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sigma_);

  // Generate the Augmented sigma points.
  AugmentedSigmaPoints(&Xsig_aug);

  // Predict the sigma points using the process function.
  for (int col = 0 ; col < n_sigma_; ++col) {
    VectorXd c = Xsig_aug.col(col);
    VectorXd pred = Xsig_pred_.col(col);

    // Extract the components.
    double px = c(0);
    double py = c(1);
    double v = c(2);
    double psi = Tools::NormalizeAngle(c(3));
    double psi_dot = c(4);
    double v_a = c(5);
    double v_psi_dot_dot = c(6);

    double sin_psi = sin(psi);
    double cos_psi = cos(psi);
    double delta_t_2 = delta_t * delta_t;

    VectorXd pred_deterministic(5);
    VectorXd pred_stochastic(5);
    pred_stochastic << 0.5 * delta_t_2 * cos_psi * v_a,
                       0.5 * delta_t_2 * sin_psi * v_a,
                       delta_t * v_a,
                       Tools::NormalizeAngle(0.5 * delta_t_2 * v_psi_dot_dot),
                       delta_t * v_psi_dot_dot;

    // Generate the prediction for the column.
    if (fabs(psi_dot) > 0.001) {
      double v_over_psi_dot = v / psi_dot;
      pred_deterministic << v_over_psi_dot * (sin(psi + (psi_dot * delta_t)) - sin_psi),
                            v_over_psi_dot * (-cos(psi + (psi_dot * delta_t)) + cos_psi),
                            0,
                            Tools::NormalizeAngle(psi_dot * delta_t),
                            0;
    } else {
      pred_deterministic << v * cos_psi * delta_t,
                            v * sin_psi * delta_t,
                            0,
                            Tools::NormalizeAngle(psi_dot * delta_t),
                            0;
    }
    pred = c.head<5>();
//    PRINT(pred);
//    PRINT(pred_stochastic);
//    PRINT(pred_deterministic);

    pred += pred_deterministic + pred_stochastic;
    pred(3) = Tools::NormalizeAngle(pred(3));
//    PRINT(pred);
    Xsig_pred_.col(col) = pred;
  }
  PRINT(Xsig_pred_);
  PredictMeanAndCovarianceSigmaPoints();
  Tools::PrintOut("PredictSigmaPoints");
}

void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out) {
  Tools::PrintIn("AugmentedSigmaPoints");

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sigma_);
  Xsig_aug.setZero();

  MatrixXd Q = MatrixXd(2,2);
  Q << std_a_ * std_a_, 0,
       0, std_yawdd_ * std_yawdd_;

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.setZero();
  x_aug.head(n_x_) = x_;
  PRINT(x_aug);

  //create augmented covariance matrix
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.setZero();
  P_aug.topLeftCorner(5,5) = P_;
  P_aug.bottomRightCorner(2,2) = Q;
  PRINT(P_aug);

  //create square root matrix
  MatrixXd A = P_aug.llt().matrixL();
  A = A * sqrt(lambda_ + n_aug_);
  //create augmented sigma points
  for (int c = 0 ; c < n_sigma_ ; ++c) {
    if (c == 0) {
      Xsig_aug.col(c) = x_aug;
    } else if (c <= n_aug_) {
      Xsig_aug.col(c) = x_aug + A.col(c - 1);
    } else {
      Xsig_aug.col(c) = x_aug - A.col(c - 1 - n_aug_) ;
    }
  }
  *Xsig_out = Xsig_aug;
  PRINT(Xsig_aug);
  Tools::PrintOut("AugmentedSigmaPoints");
}

void UKF::PredictMeanAndCovarianceSigmaPoints() {
  Tools::PrintIn("PredictMeanAndCovarianceSigmaPoints");
  // Calculate predicted state mean
  VectorXd x_pred = VectorXd(n_x_);
  x_pred.fill(0.0);
  for (int i = 0; i < n_sigma_; i++) {
    x_pred = x_pred + (weights_(i) * Xsig_pred_.col(i));
  }
  x_pred(3) = Tools::NormalizeAngle(x_pred(3));
  PRINT(x_pred);

  // Calculate predicted state covariance matrix
  MatrixXd P_pred(n_x_, n_x_);
  P_pred.fill(0.0);
  for (int i = 0; i < n_sigma_; i++) {

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_pred;
    //angle normalization
    x_diff(3) = Tools::NormalizeAngle(x_diff(3));
    P_pred = P_pred + weights_(i) * x_diff * x_diff.transpose();
  }
  PRINT(P_pred);

  x_ = x_pred;
  P_ = P_pred;
  Tools::PrintOut("PredictMeanAndCovarianceSigmaPoints");
}

void UKF::PredictRadarMeasurement(double delta_t) {
  Tools::PrintIn("PredictRadarMeasurement");
  //create matrix for sigma points in the radar measurement space
  MatrixXd Zsig = MatrixXd(n_z_, n_sigma_);

  //transform sigma points into measurement space
  for (int ii = 0; ii < Xsig_pred_.cols(); ++ii) {
    Zsig.col(ii) = Tools::PositionSpaceToRadarMeasurementSpace(Xsig_pred_.col(ii));
  }
  PRINT(Zsig);
  //calculate mean predicted measurement
  z_rad_pred_.fill(0.0);
  for (int ii = 0; ii < n_sigma_; ++ii) {  //iterate over sigma points
    z_rad_pred_ = z_rad_pred_ + (weights_(ii) * Zsig.col(ii));
  }
  PRINT(z_rad_pred_);
  z_rad_pred_(1) = Tools::NormalizeAngle(z_rad_pred_(1));

  //calculate innovation covariance matrix S
  S_rad_.fill(0.0);
  for (int ii = 0; ii < n_sigma_; ++ii) {  //iterate over sigma points
    // state difference
    VectorXd z_diff = Zsig.col(ii) - z_rad_pred_;
    //angle normalization
    z_diff(1) = Tools::NormalizeAngle(z_diff(1));
    S_rad_ = S_rad_ + (weights_(ii) * z_diff * z_diff.transpose());
  }
  S_rad_ = S_rad_ + R_rad_;
  PRINT(S_rad_);
  // Calculate the Cross correlation matrix T between state space and measurement space.
  T_rad_.fill(0.0);
  for (int ii = 0; ii < n_sigma_; ++ii) {
    // state difference
    VectorXd x_diff = Xsig_pred_.col(ii) - x_;
    VectorXd z_diff = Zsig.col(ii) - z_rad_pred_;
    //angle normalization
    x_diff(3) = Tools::NormalizeAngle(x_diff(3));
    z_diff(1) = Tools::NormalizeAngle(z_diff(1));
    T_rad_ = T_rad_ + (weights_(ii) * x_diff * z_diff.transpose());
  }
  PRINT(T_rad_);
  Tools::PrintOut("PredictRadarMeasurement");
};

