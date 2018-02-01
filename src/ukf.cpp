#include "ukf.h"
#include "Eigen/Dense"
#include "tools.h"

#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

#define PRINT(l, x) Tools::Print(l, x, #x);

// TUNABLE
// Process noise standard deviation longitudinal acceleration in m/s^2
constexpr double kStdA = 2.0;

// TUNABLE
// Process noise standard deviation yaw acceleration in rad/s^2
constexpr double kStdYawdd = M_PI / 3;

/////////////////////////////////////////////////////////////////////////
// DO NOT MODIFY measurement noise values above these are provided by the
// sensor manufacturer.

// Laser measurement noise standard deviation position1 in m
constexpr double kStdLasPx = 0.15;

// Laser measurement noise standard deviation position2 in m
constexpr double kStdLasPy = 0.15;

// Radar measurement noise standard deviation radius in m
constexpr double kStdRadR = 0.3;

// Radar measurement noise standard deviation angle in rad
constexpr double kStdRadPhi = 0.03;

// Radar measurement noise standard deviation radius change in m/s
constexpr double kStdRadRd = 0.3;
/////////////////////////////////////////////////////////////////////////

// Number of state dimensions.
constexpr int kNx = 5;

// Number of radar measurement dimensions.
constexpr int kNz = 3;

// Number of augmented dimensions.
constexpr int kNaug = 7;

// Number of sigma points.
constexpr int kNSigma = ((2 * kNaug) + 1);

// Lambda design parameter.
constexpr double kLambda =  3 - kNaug;

constexpr bool kUseRadar = true;
constexpr bool kUseLidar = true;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF()
    : x_(kNx),
      Q_(2,2),
      P_(kNx, kNx),
      Xsig_pred_(kNx, kNSigma),
      S_rad_(kNz, kNz),
      T_rad_(kNx, kNz),
      H_lid_(2, kNx),
      I_(MatrixXd::Identity(kNx, kNx)),
      z_rad_pred_(kNz),
      R_rad_(3, 3),
      R_lid_(2,2),
      weights_(kNSigma),
      radar_nis_samples_("radar_nis.csv", std::fstream::out),
      lidar_nis_samples_("lidar_nis.csv", std::fstream::out) {

  // initial state vector
  x_.fill(0.0);

  // initial covariance matrix
  // clang-format off
  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;
  // clang-format on

  // Discard all but px and py when taking LIDAR measurements.
  // clang-format off
  H_lid_ << 1, 0, 0, 0, 0,
            0, 1, 0, 0, 0;
  // clang-format on
  H_lid_t_ = H_lid_.transpose();

  // clang-format off
  R_rad_ << kStdRadR * kStdRadR, 0, 0,
            0, kStdRadPhi * kStdRadPhi ,0,
            0, 0, kStdRadRd * kStdRadRd;
  // clang-format on

  // clang-format off
  R_lid_ << kStdLasPx * kStdLasPx, 0,
            0, kStdLasPy* kStdLasPy;
  // clang-format on

  // clang-format off
  Q_ << kStdA * kStdA, 0,
      0, kStdYawdd * kStdYawdd;
  // clang-format on

  weights_(0) = kLambda / (kLambda + kNaug);
  for (int i = 1; i < kNSigma; i++) {
    double weight = 0.5 / (kNaug + kLambda);
    weights_(i) = weight;
  }
}

UKF::~UKF() = default;

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage measurement_pack) {
  counter_++;
  Tools::SetVerbosity(1);
  PRINT(3, counter_);

  if (!is_initialized_) {
    PRINT(1, kStdA);
    PRINT(1, kStdYawdd);
    PRINT(1, P_);

    // first measurement
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR &&
        kUseRadar) {
      Eigen::VectorXd x_in(4);
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      double ro = measurement_pack.raw_measurements_[0];
      double theta =
          Tools::NormalizeAngle(measurement_pack.raw_measurements_[1]);
      double px = cos(theta) * ro;
      double py = sin(theta) * ro;
      x_ << px, py, 0, 0, 0;

      is_initialized_ = true;
    } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER &&
               kUseLidar) {
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

  current_timestamp_ = measurement_pack.timestamp_ / 1000000.0;
  double delta_t = (measurement_pack.timestamp_ - time_us_) / 1000000.0;
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR &&
      kUseRadar) {
    Tools::PrintIn(3, "UKF::ProcessMeasurement RADAR Measurement");
    PRINT(2, x_);
    PRINT(2, P_);

    Prediction(delta_t);
    PredictRadarMeasurement();
    UpdateRadar(measurement_pack.raw_measurements_);
    time_us_ = measurement_pack.timestamp_;
    Tools::PrintOut(3, "UKF::ProcessMeasurement RADAR Measurement");
  } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER &&
             kUseLidar) {
    Tools::PrintIn(3, "UKF::ProcessMeasurement LIDAR Measurement");
    PRINT(2, x_);
    PRINT(2, P_);

    Prediction(delta_t);
    VectorXd measurements(2);
    measurements.fill(0.0);
    measurements(0) = measurement_pack.raw_measurements_(0);
    measurements(1) = measurement_pack.raw_measurements_(1);
    UpdateLidar(measurements);
    time_us_ = measurement_pack.timestamp_;
    Tools::PrintOut(3, "UKF::ProcessMeasurement LIDAR Measurement");
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) { PredictSigmaPoints(delta_t); }

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(const VectorXd &measurement) {
  Tools::PrintIn(3, "UpdateLidar");

  PRINT(2, measurement);
  PRINT(2, x_);

  VectorXd y = measurement - H_lid_ * x_;
  PRINT(2, y);
  PRINT(2, H_lid_);
  PRINT(2, H_lid_t_);
  PRINT(2, R_lid_);

  MatrixXd S = H_lid_ * P_ * H_lid_t_ + R_lid_;
  PRINT(2, S);

  MatrixXd K = P_ * H_lid_t_ * S.inverse();
  x_ = x_ + (K * y);
  x_(3) = Tools::NormalizeAngle(x_(3));
  MatrixXd P_prime = (I_ - K * H_lid_) * P_;
  MatrixXd P_diff = P_ - P_prime;
  P_ = P_prime;

  PRINT(2, x_);
  PRINT(2, P_);

  double nis = y.transpose() * S.inverse() * y;
  lidar_nis_samples_ << current_timestamp_ << "," << nis << std::endl;
  Tools::PrintOut(3, "UpdateLidar");
}


void UKF::PredictSigmaPoints(double delta_t) {
  Tools::PrintIn(3, "PredictSigmaPoints");

  MatrixXd Xsig_aug = MatrixXd(kNaug, kNSigma);

  // Generate the Augmented sigma points.
  GenerateAugmentedSigmaPoints(&Xsig_aug);

  // Predict the sigma points using the process function.
  for (int col = 0; col < kNSigma; ++col) {
    VectorXd sigma_point = Xsig_aug.col(col);
    VectorXd pred_sigma = Xsig_pred_.col(col);

    // Extract the components.
    double px = sigma_point(0);
    double py = sigma_point(1);
    double v = sigma_point(2);
    double psi = sigma_point(3);
    double psi_dot = sigma_point(4);
    double v_a = sigma_point(5);
    double v_psi_dot_dot = sigma_point(6);

    double sin_psi = sin(psi);
    double cos_psi = cos(psi);
    double delta_t_2 = delta_t * delta_t;

    VectorXd pred_deterministic(5);
    VectorXd pred_stochastic(5);
    pred_stochastic << 0.5 * delta_t_2 * cos_psi * v_a,
        0.5 * delta_t_2 * sin_psi * v_a, delta_t * v_a,
        0.5 * delta_t_2 * v_psi_dot_dot,
        delta_t * v_psi_dot_dot;

    // Generate the prediction for the column.
    if (fabs(psi_dot) > 0.0001) {
      double v_over_psi_dot = v / psi_dot;
      pred_deterministic << v_over_psi_dot *
                                (sin(psi + (psi_dot * delta_t)) - sin_psi),
          v_over_psi_dot * (-cos(psi + (psi_dot * delta_t)) + cos_psi), 0,
          psi_dot * delta_t, 0;
    } else {
      pred_deterministic << v * cos_psi * delta_t, v * sin_psi * delta_t, 0,
          psi_dot * delta_t, 0;
    }
    pred_sigma = sigma_point.head<5>();
    pred_sigma += (pred_deterministic + pred_stochastic);
    Xsig_pred_.col(col) = pred_sigma;
  }
  PRINT(2, Xsig_pred_);
  PredictMeanAndCovarianceSigmaPoints();
  Tools::PrintOut(3, "PredictSigmaPoints");
}

void UKF::GenerateAugmentedSigmaPoints(MatrixXd* Xsig_out) {
  Tools::PrintIn(3, "AugmentedSigmaPoints");

  // create sigma point matrix
  MatrixXd Xsig_aug(kNaug, kNSigma);
  Xsig_aug.setZero();

  // create augmented mean vector
  VectorXd x_aug(kNaug);
  x_aug.setZero();
  x_aug.head(kNx) = x_;
  PRINT(2, x_aug);

  // create augmented covariance matrix
  MatrixXd P_aug(kNaug, kNaug);
  P_aug.setZero();
  P_aug.topLeftCorner(5, 5) = P_;
  P_aug.bottomRightCorner(2, 2) = Q_;
  PRINT(2, P_aug);

  // create square root matrix
  MatrixXd A = P_aug.llt().matrixL();
  A = A * sqrt(kLambda + kNaug);
  // create augmented sigma points
  for (int c = 0; c < kNSigma; ++c) {
    if (c == 0) {
      Xsig_aug.col(c) = x_aug;
    } else if (c <= kNaug) {
      Xsig_aug.col(c) = x_aug + A.col(c - 1);
    } else {
      Xsig_aug.col(c) = x_aug - A.col(c - 1 - kNaug);
    }
  }
  *Xsig_out = Xsig_aug;
  PRINT(2, Xsig_aug);
  Tools::PrintOut(3, "AugmentedSigmaPoints");
}

void UKF::PredictMeanAndCovarianceSigmaPoints() {
  Tools::PrintIn(3, "PredictMeanAndCovarianceSigmaPoints");
  // Calculate predicted state mean
  VectorXd x_pred(kNx);
  x_pred.fill(0.0);
  for (int i = 0; i < kNSigma; i++) {
    x_pred = x_pred + weights_(i) * Xsig_pred_.col(i);
  }
  PRINT(2, x_pred);

  // Calculate predicted state covariance matrix
  MatrixXd P_pred(kNx, kNx);
  P_pred.fill(0.0);
  for (int i = 0; i < kNSigma; i++) {
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_pred;
    // angle normalization
    x_diff(3) = Tools::NormalizeAngle(x_diff(3));
    P_pred = P_pred + weights_(i) * x_diff * x_diff.transpose();
    //P_pred(3,3) = Tools::NormalizeAngle(P_pred(3,3));
  }
  PRINT(2, P_pred);

  x_pred(3) = Tools::NormalizeAngle(x_pred(3));
  x_ = x_pred;
  P_ = P_pred;
  Tools::PrintOut(3, "PredictMeanAndCovarianceSigmaPoints");
}

void UKF::PredictRadarMeasurement() {
  Tools::PrintIn(3, "PredictRadarMeasurement");
  // create matrix for sigma points in the radar measurement space
  MatrixXd Zsig = MatrixXd(kNz, kNSigma);
  Zsig.fill(0.0);
  // transform sigma points into measurement space
  for (int ii = 0; ii < Xsig_pred_.cols(); ++ii) {
    Zsig.col(ii) =
        Tools::PositionSpaceToRadarMeasurementSpace(Xsig_pred_.col(ii));
  }
  PRINT(2, Zsig);
  // calculate mean predicted measurement
  z_rad_pred_.fill(0.0);
  for (int ii = 0; ii < kNSigma; ++ii) {  // iterate over sigma points
    z_rad_pred_ = z_rad_pred_ + (weights_(ii) * Zsig.col(ii));
//    z_rad_pred_(1) = Tools::NormalizeAngle(z_rad_pred_(1));
  }
  PRINT(2, z_rad_pred_);

  // calculate innovation covariance matrix S
  S_rad_.fill(0.0);
  for (int ii = 0; ii < kNSigma; ++ii) {  // iterate over sigma points
    // state difference
    VectorXd z_diff = Zsig.col(ii) - z_rad_pred_;
    // angle normalization
    z_diff(1) = Tools::NormalizeAngle(z_diff(1));
    S_rad_ = S_rad_ + (weights_(ii) * z_diff * z_diff.transpose());
  }
  S_rad_ = S_rad_ + R_rad_;
  PRINT(2, S_rad_);

  // Calculate the Cross correlation matrix T between state space and
  // measurement space.
  T_rad_.fill(0.0);
  for (int ii = 0; ii < kNSigma; ++ii) {
    // state difference
    VectorXd x_diff = Xsig_pred_.col(ii) - x_;
    VectorXd z_diff = Zsig.col(ii) - z_rad_pred_;
    // angle normalization
    x_diff(3) = Tools::NormalizeAngle(x_diff(3));
    z_diff(1) = Tools::NormalizeAngle(z_diff(1));
    T_rad_ = T_rad_ + weights_(ii) * x_diff * z_diff.transpose();
  }
  PRINT(2, T_rad_);
  Tools::PrintOut(3, "PredictRadarMeasurement");
};

void UKF::UpdateRadar(const VectorXd &measurement) {
  Tools::PrintIn(3, "UpdateRadar");

  PRINT(2, z_rad_pred_);
  PRINT(2, measurement);

  VectorXd normalized_measurement(measurement);
  normalized_measurement(0) = measurement(0);
  normalized_measurement(1) = Tools::NormalizeAngle(measurement(1));
  normalized_measurement(2) = measurement(2);

  VectorXd z_diff = measurement - z_rad_pred_;
  z_diff(1) = Tools::NormalizeAngle(z_diff(1));

  PRINT(2, z_diff);
  MatrixXd S_rad_inverse = S_rad_.inverse();
  MatrixXd K = T_rad_ * S_rad_inverse;
  x_ = x_ + K * z_diff;
  x_(3) = Tools::NormalizeAngle(x_(3));
  P_ = P_ - K * S_rad_ * K.transpose();
  P_(3,3) = Tools::NormalizeAngle(P_(3,3));

  PRINT(2, measurement);

  PRINT(2, x_);
  PRINT(2, P_);

  double nis = z_diff.transpose() * S_rad_inverse * z_diff;
  radar_nis_samples_ << current_timestamp_ << "," << nis << std::endl;
  Tools::PrintOut(3, "UpdateRadar");
}