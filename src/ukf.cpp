#include "ukf.h"
#include "Eigen/Dense"
#include "tools.h"

#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
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

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  // initial sigma point predicted matrix
  Xsig_pred_ = MatrixXd(n_x_, (2 * n_aug_) + 1);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  // TODO 30ms2 is way way to high.
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  // TODO yaw rate is way too
  std_yawdd_ = 30;
  
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
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */


  PredictSigmaPoints(delta_t);

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
}

void UKF::PredictSigmaPoints(double delta_t) {

  MatrixXd Xsig_aug = MatrixXd(n_aug_, (2 * n_aug_) + 1);

  // Generate the Augmented sigma points.
  AugmentedSigmaPoints(&Xsig_aug);

  // Predict the sigma points using the process function.
  for (int col = 0 ; col < (2 * n_aug_) + 1; ++col) {
    VectorXd c = Xsig_aug.col(col);
    VectorXd pred = Xsig_pred_.col(col);

    // Extract the components.
    double px = c[0];
    double py = c[1];
    double v = c[2];
    double psi = c[3];
    double psi_dot = c[4];
    double v_a = c[5];
    double v_psi_dot_dot = c[6];

    double sin_psi = sin(psi);
    double cos_psi = cos(psi);
    double v_over_psi_dot = v / psi_dot;
    if (psi_dot < 0.0001) continue;

    double delta_t_2 = delta_t * delta_t;

    VectorXd pred_deterministic(5);
    VectorXd pred_stochastic(5);
    pred_stochastic << 0.5 * delta_t_2 * cos_psi * v_a,
                       0.5 * delta_t_2 * sin_psi * v_a,
                       delta_t * v_a,
                       0.5 * delta_t_2 * v_psi_dot_dot,
                       delta_t * v_psi_dot_dot;

    // Generate the prediction for the column.
    if (psi_dot != 0) {
      pred_deterministic << v_over_psi_dot * (sin(psi + (psi_dot * delta_t)) - sin_psi),
                            v_over_psi_dot * (-cos(psi + (psi_dot * delta_t)) + cos_psi),
                            0,
                            psi_dot * delta_t,
                            0;
    } else {
      pred_deterministic << v * cos_psi * delta_t,
                            v * sin_psi * delta_t,
                            0,
                            psi_dot * delta_t,
                            0;
    }
    pred = c.head<5>() + pred_deterministic + pred_stochastic;
    Xsig_pred_.col(col) = pred;
  }
  PredicMeanAndCovarianceSigmaPoints();
}

void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out) {

    //define spreading parameter
    double lambda = 3 - n_aug_;

    //create example covariance matrix

    //create augmented mean vector
    VectorXd x_aug = VectorXd(n_aug_);
    //create augmented state covariance
    MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

    //create sigma point matrix
    MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_+ 1);

    MatrixXd Q = MatrixXd(2,2);
    Q << std_a_ * std_a_, 0,
        0, std_yawdd_ * std_yawdd_;

    //create augmented mean state
    x_aug.setZero();
    x_aug.head(n_aug_) = x_;
//   std::cout  << "x_aug" << std::endl << x_aug << std::endl;

//   std::cout  << "P_aug1" << std::endl << P_aug << std::endl;
    //create augmented covariance matrix
    P_aug.topLeftCorner(5,5) = P_;
//   std::cout  << "P_aug2" << std::endl << P_aug << std::endl;
    P_aug.bottomRightCorner(2,2) = Q;
//   std::cout  << "P_aug3" << std::endl << P_aug << std::endl;

    //create square root matrix
    MatrixXd A = P_aug.llt().matrixL();
//  std::cout  << "A" << std::endl << A << std::endl;
    A = A * sqrt(lambda + n_aug_);
    //create augmented sigma points
    for (int c = 0 ; c < ((n_aug_ * 2) +1) ; ++c) {
      if (c == 0) {
        Xsig_aug.col(c) = x_aug;
      } else if (c <= n_aug_) {
        Xsig_aug.col(c) = x_aug + A.col(c - 1);
      } else {
        Xsig_aug.col(c) = x_aug - A.col(c - 1 - n_aug_) ;
      }
    }
}

void UKF::PredicMeanAndCovarianceSigmaPoints() {
  int n_sigma = 2 * n_aug_ + 1;

  VectorXd w(Xsig_pred_.cols());
  for (int ii=0; ii < w.count(); ++ii) {
    if (ii == 0) {
      w[ii] = lambda_ / lambda_ + n_x_;
    } else {
      w[ii] = 1 / 2 * (lambda_ * n_x_);
    }
  }

  VectorXd x_pred = (w * Xsig_pred_);
  x_pred.fill(0.0);
  for (int i = 0; i < n_sigma; i++) {  //iterate over sigma points
    x_pred = x_pred + w(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  MatrixXd P_pred(n_x_, n_x_);
  P_pred.fill(0.0);
  for (int i = 0; i < n_sigma; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    x_diff(3) = Tools::NormalizeAngle(x_diff(3));
    P_pred = P_pred + w(i) * x_diff * x_diff.transpose() ;
  }

  std::cout << "x_pred = " << x_pred << std::endl;
  MatrixXd X = Xsig_pred_ - x_;
  std::cout << "p_pred = " << P_pred << std::endl;
  x_ = x_pred;
  P_ = P_pred;
}


