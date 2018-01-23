#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  if (estimations.size() == 0 || estimations.size() != ground_truth.size()) {
    return rmse;
  }
  VectorXd sum(4);
  sum << 0, 0, 0, 0;
  // accumulate squared residuals
  for (int i = 0; i < estimations.size(); ++i) {
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array() * residual.array();
    sum += residual;
  }

  // calculate the mean
  VectorXd mean = sum.array() / estimations.size();

  // calculate the squared root
  rmse = mean.array().sqrt();

  // return the result
  return rmse;
}

double Tools::NormalizeAngle(double radians_in) {
  double normalized_angle = radians_in;
  double two_pi = 2 * M_PI;
  while (normalized_angle < -M_PI || normalized_angle > M_PI) {
    if (normalized_angle > M_PI) {
      normalized_angle -= two_pi;
    }
    if (normalized_angle < -M_PI) {
      normalized_angle += two_pi;
    }
  }
  return normalized_angle;
}