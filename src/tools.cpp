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

VectorXd Tools::PositionSpaceToRadarMeasurementSpace(const VectorXd &x) {
  VectorXd rm_space(3);
  rm_space.fill(0.0);
  // recover state parameters
  float px = x(0);
  float py = x(1);
  float v = x(2);
  float phi = x(3);

  float c1 = (px * px) + (py * py);
  if (c1 < 0.00001) return rm_space;
  float c2 = sqrt(c1);

  double ro = c2;
  double theta = atan2(py, px);
  double ro_dot = ((px * cos(phi) * v) + (py * sin(phi) * v)) / c2;
  rm_space << ro, theta, ro_dot;
  return rm_space;
}

void Tools::Print(const MatrixXd &matrix, const string &name) {
  Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
  std::cout << name << "=" << std::endl << matrix.format(CleanFmt) << std::endl;
}

void Tools::Print(const VectorXd &vector, const string &name) {
  Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
  std::cout << name << "=" << std::endl << vector.transpose().format(CleanFmt) << std::endl;
}

void Tools::PrintIn(const string &name) {
  std::cout << ">>>>>>>>>>>> " << name << std::endl;
}

void Tools::PrintOut(const string &name) {
  std::cout << "<<<<<<<<<<<< " << name << std::endl;
}