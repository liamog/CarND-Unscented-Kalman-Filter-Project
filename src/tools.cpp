#include "tools.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

static int current_verbosity = 1;

VectorXd Tools::CalculateRMSE(const vector /*unused*/<VectorXd> &estimations,
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
  constexpr double kTwoPI = 2.0 * M_PI;
  while (normalized_angle > M_PI) normalized_angle -=kTwoPI;
  while (normalized_angle <-M_PI) normalized_angle +=kTwoPI;
  return normalized_angle;
}

VectorXd Tools::PositionSpaceToRadarMeasurementSpace(const VectorXd &x) {
  VectorXd rm_space(3);
  rm_space.fill(0.0);
  // recover state parameters
  double px = x(0);
  double py = x(1);
  double v = x(2);
  double phi = x(3);

  double c1 = (px * px) + (py * py);
  if (c1 < 0.001) {
    return rm_space;
  }
  double c2 = sqrt(c1);

  double ro = c2;
  double theta = atan2(py, px);
  double ro_dot = ((px * cos(phi) * v) + (py * sin(phi) * v)) / c2;
  rm_space << ro, theta, ro_dot;
  return rm_space;
}

void Tools::SetVerbosity(int level) {
  current_verbosity = level;
}

void Tools::Print(int verbosity, int value, const string &name) {
  if (verbosity > current_verbosity) return;
  std::cout << name << "=" << value << std::endl;
}

void Tools::Print(int verbosity, double value, const string &name) {
  if (verbosity > current_verbosity) return;
  std::cout << name << "=" << value << std::endl;
}

void Tools::Print(int verbosity, const MatrixXd &matrix, const string &name) {
  if (verbosity > current_verbosity) return;
  Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
  std::cout << name << "=" << std::endl << matrix.format(CleanFmt) << std::endl;
}

void Tools::Print(int verbosity, const VectorXd &vector, const string &name) {
  if (verbosity > current_verbosity) return;
  Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
  std::cout << name << "=" << vector.transpose().format(CleanFmt) << std::endl;
}

void Tools::PrintIn(int verbosity, const string &name) {
  if (verbosity > current_verbosity) return;
  std::cout << ">>>>>>>>>>>> " << name << std::endl;
}

void Tools::PrintOut(int verbosity, const string &name) {
  if (verbosity > current_verbosity) return;
  std::cout << "<<<<<<<<<<<< " << name << std::endl;
}