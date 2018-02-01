#ifndef _USERS_LIAM_DEV_UDACITY_CARND_UNSCENTED_KALMAN_FILTER_PROJECT_SRC_TOOLS_H
#define _USERS_LIAM_DEV_UDACITY_CARND_UNSCENTED_KALMAN_FILTER_PROJECT_SRC_TOOLS_H
#include <vector>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

class Tools {
 public:
  /**
   * A helper method to calculate RMSE.
   */
  static VectorXd CalculateRMSE(const vector<VectorXd> &estimations,
                                const vector<VectorXd> &ground_truth);

  /**
   * A helper method to normalize an angle between -2pi and +2pi
   */
  static double NormalizeAngle(double radians_in);

  /**
   * Converts vector from state space to Radar measurement space.
   * @param x state space
   * @return Vector with 3 elements per Radar measurement.
   */
  static VectorXd PositionSpaceToRadarMeasurementSpace(const VectorXd &x);

  static void SetVerbosity(int level);
  static void Print(int verbosity, int value, const string &name);
  static void Print(int verbosity, double value, const string &name);
  static void Print(int verbosity, const MatrixXd &matrix, const string &name);
  static void Print(int verbosity, const VectorXd &vector, const string &name);
  static void PrintIn(int verbosity, const string &name);
  static void PrintOut(int verbosity, const string &name);
};

#endif /* TOOLS_H_ */