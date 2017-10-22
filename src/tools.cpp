#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

//Tools::Tools() {}

//Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
    VectorXd rmse(4);
    rmse << 0, 0, 0, 0;

    // Check inputs
    if (estimations.size() != ground_truth.size()
        || estimations.size() == 0) {
      cout << "Invalid estimation or ground_truth data" << endl;
      return rmse;
    }

    // Calculate RMSE
    for (int i = 0; i < estimations.size(); ++i) {
      VectorXd res = estimations[i] - ground_truth[i];
      res = res.array() * res.array();
      rmse += res;
    }

    rmse /= estimations.size();
    rmse = rmse.array().sqrt();

    return rmse;
}

double Tools::constrainAngle(double angle){
    angle = fmod(angle + M_PI,2*M_PI);
    if (angle < 0)
      angle += 2*M_PI;
    return angle - M_PI;
  }