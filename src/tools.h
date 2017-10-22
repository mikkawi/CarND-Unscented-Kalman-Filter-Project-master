#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

namespace Tools {
//public:
  /**
  * Constructor.
  */
  //Tools();

  /**
  * Destructor.
  */
  //virtual ~Tools();

  /**
  * A helper method to calculate RMSE.
  */
  Eigen::VectorXd CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth);

    // Constrain angle to -PI to PI
  double constrainAngle(double angle);

};

#endif /* TOOLS_H_ */