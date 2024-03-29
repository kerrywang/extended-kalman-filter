#ifndef TOOLS_H_
#define TOOLS_H_

#include <vector>
#include "Eigen/Dense"

namespace Tools {
  /**
   * A helper method to calculate RMSE.
   */
  Eigen::VectorXd CalculateRMSE(const std::vector<Eigen::VectorXd> &estimations, 
                                const std::vector<Eigen::VectorXd> &ground_truth);

  /**
   * A helper method to calculate Jacobians.
   */
  Eigen::MatrixXd CalculateJacobian(const Eigen::VectorXd& x_state);

  /**
   * A helper method to convert cartesian coordinate to polar coordinate
   */
  Eigen::VectorXd CartesianToPolar(const Eigen::VectorXd& cart_coor);
};

#endif  // TOOLS_H_
