#include "../include/tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;


VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */
   VectorXd rmse(4);
   rmse << 0, 0, 0, 0;
   size_t size = estimations.size();
   for (size_t i = 0; i < size; i++) {
       rmse += static_cast<VectorXd>((ground_truth[i] - estimations[i]).array().pow(2));
   }
   return (rmse / size).array().sqrt();
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
    MatrixXd Hj(3,4);
    // recover state parameters
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);

    // pre-compute a set of terms to avoid repeated calculation
    float c1 = px*px+py*py;
    float c2 = sqrt(c1);
    float c3 = (c1*c2);

    // check division by zero
    if (fabs(c1) < 0.0001) {
        std::cout << "CalculateJacobian () - Error - Division by Zero" << std::endl;
        return Hj;
    }

    // compute the Jacobian matrix
    Hj << (px/c2), (py/c2), 0, 0,
            -(py/c1), (px/c1), 0, 0,
            py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

    return Hj;
}

Eigen::VectorXd Tools::CartesianToPolar(const Eigen::VectorXd& cart_coor) {
    VectorXd polar(3);
    double px = cart_coor[0], py = cart_coor[1], vx = cart_coor[2], vy = cart_coor[3];
    double rho = 0, phi = 0, rho_vel = 0;
    rho = sqrt(px * px + py * py);
    if (fabs(rho) > 0.0001) {
        phi = std::atan2(py, px);
        rho_vel = (px * vx + py * vy) / rho;
    }
    polar << rho, phi, rho_vel;
    return polar;
}

