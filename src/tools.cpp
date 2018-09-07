#include <iostream>
#include "tools.h"
#include <math.h>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;
  
  // Calculate the sum of squared residuals
  for (size_t i = 0; i < estimations.size(); ++i){
//     VectorXd residual = estimations[i] - ground_truth[i];
//     residual = residual.array()*residual.array();
//     rmse += residual;
    rmse += VectorXd((estimations[i] - ground_truth[i]).array().pow(2));
  }
  // Calculate the mean sqrt
  rmse /= estimations.size();
  rmse = rmse.array().sqrt();
  
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  MatrixXd Hj(3, 4);
  
  double px = x_state(0);
  double py = x_state(1);
  double vx = x_state(2);
  double vy = x_state(3);
  
//   if (px == 0 && py ==0){
//     cout << "ERROR: Divided by zero!" << endl;
//     return Hj;
//   }
  
  double z = px * px + py * py;
  double z_sqrt = sqrt(z);
  double z_32 = pow(z, 1.5);
  double t = vx * py - vy * px;
  
  //check division by zero
  if(fabs(z) < 0.0001){
    cout << "CalculateJacobian () - Error - Division by Zero" << endl;
    return Hj;
  }
  
  Hj << px / z_sqrt, py / z_sqrt, 0, 0,
      -py / z, px / z, 0, 0,
      py * t / z_32, px * (-t) / z_32, px / z_sqrt, py / z_sqrt;
  
  return Hj;
}
