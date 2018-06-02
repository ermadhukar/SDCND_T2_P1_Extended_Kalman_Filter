#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
    x_ = x_in;
    P_ = P_in;
    F_ = F_in;
    H_ = H_in;
    R_ = R_in;
    Q_ = Q_in;
}

void KalmanFilter::Predict() {
    /* predict the state */
    
    x_ = F_ * x_; // No '+u' as there is no external motion
    P_ = F_ * P_ * F_.transpose() + Q_;

}

void KalmanFilter::Update(const VectorXd &z) {
    /* update the state by using Kalman Filter equations */

    VectorXd y = z - H_ * x_;
    MatrixXd PHt = P_ * H_.transpose();
    MatrixXd S = H_ * PHt + R_;
    MatrixXd K = PHt * S.inverse();

    x_ = x_ + (K * y);
    int x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}


void KalmanFilter::UpdateEKF(const VectorXd &z) {
    
/* update the state by using Extended Kalman Filter equations */

   
    double px = x_(0);
    double py = x_(1);
    double vx = x_(2);
    double vy = x_(3);

    /* Radar measurement space */

    double rho = sqrt(px * px + py * py);
    double theta = atan2(py, px);
    double rho_dot = (rho != 0 ? (px * vx + py * vy) / rho : 0);

    
    /* position-speed prediction */
    
    VectorXd z_pred(3);
    z_pred << rho, theta, rho_dot;

    /* update measurement */

    VectorXd y = z - z_pred;

    /* Angle normalization */

    double width = 2 * M_PI; 
    double offsetValue = y(1) + M_PI; 
    y(1) = (offsetValue - (floor(offsetValue / width) * width)) - M_PI;

    MatrixXd PHt = P_ * H_.transpose();
    MatrixXd S = H_ * PHt + R_;
    MatrixXd K = PHt * S.inverse();

    /* New State */

    x_ = x_ + (K * y);
    int x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}
