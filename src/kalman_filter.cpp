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
  /**
  TODO:
    * predict the state
  */
	x_ = F_ * x_;
	P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
	
	VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
	double sum_square = x_[0]*x_[0] + x_[1]*x_[1];
	double rho = sqrtf(sum_square);
	double phi;
	// check if px is too small
	if (fabs(x_[0]) > 0.001) {
		phi = atan2f(x_[1], x_[0]);
	}
	else {
		phi = 0.001;
	}
	
	
	double rho_dot;
	// check if rho is too small
	if (fabs(rho) > 0.001) {
		rho_dot = double(x_[0] * x_[2] + x_[1] * x_[3]) / rho;
	}
	else {
		rho = 0.0;
		rho_dot = 0.0;
	}
	

	VectorXd h = VectorXd(3);
	h << rho, phi, rho_dot;
	VectorXd y = z - h;
	
	//normalize atan
	while (y[1] > M_PI) {
		y[1] -= 2.* M_PI;
	}
	while (y[1] < -M_PI) {
		y[1] += 2.* M_PI;
	}


	//H_ is updated with Jakobian
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;

}
