#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"
#include <cmath>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  ekf_.P_ = MatrixXd(4, 4);

  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  ekf.F_ = MatrixXd(4, 4);
  ekf.F_ <<    1, 0, 0.5, 0,
          0, 1, 0, 0.5,
          0, 0, 1, 0,
          0, 0, 0, 1;


}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

  if (!is_initialized_) {


    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      float ro =  measurement_pack.raw_measurements_[0];
      float theta =  measurement_pack.raw_measurements_[1];
      float ro_dot =  measurement_pack.raw_measurements_[2];

      ekf_.P_ <<  1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1;

      float px = ro*cos(theta);
      float py = ro*sin(theta);
      float vx = ro_dot*cos(theta);
      float vy = ro_dot*sin(theta);
      if ( abs(px) < 0.0001 ) {
        px = 0.0001;
      }
      if ( abs(py) < 0.0001 ) {
        py = 0.0001;
      }
      ekf_.x_ << px, py, vx , vy;

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      ekf_.P_ <<  1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1000, 0,
            0, 0, 0, 1000;
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }

    previous_timestamp_ = measurement_pack.timestamp_;
    is_initialized_ = true;
    return;
  }

  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  double ax_noise = 9.0;
  double ay_noise = 9.0;

  double dt2 = dt * dt; //dt^2
  double dt3 = dt2 * dt; //dt^3
  double dt4 = dt3 * dt; //dt^4
  double dt4x4 = dt4 / 4; //dt^4/4
  double dt3x2 = dt3 / 2; //dt^3/2
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << dt4x4 * ax_noise, 0, dt3x2 * ax_noise, 0,
           0, dt4x4 * ay_noise, 0, dt3x2 * ay_noise,
           dt3x2 * ax_noise, 0, dt2 * ax_noise, 0,
           0, dt3x2 * ay_noise, 0, dt2 * ay_noise;

  ekf_.Predict();


  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
