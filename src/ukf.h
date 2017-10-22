#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  const double APPROX_ZERO = 0.0001;
  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;

  MatrixXd Q_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  ///* time when the state is true, in us
  long long time_us_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  /// Number of sigma points
  int n_sigma_;

  int n_z_radar_;
  int n_z_lidar_;


  ///* Sigma point spreading parameter
  double lambda_;

  float epsilon_;


  double NIS_lidar_;
  double NIS_radar_;

  // Measurement covariance matrix - radar measurement noise
  MatrixXd R_radar_; 
  // Measurement covariance matrix - laser measurement noise     
  MatrixXd R_lidar_;
  // Measurement function H matrix      
  MatrixXd H_lidar_;    


  //VectorXd z_pred_;
  //MatrixXd S_;
  //MatrixXd Zsig_;

  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  //void UpdateLidar(MeasurementPackage meas_package);
   void UpdateLidar(VectorXd z);
  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MatrixXd Xsig_pred, MatrixXd Zsig, VectorXd z_pred, MatrixXd S, VectorXd z);

  // Generates augmented sigma points
  MatrixXd AugmentedSigmaPoints();

  // Predicts sigma points using the process model
  void SigmaPointPrediction(MatrixXd Xsig_aug, double delta_t);

  // Predicts state mean and covariance
  void PredictMeanAndCovariance();

  // Predicts radar measurement mean and covariance using the measurement model
  void PredictRadarMeasurement(MatrixXd Xsig_pred, VectorXd* z_pred_out, MatrixXd* S_out, MatrixXd* Zsig_out);


  
};

#endif /* UKF_H */
