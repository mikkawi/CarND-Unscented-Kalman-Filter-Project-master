#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include "tools.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  
  n_x_ = 5;
  n_aug_ = 7;
  n_sigma_ = 2 * n_aug_ + 1;
  lambda_ = 3 - n_aug_;

  n_z_lidar_ = 2;
  n_z_radar_ = 3;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  //std_a_ = 30;
  std_a_ = 2.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  //std_yawdd_ = 30;
  std_yawdd_ = 0.7		;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  epsilon_ = 0.001;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  H_lidar_ = MatrixXd(n_z_lidar_, n_x_);
  R_lidar_ = MatrixXd(n_z_lidar_, n_z_lidar_);
  R_radar_ = MatrixXd(n_z_radar_, n_z_radar_);
  Q_ = MatrixXd::Zero(n_z_lidar_, n_z_lidar_);

  P_ = MatrixXd::Identity(n_x_, n_x_);

  H_lidar_ << 1, 0, 0, 0, 0,
              0, 1, 0, 0, 0;

  R_lidar_ << std_laspx_ * std_laspx_, 0,
              0, std_laspy_ * std_laspy_;

  R_radar_ << std_radr_ * std_radr_, 0, 0,
              0, std_radphi_ * std_radphi_, 0,
              0, 0, std_radrd_ * std_radrd_;

  Q_ << std_a_ * std_a_, 0,
        0, std_yawdd_ * std_yawdd_;

  NIS_lidar_ = 0;
  NIS_radar_ = 0;


  weights_ = VectorXd(n_sigma_);
  weights_(0)= lambda_/(lambda_+n_aug_);

  for (int i=1; i<n_sigma_; i++) {
    weights_(i) = 0.5/(n_aug_+lambda_);
  }
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  if (!is_initialized_){
      cout << "UKF: Not Initialized Yet" << endl;

      if(meas_package.sensor_type_ == MeasurementPackage::LASER){
       // cout << "UKF: Laser measurement " << endl;
        x_<< meas_package.raw_measurements_[0], meas_package.raw_measurements_[1],0,0,0;
        //cout << "UKF: x_ = " << x_ << endl;
      } else if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
        //cout << "UKF: Radar measurement" << endl;
        float px = meas_package.raw_measurements_[0] * cos(meas_package.raw_measurements_[1]);
        float py = meas_package.raw_measurements_[0] * sin(meas_package.raw_measurements_[1]);
            if (px == 0 && py == 0) {
            px = epsilon_;
            py = epsilon_;
            }
        x_ << px, py, 0,0,0;

      } 
      time_us_ = meas_package.timestamp_;
      is_initialized_ = true;
      cout << "UKF: " << endl;
      cout << "x_ initiatlized: " << is_initialized_ << endl;
      return;

  }

/* Prediction

*/
//cout << "Already initialized" << endl;
float dt = meas_package.timestamp_ - time_us_;
time_us_ = meas_package.timestamp_;

//cout << "dt = " << dt << endl;

Prediction (dt);

/* Update */

if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) {

	cout <<"UKF: Radar Measurement"<<endl;
	VectorXd z_pred;
    MatrixXd S, Zsig;

    PredictRadarMeasurement(Xsig_pred_, &z_pred, &S, &Zsig);
    //UpdateRadar(meas_package.raw_measurements_);
    //UpdateRadar(meas_package);
    UpdateRadar(Xsig_pred_, Zsig, z_pred, S, meas_package.raw_measurements_);
    NIS_lidar_ = 0;
    //return true;

} else if (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER) {
	cout<< "UKF: Lidar Measurement"<< endl;
	UpdateLidar(meas_package.raw_measurements_);
    NIS_radar_ = 0;
    //return true;

}
} 

/**
Predict Radar Measurements
**/

void UKF::PredictRadarMeasurement(MatrixXd Xsig_pred, VectorXd* z_pred_out, MatrixXd* S_out, MatrixXd* Zsig_out){
	MatrixXd Zsig 	= MatrixXd(n_z_radar_, n_sigma_);     	// Sigma points in measurement space
  	VectorXd z_pred = VectorXd(n_z_radar_);            		// Mean predicted measurement
  	MatrixXd S 		= MatrixXd(n_z_radar_, n_z_radar_);    	// Measurement covariance matrix S

  	for (int i=0;i<n_sigma_;i++) {
    	double px = Xsig_pred(0,i);
    	double py = Xsig_pred(1,i);
    	double v = Xsig_pred(2,i);
    	double yaw = Xsig_pred(3, i);

    	double vx = cos(yaw)*v;
    	double vy = sin(yaw)*v;
    	double rho = sqrt(px*px + py*py);

    	Zsig(0, i) = rho;
    	if (fabs(px) < APPROX_ZERO)
     		Zsig(1, i) = M_PI/2;       // Assume object is pointing straight up
    	else
      		Zsig(1, i) = atan2(py, px);

    	if (rho < APPROX_ZERO)
      		Zsig(2, i) = APPROX_ZERO;  // With zero rho assume rho dot is zero as well
    	else
      		Zsig(2, i) = (px * vx + py * vy) / rho;
  }

  // Calculate mean predicted measurement
  z_pred.fill(0);
  for(int i=0;i<n_sigma_;i++) {
  	z_pred+= weights_(i) * Zsig.col(i);
  }

  // Calculate measurement covariance matrix S
  S.fill(0);
  MatrixXd prod;
  for(int i=0;i<n_sigma_;i++) {
    prod = Zsig.col(i) - Zsig.col(0);

    // Normalise phi to be between -pi and pi
    prod(1) = Tools::constrainAngle(prod(1));

    S += weights_(i) * prod * prod.transpose();
  }
  S = S + R_radar_;


  // Write result
  *z_pred_out = z_pred;
  *S_out = S;
  *Zsig_out = Zsig;
}
/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  cout << "UKF: Precition Function" << endl;

  MatrixXd Xsig_aug;
  Xsig_aug = AugmentedSigmaPoints();
  SigmaPointPrediction(Xsig_aug, delta_t);
  PredictMeanAndCovariance();

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(VectorXd z) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  //cout << "UpdateLidar: H_lidar" << H_lidar_;
  VectorXd z_pred = H_lidar_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_lidar_.transpose();
  MatrixXd S = H_lidar_ * P_ * Ht + R_lidar_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  // New estimates for x and P
  x_ = x_ + (K * y);
  MatrixXd I = MatrixXd::Identity(n_x_, n_x_);
  P_ = (I - K * H_lidar_) * P_;

  // Calculate NIS
  NIS_lidar_ = y.transpose() * Si * y;



}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MatrixXd Xsig_pred, MatrixXd Zsig, VectorXd z_pred, MatrixXd S, VectorXd z) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  MatrixXd K, Tc;
  VectorXd y;

  // Create matrix for cross correlation Tc
  Tc = MatrixXd(n_x_, n_z_radar_);

  // Calculate cross correlation matrix
  Tc.fill(0);
  for (int i=0;i<n_sigma_;i++) {
    MatrixXd first = Xsig_pred.col(i) - Xsig_pred.col(0);
    first(3) = Tools::constrainAngle(first(3));

    MatrixXd second = Zsig.col(i) - Zsig.col(0);
    second(1) = Tools::constrainAngle(second(1));

    Tc += weights_(i) * first * second.transpose();
  }

  // Calculate Kalman gain K;
  K = Tc * S.inverse();

  // Update state mean and covariance matrix, normalize phi
  y = z - z_pred;
  y(1) = Tools::constrainAngle(y(1));

  // New estimates for x and P
  x_ = x_ + K * y;
  P_ = P_ - K * S * K.transpose();

  // Calculate NIS
  NIS_radar_ = y.transpose() * S.inverse() * y;
}

MatrixXd UKF::AugmentedSigmaPoints(){

  //cout << "UKF: AugmentedSigmaPoints function" << endl;
  // calculate labmda based on augmented dimention
  double lambda = 3 - n_aug_;

  //create augmented mean vector
  VectorXd x_aug = VectorXd(7);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;


  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;

  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda+n_aug_) * L.col(i);
  }
  return Xsig_aug;
}

void UKF::SigmaPointPrediction(MatrixXd Xsig_aug, double delta_t){
	//cout << "UKF: SigmaPointPrediction function" << endl;


	Xsig_pred_ = MatrixXd(n_x_,n_sigma_);
	Xsig_pred_.fill(0);

	for (int i=0;i < n_sigma_;i++) {
    double vk = Xsig_aug(2,i);
    double psi = Xsig_aug(3,i);
    double psidot = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_psidotdot = Xsig_aug(6,i);

    if (fabs(psidot) > 0.001) {
      Xsig_pred_(0, i) = Xsig_aug(0, i) + (vk / psidot) * (sin(psi + psidot * delta_t) - sin(psi)) +
                        0.5 * delta_t * delta_t * cos(psi) * nu_a;
      Xsig_pred_(1, i) = Xsig_aug(1, i) + (vk / psidot) * (-cos(psi + psidot * delta_t) + cos(psi)) +
                        0.5 * delta_t * delta_t * sin(psi) * nu_a;
    }
    else {
      Xsig_pred_(0, i) = Xsig_aug(0, i) + (vk * cos(psi)*delta_t) + 0.5 * delta_t * delta_t * cos(psi) * nu_a;
      Xsig_pred_(1, i) = Xsig_aug(1, i) + (vk * sin(psi)*delta_t) + 0.5 * delta_t * delta_t * sin(psi) * nu_a;
    }
    Xsig_pred_(2,i) = Xsig_aug(2, i) + 0 + delta_t*nu_a;
    Xsig_pred_(3,i) = Xsig_aug(3, i) + psidot*delta_t + 0.5*delta_t*delta_t*nu_psidotdot;
    Xsig_pred_(4,i) = Xsig_aug(4, i) + 0 + delta_t*nu_psidotdot;
  }



}

void UKF::PredictMeanAndCovariance(){
	//cout << "UKF: PredictMeanAndCovariance function" << endl;


	MatrixXd prod;
	x_.fill(0);
	//cout << "x_ before for loop" << x_ <<endl;
	//cout << "n_sigma_ = " <<n_sigma_ << endl;
	//cout << "weights_ = " << weights_ << endl;
	//cout << "Xsig_pred_.size(): "<<Xsig_pred_.rows()<<" X " <<Xsig_pred_.cols() <<endl;
	for(int i=0;i<n_sigma_;i++) {
    	x_ = x_ + weights_(i) * Xsig_pred_.col(i);
    }
    cout<< "x_ =" << x_ <<endl;
    P_.fill(0);
  	for(int i=0;i<n_sigma_;i++) {
    	prod = Xsig_pred_.col(i) - Xsig_pred_.col(0);

    // Normalise psi to be between -pi and pi
        prod(3) = Tools::constrainAngle(prod(3));

    P_ += weights_(i) * prod * prod.transpose(); 
  	}
  	cout << "P_ = "<< P_ << endl;

}
