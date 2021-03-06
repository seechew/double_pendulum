/* g++ velocity_verlet.cpp -o velocity_verlet.out -larmadillo */
/*
   Reference: https://scienceworld.wolfram.com/physics/DoublePendulum.html   
*/

#include <iostream>
#include <cmath>
#include <armadillo>
#include <chrono>
#include <iomanip>  //std::setprecision
#include <fstream>

const double m1 = 10;
const double m2 = 10;
const double l1 = 10;
const double l2 = 10;
const double g = 9.8067;
const double pi = 3.14159265358979323846;

void equation_of_motion(const arma::vec& theta, const arma::vec& theta_dot, arma::vec& theta_dot_dot) {

  double s = sin(theta[0]-theta[1]);
  double c = cos(theta[0]-theta[1]);

  double z1 = theta_dot[0];
  double z2 = theta_dot[1];

  double theta1 = theta[0];
  double theta2 = theta[1];

  double z1dot = (m2*g*sin(theta2)*c - m2*s*(l1*pow(z1,2)*c + l2*pow(z2,2)) -
		 (m1+m2)*g*sin(theta1)) / l1 / (m1 + m2*pow(s,2));
 
  double z2dot = ((m1+m2)*(l1*pow(z1,2)*s - g*sin(theta2) + g*sin(theta1)*c) + 
		 m2*l2*pow(z2,2)*s*c) / l2 / (m1 + m2*pow(s,2));

  theta_dot_dot[0] = z1dot;
  theta_dot_dot[1] = z2dot;
  
}

// -------------------------------------------------------------------------------
int main() {

  std::ofstream csv("data_velocity_verlet.csv");

  auto start = std::chrono::steady_clock::now();
 
  arma::vec theta(2);
  arma::vec theta_dot(2);
  arma::vec theta_dot_dot(2);

  /* INITIAL CONDITIONS */
  theta[0] = pi/2;
  theta[1] = pi;
  theta_dot[0] = 0.0;
  theta_dot[1] = 0.0;

  /* INITIAL ACCELERATION */
  theta_dot_dot[0] = 0.0;
  theta_dot_dot[1] = 0.0;
  
  double t0 = 0.0;
  double dt = 1e-2;
  double tf = 7200;
  double current_t = t0 + dt;  // INCREMENT 1
  bool firstStep = true;
  arma::vec theta_dot_mid(2,arma::fill::zeros);

  double V, T, H;
  double x1;
  double x2;
  double y1;
  double y2;

  /* WRITE HEADER TO FILE */
  csv << "time,x2,y2,H" << std::endl;
  
  while (current_t < tf+dt) {

    /* CALLING EQUATION OF MOTION */
    equation_of_motion(theta, theta_dot_mid, theta_dot_dot);
    
    theta_dot_mid = 0.5*theta_dot;
    theta_dot += theta_dot_dot * dt * (firstStep ? 0.5 : 1);
    theta_dot_mid += 0.5 * theta_dot; // get theta-dot at integer time
    theta += theta_dot *dt;
    firstStep = false;

    /* PRINT HAMILTONIAN */
    double th1  = theta[0];
    double th2  = theta[1];
    double th1d = theta_dot[0];
    double th2d = theta_dot[1];
    
    V = -(m1+m2)*l1*g*cos(th1) - m2*l2*g*cos(th2);
    
    T = 0.5*m1*pow(l1*th1d,2) + 0.5*m2*(pow(l1*th1d,2) + pow(l2*th2d,2) +
	2*l1*l2*th1d*th2d*cos(th1-th2));
    
    H = T + V;
    
    x1 =  l1*sin(theta[0]);
    y1 = -l1*cos(theta[0]);
    x2 =  x1 + l2*sin(theta[1]);
    y2 =  y1 - l2*cos(theta[1]);
    
    csv << current_t << "," << x2 << "," << y2 << "," << H << std::endl;
      
    current_t += dt;
  }
  
  auto end = std::chrono::steady_clock::now();
  
  std::cout << "DONE, ELAPSED TIME(s):  " << std::chrono::duration_cast<std::chrono::seconds>(end-start).count() << std::endl;

  return 0;
}
