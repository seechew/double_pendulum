/* g++ velocity_verlet.cpp -o velocity_verlet -larmadillo */
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

  /* DEFINE THE MASS MATRIX
  mass_matrix = [[(m1+m2)*l1, m2*l2*cos(theta_1 - theta_2)];[m2*l1*cos(theta_1 - theta_2), m2*l2]]
  
  // DEFINE THE RHS
  rhs = [[-2*g*(m1+m2)*sin(theta_1) - m2*l2*theta_dot_2^2*sin(theta_1 - theta_2];
         [  -g*m2*sin(theta_2) + m2*l1*theta_dot_1^2*sin(theta_1 - theta_2]]

  */

  /* [a b;
     c d]   */
  double a = (m1+m2)*l1;
  double b = m2*l2*cos(theta[0] - theta[1]);
  double c = m2*l1*cos(theta[0] - theta[1]);
  double d = m2*l2;

  /* HARD CODE MATRIX INVERSE */
  arma::mat inv_m(2,2);
  double det = 1/(a*d - b*c);
  
  inv_m(0,0) =  det*d;
  inv_m(0,1) = -det*b;
  inv_m(1,0) = -det*c;
  inv_m(1,1) =  det*a;
  
  arma::vec rhs(2);
  rhs[0] = -2*g*(m1+m2)*sin(theta[0]) - m2*l2*pow(theta_dot[1],2)*sin(theta[0] - theta[1]);
  rhs[1] = -g*m2*sin(theta[1]) + m2*l1*pow(theta_dot[0],2)*sin(theta[0] - theta[1]);

  /* SOLVE THE SYSTEM */
  theta_dot_dot = inv_m*rhs;  
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
  double dt = 1e-5;
  double tf = 100;
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
    //  V = -(m1+m2)*L1*g*np.cos(th1) - m2*L2*g*np.cos(th2)
    //  T = 0.5*m1*(L1*th1d)**2 + 0.5*m2*((L1*th1d)**2 + (L2*th2d)**2 +
    //          2*L1*L2*th1d*th2d*np.cos(th1-th2))
    
    V = -(m1+m2)*l1*g*cos(theta[0]) - m2*l2*g*cos(theta[1]);
    T = 0.5*m1*pow(l1*theta_dot[0],2) + 0.5*m2*(pow(l1*theta_dot[0],2) + pow(l2*theta_dot[1],2) + 2*l1*l2*theta_dot[0]*theta_dot[1]*cos(theta[0]-theta[1]));
    H = T - V;

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
