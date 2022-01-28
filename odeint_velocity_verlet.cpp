/* g++ odeint_velocity_verlet.cpp -o odeint_vv.out */
/*
   Reference: https://scienceworld.wolfram.com/physics/DoublePendulum.html   
*/

#include <iostream>
#include <cmath>
#include <vector>
#include <armadillo>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <iomanip>  //std::setprecision
#include <chrono>
#include <fstream>

const double m1 = 10;
const double m2 = 10;
const double l1 = 10;
const double l2 = 10;
const double g = 9.8067;
const double pi = 3.14159265358979323846;

typedef boost::array<double,4> state_type;
typedef boost::array<double,2> vector_type;

std::ofstream csv("data_odeint_vv.csv");

//------------------------------------------------------------------------
void myfunc(const vector_type& q, const vector_type& qdot, vector_type& qdotdot, double t) {

  
  /* UNPACK */
  double theta1 = q[0];
  double z1     = qdot[0];
  double theta2 = q[1];
  double z2     = qdot[1];

  double s = sin(theta1-theta2);
  double c = cos(theta1-theta2);

  double z1dot = (m2*g*sin(theta2)*c - m2*s*(l1*pow(z1,2)*c + l2*pow(z2,2)) -
		 (m1+m2)*g*sin(theta1)) / l1 / (m1 + m2*pow(s,2));
  double z2dot = ((m1+m2)*(l1*pow(z1,2)*s - g*sin(theta2) + g*sin(theta1)*c) + 
		 m2*l2*pow(z2,2)*s*c) / l2 / (m1 + m2*pow(s,2));
  
  qdotdot[0] = z1dot;
  qdotdot[1] = z2dot;

}

//----------------------------------------------------------------------------
void write_results(const vector_type &q, const double t)
{
  /* PRINT HAMILTONIAN */  
  double th1     = q[0];
  double th2     = q[1];
  
  double x1 =  l1*sin(th1);
  double y1 = -l1*cos(th1);
  double x2 =  x1 + l2*sin(th2);
  double y2 =  y1 - l2*cos(th2);

  csv << t << "," << x2 << "," << y2 << std::endl;
  
}

//------------------------------------------------------------------------
int main() {

  // SET IC
  std::pair< vector_type , vector_type > q;
  q.first[0]  = pi/2;
  q.first[1]  = pi;
  q.second[0] = 0.0;
  q.second[1] = 0.0;

  double t0 = 0.0;
  double tf = 7200;
  double dt = 1e-2;
  double current_t = t0 + dt; // INCREMENT 1
  double V, T, H;
  double x1, x2, y1, y2;
  
  boost::numeric::odeint::velocity_verlet<vector_type> stepper;
  
  csv << "time,x2,y2,H" << std::endl;
  
  /* START TIMING */
  auto start = std::chrono::steady_clock::now();

  while (current_t < tf) {
    
    stepper.do_step(myfunc,q,current_t,dt);

    /* PRINT HAMILTONIAN */
    double th1  = q.first[0];
    double th2  = q.first[1];
    double th1d = q.second[0];
    double th2d = q.second[1];
    
    V = -(m1+m2)*l1*g*cos(th1) - m2*l2*g*cos(th2);
    
    T = 0.5*m1*pow(l1*th1d,2) + 0.5*m2*(pow(l1*th1d,2) + pow(l2*th2d,2) +
	2*l1*l2*th1d*th2d*cos(th1-th2));
    
    H = T + V;
    
    x1 =  l1*sin(th1);
    y1 = -l1*cos(th1);
    x2 =  x1 + l2*sin(th2);
    y2 =  y1 - l2*cos(th2);
  
    csv << current_t << "," << x2 << "," << y2 << "," << H << std::endl;

    current_t += dt;

  }
  
  auto end = std::chrono::steady_clock::now();

  std::cout << "DONE, ELAPSED TIME(s):  " << std::chrono::duration_cast<std::chrono::seconds>(end-start).count() << std::endl;

}
