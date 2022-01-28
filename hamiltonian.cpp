/* g++ hamiltonian.cpp -o hamiltonian.out -larmadillo -std=c++17*/
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
std::ofstream csv("data_hamiltonian.csv");

//------------------------------------------------------------------------
void myfunc(const state_type& q, state_type& qdot, double t) {

  /* UNPACK */
  double theta1 = q[0];
  double z1     = q[1];
  double theta2 = q[2];
  double z2     = q[3];

  double s = sin(theta1-theta2);
  double c = cos(theta1-theta2);

  double theta1dot = z1;
  double z1dot = (m2*g*sin(theta2)*c - m2*s*(l1*pow(z1,2)*c + l2*pow(z2,2)) -
		 (m1+m2)*g*sin(theta1)) / l1 / (m1 + m2*pow(s,2));
  double theta2dot = z2;
  double z2dot = ((m1+m2)*(l1*pow(z1,2)*s - g*sin(theta2) + g*sin(theta1)*c) + 
		 m2*l2*pow(z2,2)*s*c) / l2 / (m1 + m2*pow(s,2));
  
  qdot[0] = theta1dot;
  qdot[1] = z1dot;
  qdot[2] = theta2dot;
  qdot[3] = z2dot;

}

//----------------------------------------------------------------------------
void write_results( const state_type &q, const double t)
{
  /* PRINT HAMILTONIAN */  
  double th1     = q[0];
  double th1d    = q[1];
  double th2     = q[2];
  double th2d    = q[3];

  double V = -(m1+m2)*l1*g*cos(th1) - m2*l2*g*cos(th2);
    
  double T = 0.5*m1*pow(l1*th1d,2) + 0.5*m2*(pow(l1*th1d,2) + pow(l2*th2d,2) +
	     2*l1*l2*th1d*th2d*cos(th1-th2));
  
  double H = T + V;
  
  double x1 =  l1*sin(th1);
  double y1 = -l1*cos(th1);
  double x2 =  x1 + l2*sin(th2);
  double y2 =  y1 - l2*cos(th2);

  csv << t << "," << x2 << "," << y2 << "," << H << std::endl;
  
}

//------------------------------------------------------------------------
int main() {

  // SET IC
  state_type q;
  std::fill(q.begin(),q.end(),0.0);
  q[0] = pi/2; 
  q[1] = 0;
  q[2] = pi; 
  q[3] = 0; 

  double t0 = 0.0;
  double tf = 7200;
  double dt = 1e-2;

  csv << "time,x2,y2,H" << std::endl;
  
  /* START TIMING */
  auto start = std::chrono::steady_clock::now();
  
  //boost::numeric::odeint::integrate(myfunc, q, t0, tf, dt, write_results);
  boost::numeric::odeint::integrate_const(boost::numeric::odeint::runge_kutta_fehlberg78<state_type>(), myfunc, q, t0, tf, dt, write_results);
  //boost::numeric::odeint::integrate_const(boost::numeric::odeint::runge_kutta_fehlberg78<state_type>(), myfunc, q, t0, tf, dt, write_results);
  //boost::numeric::odeint::integrate_n_steps(boost::numeric::odeint::euler<state_type>(),simple_mbd,q,0.0,0.14285714,7,write_simple_mbd); 
  //boost::numeric::odeint::integrate_const(boost::numeric::odeint::runge_kutta_dopri5<state_type>(),myfunc,q,t0,tf,1e-5,write_results);
  //boost::numeric::odeint::integrate_const(boost::numeric::odeint::symplectic_rkn_sb3a_mclachlan<state_type>(),myfunc,q,0.0,7200.0,1e-5,write_results);
  //boost::numeric::odeint::integrate_adaptive(boost::numeric::odeint::runge_kutta_dopri5<state_type>(),simple_mbd,q,0.0,1.0,0.0001,write_simple_mbd);
  //boost::numeric::odeint::integrate_adaptive(boost::numeric::odeint::runge_kutta_fehlberg78<state_type>(),myfunc,q,0.0,36000.0,1e-5,write_results);

  //std::cout << std::endl; // END PRINTING
  
  auto end = std::chrono::steady_clock::now();

  std::cout << "DONE, ELAPSED TIME(s):  " << std::chrono::duration_cast<std::chrono::seconds>(end-start).count() << std::endl;

}
