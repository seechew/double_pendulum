/* g++ hamiltonian.cpp -o hamiltonian -larmadillo -std=c++17*/
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
  double t1 = q[0];
  double t2 = q[1];
  double p1 = q[2];
  double p2 = q[3];

  double s = sin(t1-t2);
  double c = cos(t1-t2);
  
  double C1 = (p1*p2*s) / (l1*l2*(m1 + m2*pow(s,2)));
  double C2 = pow(l2*p1,2)*m2 + pow(l1*p2,2)*(m1+m2) - l1*l2*m2*p1*p2*c;
  C2 = C2 * sin(2*(t1-t2));
  C2 = C2 / (2*pow(l1*l2,2)*pow((m1 + m2*pow(s,2)),2));

  qdot[0] = (l2*p1 - l1*p2*c) / (pow(l1,2)*l2*(m1 + m2*pow(s,2)));
  qdot[1] = (l1*(m1+m2)*p2 - l2*m2*p1*c) / (l1*pow(l2,2)*m2*(m1 + m2*pow(s,2)));
  qdot[2] = -(m1+m2)*g*l1*sin(t1) - C1 + C2;
  qdot[3] = -m2*g*l2*sin(t2) + C1 - C2;
  
}

//----------------------------------------------------------------------------
void write_results( const state_type &q, const double t)
{
  /* PRINT HAMILTONIAN */
  
  double t1     = q[0];
  double t2     = q[1];
  double p1     = q[2];
  double p2     = q[3];

  double td1    = (l2*p1 - l1*p2*cos(t1-t2))/(pow(l1,2)*l2*(m1+m2*pow(sin(t1-t2),2)));
  double td2    = (-m2*l2*p1*cos(t1-t2)+(m1+m2)*l1*p2)/(m2*l1*pow(l2,2)*(m1+m2*pow(sin(t1-t2),2)));
  
  //  V = -(m1+m2)*L1*g*np.cos(th1) - m2*L2*g*np.cos(th2)
  //  T = 0.5*m1*(L1*th1d)**2 + 0.5*m2*((L1*th1d)**2 + (L2*th2d)**2 +
  //          2*L1*L2*th1d*th2d*np.cos(th1-th2))
  
  double V = -(m1+m2)*l1*g*cos(t1) - m2*l2*g*cos(t2);
  double T = 0.5*m1*pow(l1*td1,2) + 0.5*m2*(pow(l1*td1,2) + pow(l2*td2,2) + 2*l1*l2*td1*td2*cos(t1-t2));
  
  double H = T + V;

  double x1 =  l1*sin(t1);
  double y1 = -l1*cos(t1);
  double x2 =  x1 + l2*sin(t2);
  double y2 =  y1 - l2*cos(t2);
  
  csv << t << "," << x2 << "," << y2 << "," << H << std::endl;
  
}

//------------------------------------------------------------------------
int main() {

  // SET IC
  state_type q;
  std::fill(q.begin(),q.end(),0.0);
  q[0] = pi/2; 
  q[1] = pi;
  q[2] = 0; 
  q[3] = 0; 

  double t0 = 0.0;
  double tf = 100;
  double dt = 1e-2;

  csv << "time,x2,y2,H" << std::endl;
  
  /* START TIMING */
  auto start = std::chrono::steady_clock::now();
  
  //boost::numeric::odeint::integrate(myfunc, q, t0, tf, 1e-5, write_results);
  boost::numeric::odeint::integrate_const(boost::numeric::odeint::runge_kutta_fehlberg78<state_type>(), myfunc, q, t0, tf, dt, write_results);
  //boost::numeric::odeint::integrate_n_steps(boost::numeric::odeint::euler<state_type>(),simple_mbd,q,0.0,0.14285714,7,write_simple_mbd); 
  //boost::numeric::odeint::integrate_const(boost::numeric::odeint::runge_kutta_dopri5<state_type>(),myfunc,q,t0,tf,1e-5,write_results);
  //boost::numeric::odeint::integrate_const(boost::numeric::odeint::symplectic_rkn_sb3a_mclachlan<state_type>(),myfunc,q,0.0,7200.0,1e-5,write_results);
  //boost::numeric::odeint::integrate_adaptive(boost::numeric::odeint::runge_kutta_dopri5<state_type>(),simple_mbd,q,0.0,1.0,0.0001,write_simple_mbd);
  //boost::numeric::odeint::integrate_adaptive(boost::numeric::odeint::runge_kutta_fehlberg78<state_type>(),myfunc,q,0.0,36000.0,1e-5,write_results);

  //std::cout << std::endl; // END PRINTING
  
  auto end = std::chrono::steady_clock::now();

  std::cout << "DONE, ELAPSED TIME(s):  " << std::chrono::duration_cast<std::chrono::seconds>(end-start).count() << std::endl;

}
