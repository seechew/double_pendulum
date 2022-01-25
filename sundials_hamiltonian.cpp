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
#include <chrono>
#include <fstream>

#include <arkode/arkode_erkstep.h>
#include <nvector/nvector_serial.h>    /* serial N_Vector types, fcts., macros */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype, etc */
#include <sundials/sundials_math.h>    /* seechew: needed for the exp function */

const double m1 = 10;
const double m2 = 10;
const double l1 = 10;
const double l2 = 10;
const double g = 9.8067;
const double pi = 3.14159265358979323846;

typedef boost::array<double,4> state_type;
std::ofstream csv("data_sundials_hamiltonian.csv");

/* USER-SPECIFIED FUNCTION CALLED BY THE SUNDIALS SOLVER */
static int myf(realtype t, N_Vector q, N_Vector qdot, void *user_data);

//------------------------------------------------------------------------
void myfunc(const state_type& q, state_type& qdot, double t) {

  /* UNPACK */
  double t1 = q[0];
  double t2 = q[1];
  double p1 = q[2];
  double p2 = q[3];
  
  double C1 = (p1*p2*sin(t1-t2)) / (l1*l2*(m1 + m2*pow(sin(t1-t2),2)));
  double C2 = pow(l2,2)*pow(m2,2)*pow(p1,2) + pow(l1,2)*(m1+m2)*pow(p2,2) - l1*l2*m2*p1*p2*cos(t1-t2);
  C2 = C2 * sin(2*(t1-t2));
  C2 = C2 / (2*pow(l1,2)*pow(l2,2)*pow((m1 + m2*pow(sin(t1-t2),2)),2));

  qdot[0] = (l2*p1 - l1*p2*cos(t1-t2)) / (pow(l1,2)*l2*(m1 + m2*pow(sin(t1-t2),2)));
  qdot[1] = (l1*(m1+m2)*p2 - l2*m2*p1*cos(t1-t2)) / (l1*pow(l2,2)*m2*(m1 + m2*pow(sin(t1-t2),2)));
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

  //  V = -(m1+m2)*L1*g*np.cos(th1) - m2*L2*g*np.cos(th2)
  //  T = 0.5*m1*(L1*th1d)**2 + 0.5*m2*((L1*th1d)**2 + (L2*th2d)**2 +
  //          2*L1*L2*th1d*th2d*np.cos(th1-th2))
  
  double V = -(m1+m2)*l1*g*cos(t1) - m2*l2*g*cos(t2);
  double T = 0.5*m1*pow(l1*p1,2) + 0.5*m2*(pow(l1*p1,2) + pow(l2*p2,2) + 2*l1*l2*p1*p2*cos(t1-t2));
  
  double H = T - V;
  
  double x1 =  l1*sin(t1);
  double y1 = -l1*cos(t1);
  double x2 =  x1 + l2*sin(t2);
  double y2 =  y1 - l2*cos(t2);
  
  csv << t << "," << x2 << "," << y2 << "," << H << std::endl;
  
}

//----------------------------------------------------------------------------
static int myf(realtype t, N_Vector q, N_Vector qdot, void *user_data)
{

  state_type q_in;
  state_type qdot_out;
  
  q_in[0]  = NV_Ith_S(q,0);
  q_in[1]  = NV_Ith_S(q,1);
  q_in[2]  = NV_Ith_S(q,2);
  q_in[3]  = NV_Ith_S(q,3);
   
  /* CALLING SIMPLE MBD */
  myfunc(q_in, qdot_out, t);

  /* REASSIGN OUTPUTS */
  NV_Ith_S(qdot,0)  = qdot_out[0];
  NV_Ith_S(qdot,1)  = qdot_out[1];
  NV_Ith_S(qdot,2)  = qdot_out[2];
  NV_Ith_S(qdot,3)  = qdot_out[3];
  
  return 0;
}

//------------------------------------------------------------------------
int main() {

  realtype t0 = RCONST(0.0);             /* INITIAL TIME */
  realtype tf= RCONST(1.0);              /* FINAL TIME */
  int n_steps = 7;                       /* NUMBER OF STEPS */
  sunindextype neq = 4;                  /* NUMBER OF EQUATION */
  realtype reltol = RCONST(1.0e-6);      /* TOLERANCES */
  realtype abstol = RCONST(1.0e-10);    

  void *erkode_mem = NULL;      /* EMPTY ERKODE MEMORY STRUCTURE */

  N_Vector q = N_VNew_Serial(neq);  /* for older SUNDIALS solver only */

  // IC
  NV_Ith_S(q,0)  = pi/2;         /* SET INITIAL CONDITION */
  NV_Ith_S(q,1)  = pi;           /* SET INITIAL CONDITION */
  NV_Ith_S(q,2)  = 0.0;          /* SET INITIAL CONDITION */
  NV_Ith_S(q,3)  = 0.0;          /* SET INITIAL CONDITION */

  erkode_mem = ERKStepCreate(myf, t0, q);    /* SET METHOD TO INTEGRATE */ 
   
  int flag = ERKStepSStolerances(erkode_mem, reltol, abstol);  /* SET TOLERANCES */

  //flag = ERKStepSetTableNum(erkode_mem,ARKODE_CASH_KARP_6_4_5);
  //flag = ERKStepSetTableNum(erkode_mem,ARKODE_BOGACKI_SHAMPINE_4_2_3);
  flag = ERKStepSetTableNum(erkode_mem,BOGACKI_SHAMPINE_4_2_3);























  // SET IC
  state_type q;
  std::fill(q.begin(),q.end(),0.0);
  q[0] = pi/2; 
  q[1] = pi;
  q[2] = 0; 
  q[3] = 0; 

  double t0 = 0.0;
  double tf = 100;

  csv << "time,x2,y2,H" << std::endl;
  
  /* START TIMING */
  auto start = std::chrono::steady_clock::now();
  
  /* COMPUTATIONAL CODE */
  
  auto end = std::chrono::steady_clock::now();

  std::cout << "DONE, ELAPSED TIME(s):  " << std::chrono::duration_cast<std::chrono::seconds>(end-start).count() << std::endl;

}
