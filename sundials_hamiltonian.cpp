/* g++ sundials_hamiltonian.cpp -O3 -DNDEBUG  -I/usr/local/include -L/usr/local/lib -larmadillo -lsundials_arkode -lsundials_nvecserial -lsundials_nvecmanyvector -lm -Wl,-rpath,/usr/local/lib -o sundials.out

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
void write_results( const N_Vector q, const double t)
{
  /* PRINT HAMILTONIAN */
  
  double th1     = NV_Ith_S(q,0);
  double th1d    = NV_Ith_S(q,1);
  double th2     = NV_Ith_S(q,2);
  double th2d    = NV_Ith_S(q,3);

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
  realtype tf = RCONST(7200.0);          /* FINAL TIME */
  realtype expected_dt = RCONST(1e-2);   /* */
  int n_steps;                           /* NUMBER OF STEPS */
  sunindextype neq = 4;                  /* NUMBER OF EQUATION */
  realtype reltol = RCONST(1.0e-6);      /* TOLERANCES */
  realtype abstol = RCONST(1.0e-10);    

  void *erkode_mem = NULL;      /* EMPTY ERKODE MEMORY STRUCTURE */

  // N_Vector q = N_VNew_Serial(neq);  /* for older SUNDIALS solver only */
  sundials::Context context;
  N_Vector q = N_VNew_Serial(neq, context);

  // IC
  NV_Ith_S(q,0)  = pi/2;         /* SET INITIAL CONDITION */
  NV_Ith_S(q,1)  = 0.0;          /* SET INITIAL CONDITION */
  NV_Ith_S(q,2)  = pi;           /* SET INITIAL CONDITION */
  NV_Ith_S(q,3)  = 0.0;          /* SET INITIAL CONDITION */

  // erkode_mem = ERKStepCreate(myf, t0, q);    /* SET METHOD TO INTEGRATE */  // for older SUNDIALS solver only
  erkode_mem = ERKStepCreate(myf, t0, q, context);    /* SET METHOD TO INTEGRATE */

  int flag = ERKStepSStolerances(erkode_mem, reltol, abstol);  /* SET TOLERANCES */

  /* flag = ERKStepSetTableNum(erkode_mem,BOGACKI_SHAMPINE_4_2_3);  */ // for older SUNDIALS solver only 

  //flag = ERKStepSetTableNum(erkode_mem,ARKODE_BOGACKI_SHAMPINE_4_2_3);
  //flag = ERKStepSetTableNum(erkode_mem,ARKODE_CASH_KARP_6_4_5);
  //flag = ERKStepSetTableNum(erkode_mem,ARKODE_DORMAND_PRINCE_7_4_5);
  flag = ERKStepSetTableNum(erkode_mem,ARKODE_ARK548L2SAb_ERK_8_4_5);
  //flag = ERKStepSetTableNum(erkode_mem,ARKODE_VERNER_8_5_6);
  //flag = ERKStepSetTableNum(erkode_mem,ARKODE_FEHLBERG_13_7_8);

  realtype t = t0;
  //realtype dTout = (tf - t0)/n_steps;
  n_steps = (int) (tf - t0)/expected_dt;
  realtype dTout = expected_dt;
  realtype tout = t0 + dTout;

  csv << "time,x2,y2,H" << std::endl;

  /* START TIMING */
  auto start = std::chrono::steady_clock::now();

  /* COMPUTATIONAL CODE */
  for (int iout=0; iout<n_steps; iout++)
    {

      flag = ERKStepEvolve(erkode_mem, tout, q, &t, ARK_NORMAL);

      /* CALL FUNCTION TO PRINT OUTPUT TO SCREEN */
      write_results(q,t);
      
      /* CHECK INTEGRATOR FLAG */
      if (flag >= 0) {     /* SUCCESSFUL SOLVE: UPDATE OUTPUT TIME */
	tout += dTout;
	tout = (tout > tf) ? tf : tout;
      } else {          /* UNSUCCESSFUL SOLVE: BREAK */
	std::cerr << "+error: integrator" << std::endl;
	break;
      }
    }
  
  auto end = std::chrono::steady_clock::now();

  long int nst, nst_a, nfe, netf;
  flag = ERKStepGetNumSteps(erkode_mem, &nst);
  flag = ERKStepGetNumStepAttempts(erkode_mem, &nst_a);
  flag = ERKStepGetNumRhsEvals(erkode_mem, &nfe);
  flag = ERKStepGetNumErrTestFails(erkode_mem, &netf);

  printf("\nFinal Solver Statistics:\n");
  printf("   Internal solver steps = %li (attempted = %li)\n", nst, nst_a);
  printf("   Total RHS evals = %li\n", nfe);
  printf("   Total number of error test failures = %li\n\n", netf);
  
  std::cout << "DONE, ELAPSED TIME(s):  " << std::chrono::duration_cast<std::chrono::seconds>(end-start).count() << std::endl;

}
