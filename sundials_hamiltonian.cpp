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
void write_results( const N_Vector q, const double t)
{
  /* PRINT HAMILTONIAN */
  
  double t1     = NV_Ith_S(q,0);
  double t2     = NV_Ith_S(q,1);
  double p1     = NV_Ith_S(q,2);
  double p2     = NV_Ith_S(q,3);

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
  realtype tf= RCONST(100.0);           /* FINAL TIME */
  realtype expected_dt = RCONST(1e-2);    /* */
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
  NV_Ith_S(q,1)  = pi;           /* SET INITIAL CONDITION */
  NV_Ith_S(q,2)  = 0.0;          /* SET INITIAL CONDITION */
  NV_Ith_S(q,3)  = 0.0;          /* SET INITIAL CONDITION */

  // erkode_mem = ERKStepCreate(myf, t0, q);    /* SET METHOD TO INTEGRATE */  // for older SUNDIALS solver only
  erkode_mem = ERKStepCreate(myf, t0, q, context);    /* SET METHOD TO INTEGRATE */

  int flag = ERKStepSStolerances(erkode_mem, reltol, abstol);  /* SET TOLERANCES */

  /* flag = ERKStepSetTableNum(erkode_mem,BOGACKI_SHAMPINE_4_2_3);  */ // for older SUNDIALS solver only 
  
  //flag = ERKStepSetTableNum(erkode_mem,ARKODE_CASH_KARP_6_4_5);
  //flag = ERKStepSetTableNum(erkode_mem,ARKODE_DORMAND_PRINCE_7_4_5);
  flag = ERKStepSetTableNum(erkode_mem,ARKODE_FEHLBERG_13_7_8);
  //flag = ERKStepSetTableNum(erkode_mem,ARKODE_BOGACKI_SHAMPINE_4_2_3);

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
