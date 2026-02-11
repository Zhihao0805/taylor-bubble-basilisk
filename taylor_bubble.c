#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "maxruntime.h"

/** 
The solution can be calculated with the bubble rising in a fixed tank or we can follow the bubble as it rises. 
In the second case, the equations are solved in the frame of reference of the bubble (inertial frame of reference), 
the domain can be a little smaller and the simulations can run for longer as the bubble does not hit a wall or escape the domain.
*/
#define bubble_frame 1

#define RHO21 (0.1)
#define MU21 (0.1)
#define VelErr (1e-2)      // error tolerances in velocity
#define tsnap (0.01)

int MAXlevel;
int tint = 0;
double Oh, Bond, tmax;
//double Ma, betas, cs_init, gamma_inf;
//double Pe_s = 100.;
const double MAXTIME = 200.;

// Parameters for solving the equations in the frame moving with the bubble
double ub = 0.;
double ubold = 0.;
double sbold = 0.;
double dubdt = 0.;
double xframe = 0.;

//
char nameOut[80], dumpFile[80], restartFile[80];
//

#if bubble_frame == 1
  const double WIDTH = 8. ; 
  const double Xi = 1.505, Yi = 0, Zi = 0.;
  p[left] = dirichlet(0.);
  pf[left] = dirichlet(0.);
  u.t[right] = dirichlet(0.0);
  u.n[right] = dirichlet(-ub);
  u.n[left] = neumann(0.0);
//  uf.n[left] = neumann(0.0);
#else
  const double WIDTH = 128. ;
  const double Xi = 3.5, Yi = 0, Zi = 0.;
  p[right] = dirichlet(0.);
  pf[right] = dirichlet(0.);
#endif

u.n[top] = dirichlet(0.0);

/**
The main function can take two optional parameters: the maximum level
of adaptive refinement (as well as an optional maximum runtime). */

int main (int argc, char * argv[]) {
 
  /**
  We set the domain geometry and initial refinement. */
  
  size (WIDTH);
  origin (-L0/2, 0, 0);
  init_grid (16);

    MAXlevel = atoi(argv[1]);
//  tauy = atof(argv[2]);
  Bond = atof(argv[2]);
  Oh = atof(argv[3]);
  Ma = atof(argv[4]);
  cs_init = atof(argv[5]);

  // Ensure that all the variables were transferred properly from the terminal or job script.
  if (argc < 5){
    fprintf(ferr, "Lack of command line arguments. Check! Need %d more arguments\n",6-argc);
    return 1;
  }
  fprintf(ferr, "Level %d, Oh %2.1e, Bo %4.3f, Ma %2.1e, cs_init %2.1e\n", MAXlevel, Oh, Bond, Ma, cs_init);

  


