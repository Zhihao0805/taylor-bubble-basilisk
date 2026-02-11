/**
# Bubble rising in a large tank

We wish to study the behaviour of a single bubble rising "in a large
tank" i.e. far from any boundaries.

We use the centered Navier--Stokes solver and log performance
statistics. */
#include "axi.h"
#include "navier-stokes/centered.h"
//#include "navier-stokes/perfs.h"

/**
We have two phases e.g. air and water. For large viscosity and density
ratios, the harmonic mean for the viscosity tends to work better than
the default arithmetic mean. We "overload" the default by defining the
*mu()* macro before including the code for two phases. */

//#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
//#define FILTERED 1
#include "two-phase.h"

#define sigmavar(csi)  (1. + betas*log(1.-min(csi,0.95*gamma_inf)/gamma_inf))
scalar distLS[];

#include "../../src/aslam.h"
#include "../../src/tension_loc2.h"
#include "../../src/surfactant.h"

/**
We also need surface tension, and in 3D only we will use the
$\lambda_2$ criterion of [Jeong and Hussain,
1995](/src/references.bib#jeong1995) to display the vortices using
Basilisk View. */


/** 
The solution can be calculated with the bubble rising in a fixed tank or we can follow the bubble as it rises. 
In the second case, the equations are solved in the frame of reference of the bubble (inertial frame of reference), 
the domain can be a little smaller and the simulations can run for longer as the bubble does not hit a wall or escape the domain.
*/
#define bubble_frame 1


/**
We can control the maximum runtime. */

#include "maxruntime.h"
#include "../../src/output_vtu.h"

/**
The density ratio is 1000 and the dynamic viscosity ratio 100. */


/**
We try to replicate the results of [Cano-Lozano et al,
2016](/src/references.bib#cano2016) (obtained with Gerris). Aside from
the ratios above, there are two independent parameters which can be
described by the Galilei number
$$
Ga^2 = \frac{g D^3}{\nu^2}
$$
with $g$ the acceleration of gravity, $D$ the diameter of the bubble
and $\nu$ the kinematic viscosity of the outer fluid; and the
[Bond/Eötvös](https://en.wikipedia.org/wiki/E%C3%B6tv%C3%B6s_number)
number
$$
Bo = \frac{\rho g D^2}{\sigma}
$$
with $\rho$ the density of the outer fluid and $\sigma$ the surface
tension coefficient.

We consider two bubbles studied by Cano-Lozano et al, 2016. */

// Bubble 19 of Cano-Lozano et al, P.R.Fluids, 2016
#define RHO21 (1e-2)
#define MU21 (2e-2)

#define VelErr (1e-2)      // error tolerances in velocity

#define tsnap (0.01)



// const double Ga = 100.;
// const double Bo = 4.;
int MAXlevel;
int tint = 0;
double Oh, Bond, tmax;
double Ma, betas, cs_init, gamma_inf;
double Pe_s = 100.;
const double MAXTIME = 200.;


// Frequency of Paraview outputs (dtoutput) and frequency of "dumping" the solution for restart files (dtdump)
double dtoutput = 1.0;
double dtdump = 5.;

// Parameters for solving the equations in the frame moving with the bubble
double ub = 0.;
double ubold = 0.;
double sbold = 0.;
double dubdt = 0.;
double xframe = 0.;

//
char nameOut[80], dumpFile[80], restartFile[80];
//

/**
We choose as length unit the diameter of the bubble. The domain is
$120^3$. *Zi* is the initial position of the bubble relative to the
bottom wall. The acceleration of gravity is set to unity, which gives
a characteristic rise velocity also of order unity, which gives a
maximum time for the simulation comparable to the domain size. */

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


  /**
  We set the physical parameters: densities, viscosities and surface
  tension. */

    rho1 = 1., rho2 = RHO21;
  mu1 = Oh, mu2 = MU21*Oh;
  f.sigma0 = 1.0;
  double Uc = sqrt(Bond); //sqrt(g*d)
  double Lc = 1.;
  betas = mu1*Uc*Ma/f.sigma0;
  gamma_inf = 1.0;
  Ds = Uc*Lc/Pe_s;


  /**
  We reduce the tolerance on the divergence of the flow. This is
  important to minimise mass conservation errors for these simulations
  which are very long. */
  
  DT = 1.e-3;

  TOLERANCE = 1e-8 [*];
  run();
}

/**
For the initial conditions, we first try to restore the simulation
from a previous "restart", if this fails we refine the mesh locally to
the maximum level, in a sphere of diameter 1.5 around the bubble. We
then initialise the volume fraction for a bubble initially at (0,Yi,0)
of diameter unity. */

event init (t = 0) {
  if (!restore (file = "restart")) {
    refine (sq(x - Xi) + sq(y - Yi) + sq(z - Zi) - sq(1.) < 0 && level < MAXlevel);
    vertex scalar phi0[];
    vertex scalar phi_temp[];
     foreach_vertex() {
        phi0[] = 9999999.0;
     }   

    foreach_vertex() {
        double xc = Xi;
        double yc = Yi; 
        double zc = Zi;
        double radius = 0.5 [1];
        phi_temp[] = (sqrt(sq(x-xc) + sq(y - yc) + sq(z - zc)) - (radius));
        if (sq(phi_temp[]) < sq(phi0[])) {
            phi0[] = phi_temp[];
        }
    }
    scalar f0[];
    fractions (phi0,f0);
    f0.refine = f0.prolongation = fraction_refine;
    restriction ({f0});

    foreach() {
      f[] = f0[];
      cs_int[] = cs_init;
    }
  }
}

event properties (i++) {
    foreach() {
      sigmaloc[] = max(0.05,sigmavar(cs_int[])); //max(0.05,sigmavar(cs_int[]));
  }
}

/**
We add the acceleration of gravity (unity) in the downward (-y)
direction. */

event acceleration (i++) {
  face vector av = a;
  double sigma0 = f.sigma0;

  foreach_face(x) {
    av.x[] -= (dubdt + Bond);
  }
  scalar gradf[];
  foreach() {
    gradf[] = 0.;
    foreach_dimension() {
      gradf[] += sq((f[1]-f[-1])/(2.*Delta));
    }
    gradf[] = sqrt(gradf[]);
  }

  foreach_face(x) {
    av.x[] += 2./(rho1 + rho2)*sigma0*((sigmaloc[0,0]-sigmaloc[-1,0])/Delta)*0.5*(gradf[0,0] + gradf[-1,0]); 
  }
  foreach_face(y) {
    av.y[] += 2./(rho1 + rho2)*sigma0*((sigmaloc[0,0]-sigmaloc[0,-1])/Delta)*0.5*(gradf[] + gradf[0,-1]); 
  }

}

/**
We adapt the mesh by controlling the error on the volume fraction and
velocity field. */

event adapt (i++) {
  double uemax = VelErr;
  scalar ffilt[];
  foreach() {
    ffilt[] = pow((f[1,0]-f[-1,0])/(2.*Delta),2.) + pow((f[0,1]-f[0,-1])/(2.*Delta),2.);
  }
  adapt_wavelet ({ffilt,u}, (double[]){0.01,uemax,uemax}, MAXlevel, 3);
}

/**
## Outputs

Every ten timesteps, we output the time, volume, position, and
velocity of the bubble. */

event logfile (i++) {
#if bubble_frame == 1
  ubold = ub;
#else
  ubold = 0.;
#endif
  double xb = 0., yb = 0., zb = 0., sb = 0.;
  double vbx = 0., vby = 0., vbz = 0.;
  foreach(reduction(+:xb) reduction(+:yb) reduction(+:zb)
	  reduction(+:vbx) reduction(+:vby) reduction(+:vbz)
	  reduction(+:sb)) {
    double dv = (1.-f[])*dv();
    xb += x*dv;
    yb += y*dv;
    zb += z*dv;
    vbx += u.x[]*dv;
    vby += u.y[]*dv;
    vbz += u.z[]*dv;
    sb += dv;
  }
  ub = ubold + vbx/sb;
  sbold = sb;
#if bubble_frame == 1
  dubdt = (ub - ubold)/dt;
#else
  dubdt = 0.;
#endif
  xframe += ubold*dt;
  if ((i%10) == 0) {
    fprintf (stderr,
	   "%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n", 
	   t, sb,
	   xb/sb + xframe, yb/sb, zb/sb,
	   ub, vby/sb, vbz/sb);
    fflush (stderr);
  }
}

/**
Every time unit, we output a full snapshot of the simulation, to be
able to restart and for visualisation. In three dimensions, we compute
the value of the $\lambda_2$ field which will be used for
visualisation of vortices, as well as the streamwise vorticity
$\omega_y = \partial_x u_z - \partial_z u_x$. */

/* event snapshot (t = 0; t+=dtoutput; t <= MAXTIME)
{
  
  char fname[80]; 
  sprintf(fname,"fields_%04d.vtk",(int) tint); 
  FILE * fp = fopen(fname,"w"); 
  output_vtu((scalar*){f,p,u,cs_int},(vector*){0},fname); 
  tint += 1;

  fclose(fp); 
}

event dump_snapshot (t = dtdump; t+=dtdump; t <= MAXTIME)
{
  p.nodump = false;

  int tint = t/dtdump;
  char name[80];
  sprintf (name, "dump-%03d", (int) tint);
  dump (file = name);

}
 */

event writingFiles (t = 0; t += tsnap; t <= MAXTIME) {
  dump (file = dumpFile);
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump(file=nameOut);
}

/* char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm); */

