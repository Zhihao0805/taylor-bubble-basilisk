/**
 * taylor_clean_benchmark_axi_steady.c
 *
 * Clean Taylor bubble in a circular tube (axisymmetric r-z), using Basilisk:
 * - 2-phase incompressible Navier–Stokes (VOF)
 * - surface tension
 * - embedded boundary for tube wall (y = Rtube)
 * - periodic in axial direction (x)
 * - uniform grid (no AMR)
 *
 * Key fix for: embed-tree.h Assertion `coarse(cs)' failed
 *   -> Make wall level-set phi global and rebuild cs/fs after each adapt.
 */

#include "grid/quadtree.h"
#include "axi.h"
#undef SEPS
#define SEPS 1e-30
#include "embed.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "output.h"

/* ---------------- user / run-time parameters ---------------- */

static int    MAXLEVEL = 9;
static int    MINLEVEL = 5;
static double Re = 1.0;
static double Ca = 0.016;
static double Lz_over_R = 40.0;
static double t_end = 2.0;
static double t_steady_start = 0.0;

/* geometry */
static double Rtube = 1.0;      // tube radius (nondimensional)
static double Lx;               // domain length in x (axial)

/* bubble initialization */
static double xcm0 = 0.5;       // initial bubble centroid (axial)
static double Lb_over_R = 2.0;  // initial bubble length (in units of R), simple seed
static double cap_over_R = 1.0; // cap radius (in units of R), simple seed

/* steady statistics */
static double sum_h = 0., sum_Ub = 0., sum_xcm = 0., sum_Vgas = 0.;
static int n_steady = 0;
static double Vgas_ref = -1.;
static double xcm_unwrap = 0.;

/* ---------------- embedded boundary level-set (GLOBAL!) ---------------- */
vertex scalar phi[];   // tube wall level-set, required for AMR rebuild

/* ---------------- small helpers ---------------- */

static inline double taylor_law (double Ca_) {
  // Bretherton/Taylor-style correlation (used here as a reference)
  double c23 = pow(Ca_, 2./3.);
  return (1.34*c23)/(1. + 1.34*2.5*c23);
}

/* compute h_inf from maximum gas radius anywhere in the domain:
 * h/R = 1 - r_gas_max/R
 * This is robust when the bubble crosses periodic boundaries.
 */
static double film_thickness_global (void) {
  double rmax_g = 0.;
  foreach(reduction(max:rmax_g)) {
    if (cs[] > 0. && 1. - f[] > 0.5) {
      if (y > rmax_g) rmax_g = y;
    }
  }
  if (rmax_g <= 0.) return Rtube;
  double h = Rtube - rmax_g;
  return h;
}

/* ---------------- boundary conditions ---------------- */

/* periodic in x (axial) will be set in main() */

/* embedded wall no-slip (tube wall y=Rtube) */
u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);

/* no penetration through symmetry axis (y=0) is handled by axi.h */

/* ---------------- main ---------------- */

int main (int argc, char ** argv)
{
  if (argc > 1) MAXLEVEL   = atoi(argv[1]);
  if (argc > 2) Re         = atof(argv[2]);
  if (argc > 3) Ca         = atof(argv[3]);
  if (argc > 4) Lz_over_R  = atof(argv[4]);
  if (argc > 5) t_end      = atof(argv[5]);

  Lx = Lz_over_R*Rtube;
  t_steady_start = 0.7*t_end;
  fprintf (stderr, "[args] MAXLEVEL=%d Re=%g Ca=%g Lz/R=%g t_end=%g\n",
           MAXLEVEL, Re, Ca, Lz_over_R, t_end);

  /* domain */
  size (Lx);
  origin (0., 0.);         // axis at y=0
  N = 1 << MAXLEVEL;
  periodic (right);

  /* embedded fractions use default refine/prolongation from embed.h */
  CFL = 0.5;
  DT = 1e-3;


  /* physical properties (nondimensional) */
  // Reference choices:
  //  Uref = 1, Rref = Rtube = 1
  //  Re = rho*U*R/mu  => mu1 = 1/Re (if rho1=1)
  //  Ca = mu*U/sigma  => sigma = mu1/Ca
  rho1 = 1.0;          // liquid
  rho2 = 0.01;         // gas (you can tune)
  mu1  = 1.0/Re;
  mu2  = 0.01*mu1;     // gas viscosity ratio (tune if needed)
  f.sigma = mu1/Ca;

  fprintf (stderr,
           "ARGS: MAXLEVEL=%d Re=%g Ca=%g Lz/R=%g TMAX=%g\n",
           MAXLEVEL, Re, Ca, Lz_over_R, t_end);

  fprintf (stderr, "[main] before run: iter=%d t=%g dt=%g t_end=%g\n", iter, t, dt, t_end);
  run();
  fprintf (stderr, "[main] after run: iter=%d t=%g dt=%g t_end=%g\n", iter, t, dt, t_end);
}

/* ---------------- geometry: tube wall (only needs to set phi once) ---------------- */

event geometry (i = 0)
{
  foreach_vertex()
    phi[] = Rtube - y;  // phi>0 fluid region, phi<0 solid (outside tube)
  boundary ({phi});

  fractions (phi, cs, fs);
  restriction ({cs, fs});
  boundary ({cs, fs});
  return 0;
}

/* ---------------- initialization ---------------- */

event init (t = 0)
{
  // periodic already set in main
  // ensure wall fractions exist
  if (!cs.boundary) boundary({cs, fs});

  // simple seed: a "capsule" bubble centered at xcm0, radius ~ Rtube
  // gas region where f=0; liquid where f=1
  double Lb = Lb_over_R*Rtube;
  double Rc = cap_over_R*Rtube;

  foreach() {
    // inside tube only
    if (cs[] > 0.) {
      double dx = x - xcm0;
      // capsule: cylinder of half-length (Lb/2 - Rc) + hemispherical caps
      double half_cyl = max(0., 0.5*Lb - Rc);

      if (fabs(dx) <= half_cyl) {
        // cylinder part
      } else {
        // cap part (circle in r-z)
      }

      // make a central bubble by cutting radius slightly smaller than wall
      // so initial film is not exactly zero:
      double r_bub = 0.95*Rtube;
      double inside = (fabs(dx) <= 0.5*Lb) && (y <= r_bub);

      f[] = inside ? 0.0 : 1.0; // gas:0, liquid:1
      u.x[] = 0.;
      u.y[] = 0.;
    } else {
      // solid region
      f[] = 1.0;
      u.x[] = u.y[] = 0.;
    }
  }
  boundary ({f, u.x, u.y});

  // optional: quick init-check (gas volume + centroid)
  double Vgas = 0., xmom = 0.;
  foreach(reduction(+:Vgas) reduction(+:xmom)) {
    double g = (1. - f[])*cs[];
    if (g > 0.) {
      double dV = g*dv();
      Vgas += dV;
      xmom += x*dV;
    }
  }
  double xcm = (Vgas > 0 ? xmom/Vgas : 0.);
  fprintf (stderr, "[init-check] Vgas=%g xcm=%g ok=%d\n", Vgas, xcm, (Vgas > 0));
  return 0;
}

/* ---------------- diagnostics ---------------- */

event diagnostics (i++; t <= t_end)
{
  static double xcm_prev = 0., t_prev = 0.;
  double Vgas = 0., xmom = 0.;

  foreach(reduction(+:Vgas) reduction(+:xmom)) {
    double g = (1. - f[])*cs[];
    if (g > 0.) {
      double dV = g*dv();
      Vgas += dV;
      xmom += x*dV;
    }
  }
  double xcm = (Vgas > 0 ? xmom/Vgas : 0.);
  double Ub = 0.;
  if (i > 1) {
    double dx = xcm - xcm_prev;
    if (dx >  0.5*Lx) dx -= Lx;
    if (dx < -0.5*Lx) dx += Lx;
    Ub = dx/max(1e-30, (t - t_prev));
    xcm_unwrap += dx;
  } else
    xcm_unwrap = xcm;

  double h = film_thickness_global ();
  double h_over_R = h/Rtube;

  double h_taylor = taylor_law(Ca);
  double rel_error = fabs(h_over_R - h_taylor)/max(1e-30, h_taylor);

  if (i == 0 || i % 20 == 0)
    fprintf (stderr,
             "t=%8.4f Ca=%g h/R=% .8e Ub=% .6e xcm=% .6e Vgas=% .6e Taylor=% .6e\n",
             t, Ca, h_over_R, Ub, xcm, Vgas, h_taylor);

  static int wrote_header = 0;
  if (i == 0 || i % 10 == 0) {
    FILE * fp = fopen ("intermediate/summary.csv", (wrote_header ? "a" : "w"));
    if (fp) {
      if (!wrote_header) {
        fprintf(fp, "t,Ca,h_inf_over_R,Ub,xcm,Vgas,rel_err_to_Taylor\n");
        wrote_header = 1;
      }
      fprintf(fp, "%.8g,%.8g,%.12g,%.12g,%.12g,%.12g,%.12g\n",
              t, Ca, h_over_R, Ub, xcm, Vgas, rel_error);
      fclose(fp);
    }
  }

  if (Vgas_ref < 0.)
    Vgas_ref = Vgas;

  if (t >= t_steady_start) {
    sum_h += h_over_R;
    sum_Ub += Ub;
    sum_xcm += xcm;
    sum_Vgas += Vgas;
    n_steady++;
  }

  xcm_prev = xcm;
  t_prev = t;
  return 0;
}

/* ---------------- summary output ---------------- */

event summary (t = t_end)
{
  double h_mean = (n_steady > 0 ? sum_h/n_steady : 0.);
  double Ub_mean = (n_steady > 0 ? sum_Ub/n_steady : 0.);
  double Vgas_mean = (n_steady > 0 ? sum_Vgas/n_steady : 0.);
  double h_taylor = taylor_law(Ca);
  double rel_error = fabs(h_mean - h_taylor)/max(1e-30, h_taylor);
  double vgas_drift = (Vgas_ref > 0. ? (Vgas_mean - Vgas_ref)/Vgas_ref : 0.);

  FILE * fp = fopen ("intermediate/final.csv", "w");
  if (fp) {
    fprintf(fp, "Ca,h_inf_over_R,h_Taylor_over_R,rel_error,Ub_mean,Vgas_mean,Vgas_drift\n");
    fprintf(fp, "%.8g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g\n",
            Ca, h_mean, h_taylor, rel_error, Ub_mean, Vgas_mean, vgas_drift);
    fclose(fp);
  }

  FILE * ft = fopen ("intermediate/final.txt", "w");
  if (ft) {
    fprintf(ft, "Ca=%g h_mean/R=%g Taylor=%g rel_error=%g Ub_mean=%g Vgas_mean=%g Vgas_drift=%g\n",
            Ca, h_mean, h_taylor, rel_error, Ub_mean, Vgas_mean, vgas_drift);
    fclose(ft);
  }

  fprintf(stderr,
          "[final] Ca=%g h_mean/R=% .8e Taylor=% .8e rel_error=% .8e Ub_mean=% .6e Vgas_drift=% .3e\n",
          Ca, h_mean, h_taylor, rel_error, Ub_mean, vgas_drift);
  return 0;
}

/* ---------------- output (optional quick-look images) ---------------- */

event snapshots (t += 0.05; t <= t_end)
{
  // write dumps if you want restart/postprocess
  char name[256];
  sprintf(name, "intermediate/dump-%g", t);
  dump (file = name);
  return 0;
}

event movie (t += 0.02; t <= t_end)
{
  output_ppm (f, file = "intermediate/clean_case.mp4",
              min = 0., max = 1., n = 600, linear = true);
  return 0;
}

/* ---------------- no AMR ---------------- */
