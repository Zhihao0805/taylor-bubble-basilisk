/**
 * taylor_benchmark_2Dpaper.c
 *
 * Section-4-style 2D Taylor bubble benchmark driver:
 * - Inlet: constant velocity
 * - Outlet: fixed pressure
 * - Walls: no-slip
 * - Re fixed by case, Ca swept
 * - AMR on velocity + VOF fraction (paper also used cs; see report for deviation)
 *
 * CLI:
 *   ./run2d LEV Re Ca TMAX POISSON_TOL TAG [MAKE_MOVIE] [EPS_U] [EPS_F] [WRITE_DUMPS]
 */

#include "grid/quadtree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "output.h"

/* ---------------- run parameters ---------------- */
static int LEV = 10;
static int MINLEV = 6;
static double Re = 0.1;
static double Ca = 1e-3;
static double t_end = 0.6;
static double poisson_tol = 1e-5;
static double eps_u = 7e-3;
static double eps_f = 1e-1;
static int make_movie = 0;
static int write_dumps = 0;

/* Geometry/nondimensionalization */
static const double H = 1.0;      /* channel height */
static const double Lx = 10.0;    /* streamwise length */
static const double Uin = 0.1;    /* reference/inlet velocity */

/* Output tags/paths */
static char tag[128] = "case";
static char ts_path[256], final_path[256], movie_path[256];

/* Helpers for moving-centered movie */
scalar f_centered[];

/* Steady-window accumulation (last 30% of runtime) */
static double sum_uduf = 0., sum_ud = 0., sum_vdrop = 0.;
static int n_steady = 0;

/* ---------------- boundaries ---------------- */
u.n[left]  = dirichlet(Uin);
u.t[left]  = dirichlet(0.);
f[left]    = dirichlet(0.); /* continuous phase enters at inlet */

p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);
u.n[right] = neumann(0.);
u.t[right] = neumann(0.);
f[right]   = neumann(0.);

u.n[bottom] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);
u.n[top]    = dirichlet(0.);
u.t[top]    = dirichlet(0.);

/* ---------------- theory relations (paper Eq. 7-9) ---------------- */
static inline double eq7_two_t_over_H (double Ca_) {
  return 0.643*pow(3.*Ca_, 2./3.);
}

static inline double eq8_two_t_over_H (double Ca_) {
  double num = eq7_two_t_over_H(Ca_);
  return num/(1. + 0.643*2.5*pow(3.*Ca_, 2./3.));
}

static inline double eq9_ud_over_uf_from_two_t_over_H (double two_t_over_H) {
  double denom = 1. - two_t_over_H;
  if (denom < 1e-12) denom = 1e-12;
  return 1./denom;
}

/* ---------------- diagnostics helpers ---------------- */
static double bubble_stats (double * Ud, double * Vdrop) {
  double v = 0., umom = 0., xmom = 0.;

  foreach (reduction(+:v) reduction(+:umom) reduction(+:xmom)) {
    if (y <= H && f[] > 1e-12) {
      double dV = f[]*dv();
      v += dV;
      umom += u.x[]*dV;
      xmom += x*dV;
    }
  }

  *Vdrop = v;
  *Ud = (v > 0. ? umom/v : 0.);
  return (v > 0. ? xmom/v : 0.);
}

static int fluid_leaf_cells (void) {
  int n = 0;
  foreach (reduction(+:n))
    if (is_leaf(cell) && y <= H)
      n++;
  return n;
}

/* ---------------- main ---------------- */
int main (int argc, char ** argv)
{
  if (argc > 1) LEV = atoi(argv[1]);
  if (argc > 2) Re = atof(argv[2]);
  if (argc > 3) Ca = atof(argv[3]);
  if (argc > 4) t_end = atof(argv[4]);
  if (argc > 5) poisson_tol = atof(argv[5]);
  if (argc > 6) snprintf(tag, sizeof(tag), "%s", argv[6]);
  if (argc > 7) make_movie = atoi(argv[7]);
  if (argc > 8) eps_u = atof(argv[8]);
  if (argc > 9) eps_f = atof(argv[9]);
  if (argc > 10) write_dumps = atoi(argv[10]);

  MINLEV = max(4, LEV - 4);
  TOLERANCE = poisson_tol;

  snprintf(ts_path, sizeof(ts_path), "intermediate/%s_timeseries.csv", tag);
  snprintf(final_path, sizeof(final_path), "intermediate/%s_final.csv", tag);
  snprintf(movie_path, sizeof(movie_path), "intermediate/%s_centered.mp4", tag);

  size (Lx);
  origin (0., 0.);
  init_grid (1 << MINLEV);

  /* Keep channel height = H inside a square computational box. */
  mask (y > H ? top : none);

  rho1 = 1.; rho2 = 1.;             /* density ratio = 1 */
  mu1 = Uin/Re; mu2 = Uin/Re;       /* viscosity ratio = 1 */
  f.sigma = mu1*Uin/Ca;             /* Ca = mu*U/sigma */

  CFL = 0.4;
  DT = 5e-4;

  fprintf (stderr,
           "[args] TAG=%s LEV=%d MINLEV=%d Re=%g Ca=%g TMAX=%g TOL=%g eps_u=%g eps_f=%g movie=%d dumps=%d\n",
           tag, LEV, MINLEV, Re, Ca, t_end, poisson_tol, eps_u, eps_f, make_movie, write_dumps);

  run();
}

/* ---------------- init ---------------- */
event init (t = 0)
{
  double a = 0.75*H;     /* streamwise semi-axis */
  double b = 0.45*H;     /* cross-stream semi-axis */
  double x0 = 2.5*H;
  double y0 = 0.5*H;

  foreach() {
    if (y <= H) {
      double xi = (x - x0)/a;
      double yi = (y - y0)/b;
      f[] = (xi*xi + yi*yi <= 1.) ? 1. : 0.;
      u.x[] = Uin;
      u.y[] = 0.;
    } else {
      f[] = 0.;
      u.x[] = u.y[] = 0.;
    }
  }
  boundary ({f, u.x, u.y});
}

/* ---------------- AMR ---------------- */
event adapt_mesh (i++)
{
  adapt_wavelet ((scalar *){u.x, u.y, f},
                 (double[]){eps_u, eps_u, eps_f},
                 LEV, MINLEV);
}

/* ---------------- time diagnostics ---------------- */
event diagnostics (t += 0.01; t <= t_end + 1e-12)
{
  static int wrote_header = 0;
  double Ud = 0., Vdrop = 0.;
  double xcm = bubble_stats (&Ud, &Vdrop);
  double Ud_over_Uf = Ud/Uin;
  double two_tH_ud = 1. - 1./max(Ud_over_Uf, 1e-12);

  double two_tH_b = eq7_two_t_over_H(Ca);
  double two_tH_a = eq8_two_t_over_H(Ca);
  double uduf_b = eq9_ud_over_uf_from_two_t_over_H(two_tH_b);
  double uduf_a = eq9_ud_over_uf_from_two_t_over_H(two_tH_a);
  int ncell = fluid_leaf_cells();

  FILE * fp = fopen(ts_path, wrote_header ? "a" : "w");
  if (fp) {
    if (!wrote_header) {
      fprintf(fp, "t,Re,Ca,LEV,TOL,Uf,Ud,Ud_over_Uf,two_t_over_H_from_ud,two_t_over_H_bretherton,two_t_over_H_aussillous,Ud_over_Uf_bretherton,Ud_over_Uf_aussillous,Ncells,Vdrop,xcm\n");
      wrote_header = 1;
    }
    fprintf(fp, "%.8g,%.8g,%.8g,%d,%.8g,%.8g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%d,%.12g,%.12g\n",
            t, Re, Ca, LEV, poisson_tol, Uin, Ud, Ud_over_Uf, two_tH_ud, two_tH_b, two_tH_a, uduf_b, uduf_a, ncell, Vdrop, xcm);
    fclose(fp);
  }

  if (t >= 0.7*t_end) {
    sum_uduf += Ud_over_Uf;
    sum_ud += Ud;
    sum_vdrop += Vdrop;
    n_steady++;
  }

  if (i == 0 || i % 25 == 0)
    fprintf(stderr, "t=%7.4f Ca=%g Ud/Uf=%.8g Ud=%.8g Vdrop=%.8g N=%d xcm=%.8g\n",
            t, Ca, Ud_over_Uf, Ud, Vdrop, ncell, xcm);
}

/* ---------------- centered movie ---------------- */
event movie (t += 0.02; t <= t_end + 1e-12)
{
  if (!make_movie)
    return 0;

  double Ud = 0., Vdrop = 0.;
  double xcm = bubble_stats (&Ud, &Vdrop);
  double shift = xcm - 0.5*Lx;

  foreach() {
    if (y <= H) {
      double xp = x + shift;
      if (xp < 0.) xp = 0.;
      if (xp > Lx) xp = Lx;
      f_centered[] = interpolate (f, xp, y);
    } else
      f_centered[] = 0.;
  }
  boundary ({f_centered});

  output_ppm (f_centered, file = movie_path,
              min = 0., max = 1., n = 800, linear = true);
}

/* ---------------- optional dumps ---------------- */
event snapshots (t += 0.05; t <= t_end + 1e-12)
{
  if (!write_dumps)
    return 0;

  char name[256];
  sprintf(name, "intermediate/%s_dump-%g", tag, t);
  dump (file = name);
}

/* ---------------- final summary ---------------- */
event summary (t = end)
{
  double Ud = 0., Vdrop = 0.;
  double xcm = bubble_stats (&Ud, &Vdrop);
  double Ud_over_Uf = Ud/Uin;
  double two_tH_ud = 1. - 1./max(Ud_over_Uf, 1e-12);
  double two_tH_b = eq7_two_t_over_H(Ca);
  double two_tH_a = eq8_two_t_over_H(Ca);
  double uduf_b = eq9_ud_over_uf_from_two_t_over_H(two_tH_b);
  double uduf_a = eq9_ud_over_uf_from_two_t_over_H(two_tH_a);

  double UdUf_mean = (n_steady > 0 ? sum_uduf/n_steady : Ud_over_Uf);
  double Vdrop_mean = (n_steady > 0 ? sum_vdrop/n_steady : Vdrop);
  double rel_err_b = fabs(UdUf_mean - uduf_b)/max(1e-12, uduf_b);
  double rel_err_a = fabs(UdUf_mean - uduf_a)/max(1e-12, uduf_a);

  FILE * fp = fopen(final_path, "w");
  if (fp) {
    fprintf(fp, "tag,Re,Ca,LEV,TOL,Uf,Ud_last,Ud_over_Uf_last,Ud_over_Uf_mean_last30,two_t_over_H_from_ud_last,two_t_over_H_bretherton,two_t_over_H_aussillous,Ud_over_Uf_bretherton,Ud_over_Uf_aussillous,rel_err_bretherton,rel_err_aussillous,Vdrop_last,Vdrop_mean_last30,xcm_last,n_steady_samples\n");
    fprintf(fp, "%s,%.8g,%.8g,%d,%.8g,%.8g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%d\n",
            tag, Re, Ca, LEV, poisson_tol, Uin, Ud, Ud_over_Uf, UdUf_mean,
            two_tH_ud, two_tH_b, two_tH_a, uduf_b, uduf_a, rel_err_b, rel_err_a,
            Vdrop, Vdrop_mean, xcm, n_steady);
    fclose(fp);
  }

  fprintf(stderr,
          "[final] tag=%s Re=%g Ca=%g LEV=%d TOL=%g Ud/Uf_mean(last30%%)=%.8g err_b=%.4g err_a=%.4g\n",
          tag, Re, Ca, LEV, poisson_tol, UdUf_mean, rel_err_b, rel_err_a);
}
