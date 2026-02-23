#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_eqn_type.h>
#include <gkyl_fem_poisson_bctype.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_gyrokinetic_run.h>
#include <gkyl_math.h>

#include <rt_arg_parse.h>

// State of the pseudo orbit-averaged integrator.
enum gk_multirate_state {
  GK_MULTIRATE_NONE = 0, // Haven't started.
  GK_MULTIRATE_SEP, // Slow electrons phase.
  GK_MULTIRATE_FDP, // Full dynamics phase.
  GK_MULTIRATE_COMPLETED, // Finished simulation.
};

struct gk_multirate_phase_params {
  enum gk_multirate_state phase; // Type of phase.
  int num_frames; // Number of frames.
  double duration; // Duration.
  double alpha; // Factor multiplying electron df/dt.
  bool is_static_field; // Whether to evolve the field.
  enum gkyl_gyrokinetic_fdot_multiplier_type fdot_mult_type; // Type of df/dt multipler.
};

// Define the context of the simulation. This is basically all the globals
struct gk_ltx_ctx {
  int cdim, vdim; // Dimensionality.

  // Geometry and magnetic field.
  double Lz;        // Domain size along magnetic field.
  double z_min;  double z_max;
  double psi_min;  double psi_max;
  double psi_LCFS;    // Radial location of the last closed flux surface.
  double Lx_core; // Radial extent of core in psi
  double Lx; // Total radial extent in psi
  // Plasma parameters.
  double me;  double qe;
  double mi;  double qi;
  double n0;  double Te0;  double Ti0; 

  // Collisions.
  double nuFrac;
  double nuElc;  double nuIon;
  double nuElcIon;  double nuIonElc;

  double anom_diff_D; // Anomalous particle diffusivity.

  // Source parameters.
  double n_srcOMP;        // Amplitude of the OMP source
  double x_srcOMP;        // Radial location of the OMP source.
  double Te_srcOMP;       // Te for the OMP source.
  double Ti_srcOMP;       // Ti for the OMP source.
  double sigma_srcOMP;    // Radial spread of the OMP source.
  double floor_src;       // Source floor.

  // Grid parameters.
  int Nz;
  int Nvpar;
  int Nmu;
  int cells[GKYL_MAX_DIM]; // Number of cells in all directions.
  int poly_order;

  // Physical velocity space limits
  double vpar_max_elc; // Parallel velocity extents for electrons.
  double mu_max_elc; // Maximum magnetic moment for electrons.
  double vpar_max_ion; // Parallel velocity extents for ions.
  double mu_max_ion; // Maximum magnetic moment for ions.

  // Computational velocity space limits
  double vpar_min_elc_c, vpar_max_elc_c;
  double mu_min_elc_c, mu_max_elc_c;
  double vpar_min_ion_c, vpar_max_ion_c;
  double mu_min_ion_c, mu_max_ion_c;

  double alpha_sep; // Multiplier for the electron df/dt.

  double t_end; // End time.
  int num_frames; // Number of output frames.
  int num_phases; // Number of phases.
  struct gk_multirate_phase_params *mr_phases; // Phases to run.  
  double write_phase_freq; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  double int_diag_calc_freq; // Frequency of calculating integrated diagnostics (as a factor of num_frames).  
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

// Z is constant at -0.8
// R goes from 1.5 to 1.75
void pfunc_upper(double s, double* RZ){
//  // Vertical plate.
//  RZ[0] = 0.14047;
//  RZ[1] = -(s-0.061)*0.6;
  // Tilted plate.
  double p0[2] = {0.14047, 0.0};
  double p1[2] = {0.182, -0.327};
  RZ[0] = (1-s)*p0[0]+s*p1[0];
  RZ[1] = (1-s)*p0[1]+s*p1[1];
  // Parametrized curve.
  const int npts = 100;
  double t[] = {
    0.        , 0.01010101, 0.02020202, 0.03030303, 0.04040404, 0.05050505,
    0.06060606, 0.07070707, 0.08080808, 0.09090909, 0.1010101 , 0.11111111,
    0.12121212, 0.13131313, 0.14141414, 0.15151515, 0.16161616, 0.17171717,
    0.18181818, 0.19191919, 0.2020202 , 0.21212121, 0.22222222, 0.23232323,
    0.24242424, 0.25252525, 0.26262626, 0.27272727, 0.28282828, 0.29292929,
    0.3030303 , 0.31313131, 0.32323232, 0.33333333, 0.34343434, 0.35353535,
    0.36363636, 0.37373737, 0.38383838, 0.39393939, 0.4040404 , 0.41414141,
    0.42424242, 0.43434343, 0.44444444, 0.45454545, 0.46464646, 0.47474747,
    0.48484848, 0.49494949, 0.50505051, 0.51515152, 0.52525253, 0.53535354,
    0.54545455, 0.55555556, 0.56565657, 0.57575758, 0.58585859, 0.5959596 ,
    0.60606061, 0.61616162, 0.62626263, 0.63636364, 0.64646465, 0.65656566,
    0.66666667, 0.67676768, 0.68686869, 0.6969697 , 0.70707071, 0.71717172,
    0.72727273, 0.73737374, 0.74747475, 0.75757576, 0.76767677, 0.77777778,
    0.78787879, 0.7979798 , 0.80808081, 0.81818182, 0.82828283, 0.83838384,
    0.84848485, 0.85858586, 0.86868687, 0.87878788, 0.88888889, 0.8989899 ,
    0.90909091, 0.91919192, 0.92929293, 0.93939394, 0.94949495, 0.95959596,
    0.96969697, 0.97979798, 0.98989899, 1.        
  };
  double R[] = {
    0.14047   , 0.14174962, 0.14273474, 0.14346879, 0.14399523, 0.14435747,
    0.14459896, 0.14476313, 0.14489342, 0.14503325, 0.14522608, 0.14551532,
    0.14594398, 0.14652831, 0.14724303, 0.14805929, 0.14894823, 0.14988067,
    0.15082262, 0.15173639, 0.15258417, 0.15332818, 0.15393485, 0.15444675,
    0.15497222, 0.1556218 , 0.15649171, 0.15757559, 0.15882221, 0.16018015,
    0.16159799, 0.16303165, 0.16448183, 0.16596607, 0.16750196, 0.16910707,
    0.17079897, 0.17259525, 0.17451335, 0.17656052, 0.17872727, 0.18100255,
    0.1833753 , 0.1858363 , 0.18838289, 0.19101384, 0.19372794, 0.19652279,
    0.1993922 , 0.20232915, 0.20532667, 0.20839436, 0.21160025, 0.21502521,
    0.21875013, 0.22284552, 0.22728356, 0.23198191, 0.23685765, 0.24182934,
    0.2468618 , 0.25197457, 0.25719035, 0.2625318 , 0.26801411, 0.27363332,
    0.27938251, 0.28525469, 0.29123681, 0.29730416, 0.30343066, 0.3095913 ,
    0.31577417, 0.32197626, 0.32819474, 0.33442466, 0.34065371, 0.34686797,
    0.35305417, 0.35921492, 0.36536894, 0.37153565, 0.37772614, 0.38390769,
    0.39003327, 0.3960564 , 0.40196738, 0.40781674, 0.41366055, 0.41955084,
    0.42546712, 0.43132581, 0.43704124, 0.44254048, 0.44785751, 0.45307983,
    0.45829535, 0.46359198, 0.46905762, 0.47478018 
  };
  double Z[] = {
     0.        , -0.00614024, -0.01230155, -0.01848084, -0.02467501, -0.03088095,
    -0.03709557, -0.04331577, -0.04953845, -0.05576051, -0.06197885, -0.06819037,
    -0.07439199, -0.08058222, -0.08676199, -0.09293243, -0.09909469, -0.10525007,
    -0.11140214, -0.11755623, -0.12371771, -0.12989196, -0.13608376, -0.14228739,
    -0.14848797, -0.15467038, -0.16082139, -0.16694179, -0.17303848, -0.17911836,
    -0.18518835, -0.19125409, -0.19731337, -0.20336104, -0.20939197, -0.21540098,
    -0.22138293, -0.22733267, -0.23324507, -0.23911718, -0.24494961, -0.2507433 ,
    -0.2564992 , -0.26221759, -0.2678964 , -0.27353307, -0.27912502, -0.28467231,
    -0.29018372, -0.29566975, -0.30114095, -0.30659113, -0.31195538, -0.31715586,
    -0.32211473, -0.32676266, -0.33111144, -0.33521774, -0.33913873, -0.34293037,
    -0.34661107, -0.35015479, -0.35353293, -0.3567169 , -0.35967613, -0.36237505,
    -0.3647773 , -0.36684671, -0.36857033, -0.36998013, -0.37111323, -0.37200621,
    -0.37268893, -0.37318668, -0.37352468, -0.37371327, -0.37371068, -0.37346375,
    -0.37292306, -0.3721274 , -0.37120513, -0.37028866, -0.36947045, -0.3686337 ,
    -0.36759322, -0.36616572, -0.36430208, -0.3621725 , -0.35996742, -0.35786656,
    -0.35585573, -0.35375228, -0.35136798, -0.34854121, -0.34533453, -0.34192265,
    -0.33848114, -0.33518554, -0.33221142, -0.32973434
  };
  // Find indices in t that bound s.
  int idx_tlo, idx_tup;
  if (s < 1e-8) {
    idx_tlo = 0;
    idx_tup = 0;
  }
  else if (fabs(s-1.0) < 1e-8) {
    idx_tlo = npts-1;
    idx_tup = npts-1;
  }
  else {
    for (int i=0; i<npts-1; i++) {
      if (t[i] <= s && s < t[i+1]) {
        idx_tlo = i;
        idx_tup = i+1;
        break;
      }
    }
  }
  // Interpolate the value of R and Z.
  double Dt = t[idx_tup]-t[idx_tlo];
  if (idx_tlo == idx_tup) {
    RZ[0] = R[idx_tlo];
    RZ[1] = Z[idx_tlo];
  }
  else {
    RZ[0] = ((s-t[idx_tlo])/Dt)*R[idx_tup] + ((t[idx_tup]-s)/Dt)*R[idx_tlo];
    RZ[1] = ((s-t[idx_tlo])/Dt)*Z[idx_tup] + ((t[idx_tup]-s)/Dt)*Z[idx_tlo];
  }
}

void pfunc_lower(double s, double* RZ){
//  // Vertical plate.
//  RZ[0] = 0.14047;
//  RZ[1] = (s-0.061)*0.6;
//  // Tilted plate.
//  double p0[2] = {0.14047, 0.0};
//  double p1[2] = {0.182, 0.327};
//  RZ[0] = (1-s)*p0[0]+s*p1[0];
//  RZ[1] = (1-s)*p0[1]+s*p1[1];
  // Parametrized curve.
  const int npts = 100;
  double t[] = {
    0.        , 0.01010101, 0.02020202, 0.03030303, 0.04040404, 0.05050505,
    0.06060606, 0.07070707, 0.08080808, 0.09090909, 0.1010101 , 0.11111111,
    0.12121212, 0.13131313, 0.14141414, 0.15151515, 0.16161616, 0.17171717,
    0.18181818, 0.19191919, 0.2020202 , 0.21212121, 0.22222222, 0.23232323,
    0.24242424, 0.25252525, 0.26262626, 0.27272727, 0.28282828, 0.29292929,
    0.3030303 , 0.31313131, 0.32323232, 0.33333333, 0.34343434, 0.35353535,
    0.36363636, 0.37373737, 0.38383838, 0.39393939, 0.4040404 , 0.41414141,
    0.42424242, 0.43434343, 0.44444444, 0.45454545, 0.46464646, 0.47474747,
    0.48484848, 0.49494949, 0.50505051, 0.51515152, 0.52525253, 0.53535354,
    0.54545455, 0.55555556, 0.56565657, 0.57575758, 0.58585859, 0.5959596 ,
    0.60606061, 0.61616162, 0.62626263, 0.63636364, 0.64646465, 0.65656566,
    0.66666667, 0.67676768, 0.68686869, 0.6969697 , 0.70707071, 0.71717172,
    0.72727273, 0.73737374, 0.74747475, 0.75757576, 0.76767677, 0.77777778,
    0.78787879, 0.7979798 , 0.80808081, 0.81818182, 0.82828283, 0.83838384,
    0.84848485, 0.85858586, 0.86868687, 0.87878788, 0.88888889, 0.8989899 ,
    0.90909091, 0.91919192, 0.92929293, 0.93939394, 0.94949495, 0.95959596,
    0.96969697, 0.97979798, 0.98989899, 1.        
  };
  double R[] = {
    0.14047   , 0.14174962, 0.14273474, 0.14346879, 0.14399523, 0.14435747,
    0.14459896, 0.14476313, 0.14489342, 0.14503325, 0.14522608, 0.14551532,
    0.14594398, 0.14652831, 0.14724303, 0.14805929, 0.14894823, 0.14988067,
    0.15082262, 0.15173639, 0.15258417, 0.15332818, 0.15393485, 0.15444675,
    0.15497222, 0.1556218 , 0.15649171, 0.15757559, 0.15882221, 0.16018015,
    0.16159799, 0.16303165, 0.16448183, 0.16596607, 0.16750196, 0.16910707,
    0.17079897, 0.17259525, 0.17451335, 0.17656052, 0.17872727, 0.18100255,
    0.1833753 , 0.1858363 , 0.18838289, 0.19101384, 0.19372794, 0.19652279,
    0.1993922 , 0.20232915, 0.20532667, 0.20839436, 0.21160025, 0.21502521,
    0.21875013, 0.22284552, 0.22728356, 0.23198191, 0.23685765, 0.24182934,
    0.2468618 , 0.25197457, 0.25719035, 0.2625318 , 0.26801411, 0.27363332,
    0.27938251, 0.28525469, 0.29123681, 0.29730416, 0.30343066, 0.3095913 ,
    0.31577417, 0.32197626, 0.32819474, 0.33442466, 0.34065371, 0.34686797,
    0.35305417, 0.35921492, 0.36536894, 0.37153565, 0.37772614, 0.38390769,
    0.39003327, 0.3960564 , 0.40196738, 0.40781674, 0.41366055, 0.41955084,
    0.42546712, 0.43132581, 0.43704124, 0.44254048, 0.44785751, 0.45307983,
    0.45829535, 0.46359198, 0.46905762, 0.47478018 
  };
  double Z[] = {
     0.        , 0.00614024, 0.01230155, 0.01848084, 0.02467501, 0.03088095,
     0.03709557, 0.04331577, 0.04953845, 0.05576051, 0.06197885, 0.06819037,
     0.07439199, 0.08058222, 0.08676199, 0.09293243, 0.09909469, 0.10525007,
     0.11140214, 0.11755623, 0.12371771, 0.12989196, 0.13608376, 0.14228739,
     0.14848797, 0.15467038, 0.16082139, 0.16694179, 0.17303848, 0.17911836,
     0.18518835, 0.19125409, 0.19731337, 0.20336104, 0.20939197, 0.21540098,
     0.22138293, 0.22733267, 0.23324507, 0.23911718, 0.24494961, 0.2507433 ,
     0.2564992 , 0.26221759, 0.2678964 , 0.27353307, 0.27912502, 0.28467231,
     0.29018372, 0.29566975, 0.30114095, 0.30659113, 0.31195538, 0.31715586,
     0.32211473, 0.32676266, 0.33111144, 0.33521774, 0.33913873, 0.34293037,
     0.34661107, 0.35015479, 0.35353293, 0.3567169 , 0.35967613, 0.36237505,
     0.3647773 , 0.36684671, 0.36857033, 0.36998013, 0.37111323, 0.37200621,
     0.37268893, 0.37318668, 0.37352468, 0.37371327, 0.37371068, 0.37346375,
     0.37292306, 0.3721274 , 0.37120513, 0.37028866, 0.36947045, 0.3686337 ,
     0.36759322, 0.36616572, 0.36430208, 0.3621725 , 0.35996742, 0.35786656,
     0.35585573, 0.35375228, 0.35136798, 0.34854121, 0.34533453, 0.34192265,
     0.33848114, 0.33518554, 0.33221142, 0.32973434
  };
  // Find indices in t that bound s.
  int idx_tlo, idx_tup;
  if (s < 1e-8) {
    idx_tlo = 0;
    idx_tup = 0;
  }
  else if (fabs(s-1.0) < 1e-8) {
    idx_tlo = npts-1;
    idx_tup = npts-1;
  }
  else {
    for (int i=0; i<npts-1; i++) {
      if (t[i] <= s && s < t[i+1]) {
        idx_tlo = i;
        idx_tup = i+1;
        break;
      }
    }
  }
  // Interpolate the value of R and Z.
  double Dt = t[idx_tup]-t[idx_tlo];
  if (idx_tlo == idx_tup) {
    RZ[0] = R[idx_tlo];
    RZ[1] = Z[idx_tlo];
  }
  else {
    RZ[0] = ((s-t[idx_tlo])/Dt)*R[idx_tup] + ((t[idx_tup]-s)/Dt)*R[idx_tlo];
    RZ[1] = ((s-t[idx_tlo])/Dt)*Z[idx_tup] + ((t[idx_tup]-s)/Dt)*Z[idx_tlo];
  }
}

// Source profiles.
void density_srcOMP(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[1];

  struct gk_ltx_ctx *app = ctx;
  double n_srcOMP = app->n_srcOMP;
  double x_srcOMP = app->x_srcOMP;
  double sigma_srcOMP = app->sigma_srcOMP;
  double floor_src = app->floor_src;

  fout[0] = n_srcOMP * fmax(0.5*(1.0-tanh(-(x-x_srcOMP)/(2.0*sigma_srcOMP))), floor_src);
}

void zero_func(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void temp_elc_srcOMP(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[1];

  struct gk_ltx_ctx *app = ctx;
  double x_srcOMP = app->x_srcOMP;
  double sigma_srcOMP = app->sigma_srcOMP;
  double Te_srcOMP = app->Te_srcOMP;

  if (x_srcOMP - 5*sigma_srcOMP < x) {
    fout[0] = Te_srcOMP;
  } else {
    fout[0] = Te_srcOMP/10.0;
  }
}

void temp_ion_srcOMP(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[1];

  struct gk_ltx_ctx *app = ctx;
  double x_srcOMP = app->x_srcOMP;
  double sigma_srcOMP = app->sigma_srcOMP;
  double Ti_srcOMP = app->Ti_srcOMP;

  if (x_srcOMP - 5*sigma_srcOMP < x) {
    fout[0] = Ti_srcOMP;
  } else {
    fout[0] = Ti_srcOMP/10.0;
  }
}

// Ion initial conditions
void density_init(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[1];

  struct gk_ltx_ctx *app = ctx;
  double psi_max = app->psi_max;
  double Lx = app->Lx;
  double n0 = app->n0;

  fout[0] = 0.8*n0*(1.5 - tanh((9.0/Lx)*(psi_max-0.646767778*Lx-x)));
}

void upar_init(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[1];

  struct gk_ltx_ctx *app = ctx;

  fout[0] = 0.0;
}

void temp_init_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[1];

  struct gk_ltx_ctx *app = ctx;
  double psi_max = app->psi_max;
  double Lx = app->Lx;
  double Ti0 = app->Ti0;

  fout[0] = 0.8*Ti0*(1.1 - tanh((9.0/Lx)*(psi_max-0.646767778*Lx-x)));
}

void temp_init_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[1];

  struct gk_ltx_ctx *app = ctx;
  double psi_max = app->psi_max;
  double Lx = app->Lx;
  double Te0 = app->Te0;

  fout[0] = 0.8*Te0*(1.1 - tanh((9.0/Lx)*(psi_max-0.646767778*Lx-x)));
}

void density_bcx_lo(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  density_init(t, xn, fout, ctx);
}

void upar_bcx_lo(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  upar_init(t, xn, fout, ctx);
}

void temp_bcx_lo_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  temp_init_elc(t, xn, fout, ctx);
}

void temp_bcx_lo_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  temp_init_ion(t, xn, fout, ctx);
}

void
diffusion_D_func(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0], z = xn[1];

  struct gk_ltx_ctx *app = ctx;

  fout[0] = app->anom_diff_D; // Diffusivity [m^2/s].
}

// Velocity space mappings.
void mapc2p_vel_elc(double t, const double *vc, double* GKYL_RESTRICT vp, void *ctx)
{
  struct gk_ltx_ctx *app = ctx;
  double vpar_max_elc = app->vpar_max_elc;
  double mu_max_elc = app->mu_max_elc;

  double cvpar = vc[0], cmu = vc[1];
  // Linear map up to vpar_max/2, then quadratic.
  if (fabs(cvpar) <= 0.5)
    vp[0] = vpar_max_elc*cvpar;
  else if (cvpar < -0.5)
    vp[0] = -vpar_max_elc*2.0*pow(cvpar,2);
  else
    vp[0] =  vpar_max_elc*2.0*pow(cvpar,2);

  // Quadratic map in mu.
  vp[1] = mu_max_elc*pow(cmu,2);
}

void mapc2p_vel_ion(double t, const double *vc, double* GKYL_RESTRICT vp, void *ctx)
{
  struct gk_ltx_ctx *app = ctx;
  double vpar_max_ion = app->vpar_max_ion;
  double mu_max_ion = app->mu_max_ion;

  double cvpar = vc[0], cmu = vc[1];
  // Linear map up to vpar_max/2, then quadratic.
  if (fabs(cvpar) <= 0.5)
    vp[0] = vpar_max_ion*cvpar;
  else if (cvpar < -0.5)
    vp[0] = -vpar_max_ion*2.0*pow(cvpar,2);
  else
    vp[0] =  vpar_max_ion*2.0*pow(cvpar,2);

  // Quadratic map in mu.
  vp[1] = mu_max_ion*pow(cmu,2);
}

void
fdot_multiplier_profile_elc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0], z = xn[1];

  struct gk_ltx_ctx *app = ctx;

  fout[0] = app->alpha_sep;
}

struct gk_ltx_ctx
create_ctx(void)
{
  int cdim = 2, vdim = 2; // Dimensionality.

  // Universal constant parameters.
  double eps0 = GKYL_EPSILON0, eV = GKYL_ELEMENTARY_CHARGE;
  double mp = GKYL_PROTON_MASS, me = GKYL_ELECTRON_MASS;
  double qi = eV; // ion charge
  double qe = -eV; // electron charge

  // Grid parameters
  int Nx = 16;
  double Nx_core_frac = (5.0/8.0);
  int Nx_core = floor(Nx_core_frac*Nx);
  int Nx_sol = Nx-Nx_core;
  int Nz = 16;
  int Nvpar = 12;
  int Nmu = 8;
  int poly_order = 1;

  // Geometry and magnetic field.
  double Lz = 2.0*(M_PI-1e-14);    // Domain size along magnetic field.
  double psi_LCFS = -5.4760172700000003e-03; // psi at LCFS. Taken from efit.
  double psi_max = psi_LCFS + 0.003; // Inner core boundary (psi increases towards core).
  double delta_psi = (psi_max-psi_LCFS)/Nx_core; // Cell length in psi.
  double psi_min = psi_LCFS - Nx_sol*delta_psi; // Outer SOL boundary (psi increases towards core).
  double Lx = psi_max - psi_min;
  double Lx_core = psi_LCFS - psi_min;
  printf("Nx_core          = %3d\n",Nx_core);
  printf("Nx_sol           = %3d\n",Nx_sol);
  printf("Nx               = %3d\n",Nx);
  printf("psi_LCFS         = %.9e\n",psi_LCFS        );
  printf("psi_min          = %.9e\n",psi_min         );
  printf("psi_max          = %.9e\n",psi_max         );
  printf("psi_LCFS-psi_min = %.9e\n",psi_LCFS-psi_min);
  printf("psi_max-psi_LCFS = %.9e\n",psi_max-psi_LCFS);
  printf("\n");

  // Plasma parameters. Chosen based on the value of a cubic sline
  // between the last TS data inside the LCFS and the probe data in
  // in the far SOL, near R=0.475 m.
  double B0  = 0.5*(1.666901e-01+6.814593e-01);
  double mi  = mp;   // Hydrogen ions.
  double Te_max = 320*eV;
  double Ti_max = Te_max/2.0;
  double n_max  = 1.6e19;   // [1/m^3]

  double Te0 = 0.5*Te_max;
  double Ti0 = 0.5*Ti_max;
  double n0  = 0.5*n_max;   // [1/m^3]

  double nuFrac = 1.0;
  // Electron-electron collision freq.
  double logLambdaElc = 6.6 - 0.5 * log(n0/1e20) + 1.5 * log(Ti0/eV);
  double nuElc = nuFrac * logLambdaElc * pow(eV, 4) * n0 /
    (6*sqrt(2.) * pow(M_PI,3./2.) * pow(eps0,2) * sqrt(me) * pow(Te0,3./2.));
  // Ion-ion collision freq.
  double logLambdaIon = 6.6 - 0.5 * log(n0/1e20) + 1.5 * log(Ti0/eV);
  double nuIon = nuFrac * logLambdaIon * pow(eV, 4) * n0 /
    (12 * pow(M_PI,3./2.) * pow(eps0,2) * sqrt(mi) * pow(Ti0,3./2.));
  double nuElcIon = nuElc*sqrt(2.0);
  double nuIonElc = nuElcIon*(me/mi);

  double anom_diff_D = 30.0; // Anomalous particle diffusivity [m^2/s].

  double vte = sqrt(Te0/me), vti = sqrt(Ti0/mi); // Thermal speeds.

  double c_s = sqrt(Te0/mi);
  double omega_ci = fabs(qi*B0/mi);
  double rho_s = c_s/omega_ci;

  // Source parameters
  // TRANSP estimated P_OH=179.19e5 kW for this shot.
  double n_LCFS = 0.3e19;
  double Lc_LCFS = 6.0; // Approx. connection length just outside LCFS at OMP.
  double Te_LCFS = 210.0*eV;
  double cs_LCFS = sqrt(Te_LCFS/mi);
  double n_srcOMP = n_LCFS*cs_LCFS/Lc_LCFS;
  double x_srcOMP = psi_max-0.15*Lx;
  double Te_srcOMP = Te_max;
  double Ti_srcOMP = Ti_max;
  double sigma_srcOMP = 0.03*Lx;
  double floor_src = 1e-2;

  // Physical velocity space limits
  double vpar_max_elc = 5.*vte;
  double mu_max_elc = me*pow(4*vte,2)/(2*B0);
  double vpar_max_ion = 5.*vti;
  double mu_max_ion = mi*pow(4*vti,2)/(2*B0);

  // Computational velocity space limits.
  double vpar_min_ion_c = -1.0/sqrt(2.0);
  double vpar_max_ion_c =  1.0/sqrt(2.0);
  double mu_min_ion_c = 0.;
  double mu_max_ion_c = 1.;
  double vpar_min_elc_c = -1.0/sqrt(2.0);
  double vpar_max_elc_c =  1.0/sqrt(2.0);
  double mu_min_elc_c = 0.;
  double mu_max_elc_c = 1.;

  // Factor multiplying collisionless terms.
  double alpha_sep = 0.1;
  double alpha_fdp = 1.0;
  // Duration of each phase.
  double tau_sep = 500.0e-6;
  double tau_fdp = 50.0e-6;
  double tau_fdp_extra = 2*tau_fdp;
  int num_cycles = 2; // Number of SEP+FDP cycles to run.

  // Frame counts for each phase type (specified independently)
  int num_frames_sep = 10; // Frames per SEP phase.
  int num_frames_fdp = 1; // Frames per FDP phase.
  int num_frames_fdp_extra = 2*num_frames_fdp;  // Frames for the extra FDP phase

  // Whether to evolve the field.
  bool is_static_field_sep = false;
  bool is_static_field_fdp = false;
  // Type of df/dt multipler.
  enum gkyl_gyrokinetic_fdot_multiplier_type fdot_mult_type_sep = GKYL_GK_FDOT_MULTIPLIER_USER_INPUT;
  enum gkyl_gyrokinetic_fdot_multiplier_type fdot_mult_type_fdp = GKYL_GK_FDOT_MULTIPLIER_NONE;

  // Calculate phase structure
  double t_end = (tau_sep + tau_fdp)*num_cycles + tau_fdp_extra;
  double tau_pair = tau_sep+tau_fdp; // Duration of an sep+FDP pair.
  int num_phases = 2*num_cycles + 1;
  int num_frames = num_cycles * (num_frames_sep + num_frames_fdp) + num_frames_fdp_extra;

  struct gk_multirate_phase_params *mr_phases = gkyl_malloc(num_phases * sizeof(struct gk_multirate_phase_params));
  for (int i=0; i<(num_phases-1)/2; i++) {
    // SEPs.
    mr_phases[2*i].phase = GK_MULTIRATE_SEP;
    mr_phases[2*i].num_frames = num_frames_sep;
    mr_phases[2*i].duration = tau_sep;
    mr_phases[2*i].alpha = alpha_sep;
    mr_phases[2*i].is_static_field = is_static_field_sep;
    mr_phases[2*i].fdot_mult_type = fdot_mult_type_sep;

    // FDPs.
    mr_phases[2*i+1].phase = GK_MULTIRATE_FDP;
    mr_phases[2*i+1].num_frames = num_frames_fdp;
    mr_phases[2*i+1].duration = tau_fdp;
    mr_phases[2*i+1].alpha = alpha_fdp;
    mr_phases[2*i+1].is_static_field = is_static_field_fdp;
    mr_phases[2*i+1].fdot_mult_type = fdot_mult_type_fdp;
  }
  // Add an extra, longer FDP.
  mr_phases[num_phases-1].phase = GK_MULTIRATE_FDP;
  mr_phases[num_phases-1].num_frames = num_frames_fdp_extra;
  mr_phases[num_phases-1].duration = tau_fdp_extra;
  mr_phases[num_phases-1].alpha = alpha_fdp;
  mr_phases[num_phases-1].is_static_field = is_static_field_fdp;
  mr_phases[num_phases-1].fdot_mult_type = fdot_mult_type_fdp;

  double write_phase_freq = 1.0; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  double int_diag_calc_freq = 5; // Frequency of calculating integrated diagnostics (as a factor of num_frames). 
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct gk_ltx_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,
    .Lz    = Lz    ,
    .z_min = -Lz/2.,  .z_max = Lz/2.,
    .psi_min = psi_min,  .psi_max = psi_max,
    .psi_LCFS = psi_LCFS,
    .Lx = Lx,
    .Lx_core = Lx_core,
  
    .me = me,  .qe = qe,
    .mi = mi,  .qi = qi,
    .n0 = n0,  .Te0 = Te0,  .Ti0 = Ti0,
  
    .nuFrac = nuFrac,
    .nuElc = nuElc,  .nuIon = nuIon,
    .nuElcIon = nuElcIon,  .nuIonElc = nuIonElc,
  
    .anom_diff_D = anom_diff_D,

    .n_srcOMP     = n_srcOMP    ,
    .Te_srcOMP    = Te_srcOMP   ,
    .Ti_srcOMP    = Ti_srcOMP   ,
    .x_srcOMP     = x_srcOMP    ,
    .sigma_srcOMP = sigma_srcOMP,
    .floor_src    = floor_src   ,
  
    .Nz = Nz,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .cells = {Nx, Nz, Nvpar, Nmu},
    .poly_order = poly_order,

    // Physical velocity space limits
    .vpar_max_elc = vpar_max_elc,
    .mu_max_elc = mu_max_elc,
    .vpar_max_ion = vpar_max_ion,
    .mu_max_ion = mu_max_ion,
    // Computational velocity space limits
    .vpar_min_elc_c = vpar_min_elc_c,
    .vpar_max_elc_c = vpar_max_elc_c,
    .mu_min_elc_c = mu_min_elc_c,
    .mu_max_elc_c = mu_max_elc_c,
    .vpar_min_ion_c = vpar_min_ion_c,
    .vpar_max_ion_c = vpar_max_ion_c,
    .mu_min_ion_c = mu_min_ion_c,
    .mu_max_ion_c = mu_max_ion_c,

    .alpha_sep = alpha_sep,

    .t_end = t_end,  .num_frames = num_frames,
    .num_phases = num_phases,
    .mr_phases = mr_phases,
    .write_phase_freq = write_phase_freq,
    .int_diag_calc_freq   = int_diag_calc_freq,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };
  return ctx;
}

void
release_ctx(struct gk_ltx_ctx *ctx)
{
  gkyl_free(ctx->mr_phases);
}

void
calc_integrated_diagnostics(struct gkyl_tm_trigger* iot, gkyl_gyrokinetic_app* app,
  double t_curr, bool force_calc, double dt)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_calc) {
    gkyl_gyrokinetic_app_calc_field_energy(app, t_curr);
    gkyl_gyrokinetic_app_calc_integrated_mom(app, t_curr);

    if ( !(dt < 0.0) )
      gkyl_gyrokinetic_app_save_dt(app, t_curr, dt);
  }
}

void
write_data(struct gkyl_tm_trigger* iot_conf, struct gkyl_tm_trigger* iot_phase,
  gkyl_gyrokinetic_app* app, double t_curr, bool force_write)
{
  bool trig_now_conf = gkyl_tm_trigger_check_and_bump(iot_conf, t_curr);
  if (trig_now_conf || force_write) {
    int frame = (!trig_now_conf) && force_write? iot_conf->curr : iot_conf->curr-1;
    gkyl_gyrokinetic_app_write_conf(app, t_curr, frame);

    gkyl_gyrokinetic_app_write_field_energy(app);
    gkyl_gyrokinetic_app_write_integrated_mom(app);
    gkyl_gyrokinetic_app_write_dt(app);
  }

  bool trig_now_phase = gkyl_tm_trigger_check_and_bump(iot_phase, t_curr);
  if (trig_now_phase || force_write) {
    int frame = (!trig_now_conf) && force_write? iot_conf->curr : iot_conf->curr-1;

    gkyl_gyrokinetic_app_write_phase(app, t_curr, frame);
  }
}

struct time_frame_state {
  double t_curr; // Current simulation time.
  double t_end; // End time of current phase.
  int frame_curr; // Current frame.
  int num_frames; // Number of frames at the end of current phase.
};

void reset_io_triggers(struct gk_ltx_ctx *ctx, struct time_frame_state *tfs,
  struct gkyl_tm_trigger *trig_write_conf, struct gkyl_tm_trigger *trig_write_phase,
  struct gkyl_tm_trigger *trig_calc_intdiag)
{
  // Reset I/O triggers:
  double t_curr = tfs->t_curr;
  double t_end = tfs->t_end;
  int frame_curr = tfs->frame_curr;
  int num_frames = tfs->num_frames;
  int num_int_diag_calc = ctx->int_diag_calc_freq*num_frames;

  // Prevent division by zero when frame_curr equals num_frames
  int frames_remaining = num_frames - frame_curr;
  double time_remaining = t_end - t_curr;

  trig_write_conf->dt = time_remaining / frames_remaining;
  trig_write_conf->tcurr = t_curr;
  trig_write_conf->curr = frame_curr;

  trig_write_phase->dt = time_remaining / (ctx->write_phase_freq * frames_remaining);
  trig_write_phase->tcurr = t_curr;
  trig_write_phase->curr = frame_curr;

  int diag_frames = GKYL_MAX2(frames_remaining, (num_int_diag_calc/num_frames) * frames_remaining);
  trig_calc_intdiag->dt = time_remaining / diag_frames;
  trig_calc_intdiag->tcurr = t_curr;
  trig_calc_intdiag->curr = frame_curr;
}

void run_phase(gkyl_gyrokinetic_app* app, struct gk_ltx_ctx *ctx, double num_steps,
  struct gkyl_tm_trigger *trig_write_conf, struct gkyl_tm_trigger *trig_write_phase,
  struct gkyl_tm_trigger *trig_calc_intdiag,  struct time_frame_state *tfs,
  struct gk_multirate_phase_params *pparams)
{
  tfs->t_end = tfs->t_curr + pparams->duration;
  tfs->num_frames = tfs->frame_curr + pparams->num_frames;

  // Run an SEP or FDP.
  double t_curr = tfs->t_curr;
  double t_end = tfs->t_end;

  // Reset I/O triggers:
  reset_io_triggers(ctx, tfs, trig_write_conf, trig_write_phase, trig_calc_intdiag);

  // Reset simulation parameters and function pointers.
  struct gkyl_gyrokinetic_fdot_multiplier fdot_mult_inp = {
    .type = pparams->fdot_mult_type,
    .profile = fdot_multiplier_profile_elc,
    .profile_ctx = &ctx,
    .cellwise_const = true,
    .write_diagnostics = true,
  };
  gkyl_gyrokinetic_app_reset_species_fdot_multiplier(app, t_curr, "elc", fdot_mult_inp);

  // Compute initial guess of maximum stable time-step.
  double dt = t_end - t_curr;

  // Initialize small time-step check.
  double dt_init = -1.0, dt_failure_tol = ctx->dt_failure_tol;
  int num_failures = 0, num_failures_max = ctx->num_failures_max;

  long step = 1;
  while ((t_curr < t_end) && (step <= num_steps))
  {
    if (step == 1 || step % 10 == 0)
      gkyl_gyrokinetic_app_cout(app, stdout, "Taking time-step at t = %g ...", t_curr);

    dt = fmin(dt, t_end - t_curr); // Don't step beyond t_end.
    struct gkyl_update_status status = gkyl_gyrokinetic_update(app, dt);

    if (step == 1 || step % 10 == 0)
      gkyl_gyrokinetic_app_cout(app, stdout, " dt = %g\n", status.dt_actual);

    if (!status.success)
    {
      gkyl_gyrokinetic_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }
    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    calc_integrated_diagnostics(trig_calc_intdiag, app, t_curr, t_curr > t_end, status.dt_actual);
    write_data(trig_write_conf, trig_write_phase, app, t_curr, t_curr > t_end);

    if (dt_init < 0.0) {
      dt_init = status.dt_actual;
    }
    else if (status.dt_actual < dt_failure_tol * dt_init) {
      num_failures += 1;

      gkyl_gyrokinetic_app_cout(app, stdout, "WARNING: Time-step dt = %g", status.dt_actual);
      gkyl_gyrokinetic_app_cout(app, stdout, " is below %g*dt_init ...", dt_failure_tol);
      gkyl_gyrokinetic_app_cout(app, stdout, " num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        gkyl_gyrokinetic_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
        gkyl_gyrokinetic_app_cout(app, stdout, "%d consecutive times. Aborting simulation ....\n", num_failures_max);
        calc_integrated_diagnostics(trig_calc_intdiag, app, t_curr, true, status.dt_actual);
        write_data(trig_write_conf, trig_write_phase, app, t_curr, true);
        break;
      }
    }
    else {
      num_failures = 0;
    }

    step += 1;
  }

  tfs->t_curr = t_curr;
  tfs->frame_curr = tfs->frame_curr+pparams->num_frames;
}

int main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) MPI_Init(&argc, &argv);
#endif

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  struct gk_ltx_ctx ctx = create_ctx(); // Context for init functions.

  int cells_x[ctx.cdim], cells_v[ctx.vdim];
  for (int d=0; d<ctx.cdim; d++)
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  for (int d=0; d<ctx.vdim; d++)
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim+d]);

  // Construct communicator for use in app.
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);

  // Electrons.
  struct gkyl_gyrokinetic_projection elc_bcx_lo = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
    .ctx_density = &ctx,
    .density = density_bcx_lo,
    .ctx_upar = &ctx,
    .upar = upar_bcx_lo,
    .ctx_temp = &ctx,
    .temp = temp_bcx_lo_elc,
  };

  struct gkyl_gyrokinetic_species elc = {
    .name = "elc",
    .charge = ctx.qe, .mass = ctx.me,
    .vdim = ctx.vdim,
    .lower = { ctx.vpar_min_elc_c, ctx.mu_min_elc_c},
    .upper = { ctx.vpar_max_elc_c, ctx.mu_max_elc_c},
    .cells = { cells_v[0], cells_v[1] },

    .mapc2p = {
      .mapping = mapc2p_vel_elc,
      .ctx = &ctx,
    },

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = density_init,
      .ctx_upar = &ctx,
      .upar= upar_init,
      .ctx_temp = &ctx,
      .temp = temp_init_elc,      
    },

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .num_cross_collisions = 1,
      .collide_with = { "ion" },
      .den_ref = ctx.n0,
      .temp_ref = ctx.Te0,
    },

    .anomalous_diffusion = {
      .anomalous_diff_id = GKYL_GK_ANOMALOUS_DIFF_D,
      .D_profile = diffusion_D_func,
      .D_profile_ctx = &ctx,
//      .write_diagnostics = true,
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
        .ctx_density = &ctx,
        .ctx_upar = &ctx,
        .ctx_temp = &ctx,
        .density = density_srcOMP,
        .upar = zero_func,
        .temp = temp_elc_srcOMP,
      }, 
      .diagnostics = {
        .num_diag_moments = 5,
        .diag_moments = { GKYL_F_MOMENT_MAXWELLIAN, },
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = { GKYL_F_MOMENT_M0M1M2, },
      },
    },

    .time_rate_multiplier = {
      .type = GKYL_GK_FDOT_MULTIPLIER_USER_INPUT, // So solvers are allocated.
      .profile = fdot_multiplier_profile_elc,
      .profile_ctx = &ctx,
      .cellwise_const = true,
      .write_diagnostics = true,
    },

    .bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB, },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_FIXED_FUNC, .projection = elc_bcx_lo, },
      { .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_IWL, },
      { .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_IWL, },
    },

    .num_diag_moments = 4,
    .diag_moments = { GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_BIMAXWELLIAN, },
  };

  // Ions.
  struct gkyl_gyrokinetic_projection ion_bcx_lo = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
    .ctx_density = &ctx,
    .density = density_bcx_lo,
    .ctx_upar = &ctx,
    .upar = upar_bcx_lo,
    .ctx_temp = &ctx,
    .temp = temp_bcx_lo_ion,
  };

  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.qi, .mass = ctx.mi,
    .vdim = ctx.vdim,
    .lower = { ctx.vpar_min_ion_c, ctx.mu_min_ion_c},
    .upper = { ctx.vpar_max_ion_c, ctx.mu_max_ion_c},
    .cells = { cells_v[0], cells_v[1] },

    .mapc2p = {
      .mapping = mapc2p_vel_ion,
      .ctx = &ctx,
    },

    .polarization_density = ctx.n0,

//    .positivity = {
//      .type = GKYL_GK_POSITIVITY_SHIFT,
//      .quasineutrality_rescale = true,
//    },

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = density_init,
      .ctx_upar = &ctx,
      .upar= upar_init,
      .ctx_temp = &ctx,
      .temp = temp_init_ion,      
    },

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .num_cross_collisions = 1,
      .collide_with = { "elc" },
      .den_ref = ctx.n0,
      .temp_ref = ctx.Ti0,
    },

    .anomalous_diffusion = {
      .anomalous_diff_id = GKYL_GK_ANOMALOUS_DIFF_D,
      .D_profile = diffusion_D_func,
      .D_profile_ctx = &ctx,
//      .write_diagnostics = true,
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
        .ctx_density = &ctx,
        .ctx_upar = &ctx,
        .ctx_temp = &ctx,
        .density = density_srcOMP,
        .upar = zero_func,
        .temp = temp_ion_srcOMP,
      }, 
      .diagnostics = {
        .num_diag_moments = 5,
        .diag_moments = { GKYL_F_MOMENT_MAXWELLIAN, },
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = { GKYL_F_MOMENT_M0M1M2, },
      },
    },

    .bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB, },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_FIXED_FUNC, .projection = ion_bcx_lo, },
      { .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_IWL, },
      { .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_IWL, },
    },
    
    .num_diag_moments = 4,
    .diag_moments = { GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_BIMAXWELLIAN, },
  };

//  struct gkyl_poisson_bias_line target_corner_bcs[] = {
//    {
//     .perp_dirs = {0, 1}, // Directions perpendicular to line.
//     .perp_coords = {ctx.psi_LCFS, ctx.z_min}, // Coordinates of the line in perpendicular directions.
//     .val = 0.0, // Biasing value.
//    },
//    {
//     .perp_dirs = {0, 1}, // Directions perpendicular to line.
//     .perp_coords = {ctx.psi_LCFS, ctx.z_max}, // Coordinates of the line in perpendicular directions.
//     .val = 0.0, // Biasing value.
//    },
//  };
//
//  struct gkyl_poisson_bias_line_list bias_line_list = {
//    .num_bias_line = 2,
//    .bl = target_corner_bcs,
//  };

  // Field.
  struct gkyl_gyrokinetic_field field = {
    .gkfield_id = GKYL_GK_FIELD_ES_IWL,
    .poisson_bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0}, },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0}, },
    },
//    .bias_line_list = &bias_line_list,
  };

  struct gkyl_efit_inp efit_inp = {
    // psiRZ and related inputs
    .filepath = "../../../experiment/li_863mg_103955_03/LTX_103955_03.eqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
    .reflect = true,
  };

  struct gkyl_tok_geo_grid_inp grid_inp = {
    .ftype = GKYL_GEOMETRY_TOKAMAK_IWL,
    .rclose = 0.7,
    .rleft  = 0.1,
    .rright = 0.7,
    .rmin   = 0.1,
    .rmax   = 0.7,
    .zmin   = -0.3675,
    .zmax   =  0.3675,
    .plate_spec = true,
    .plate_func_lower = pfunc_lower,
    .plate_func_upper = pfunc_upper,
  }; 

  // GK app
  struct gkyl_gk app_inp = {
    .name = "gk_ltx_iwl_2x2v_p1",

    .cdim = ctx.cdim,
    .lower = { ctx.psi_min, ctx.z_min },
    .upper = { ctx.psi_max, ctx.z_max },
    .cells = { cells_x[0], cells_x[1] },
    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,

    .geometry = {
      .world = {0.0},
      .geometry_id = GKYL_GEOMETRY_TOKAMAK,
      .efit_info = efit_inp,
      .tok_grid_info = grid_inp,
      .has_LCFS = true,
      .x_LCFS = ctx.psi_LCFS, // Location of last closed flux surface.
    },

    .num_periodic_dir = 0,
    .periodic_dirs = {  },

    .num_species = 2,
    .species = { elc, ion },
    .field = field,

    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0], app_args.cuts[1] },
      .comm = comm,
    },
  };

  // Create app object.
  gkyl_gyrokinetic_app *app = gkyl_gyrokinetic_app_new(&app_inp);

  // Triggers for IO.
  struct gkyl_tm_trigger trig_write_conf, trig_write_phase, trig_calc_intdiag;

  struct time_frame_state tfs = {
    .t_curr = 0.0, // Initial simulation time.
    .frame_curr = 0, // Initial frame.
    .t_end = ctx.mr_phases[0].duration, // Final time of 1st phase.
    .num_frames = ctx.mr_phases[0].num_frames, // Number of frames in 1st phase.
  };

  int phase_idx_init = 0, phase_idx_end = ctx.num_phases; // Initial and final phase index.
  if (app_args.is_restart) {
    struct gkyl_app_restart_status status = gkyl_gyrokinetic_app_read_from_frame(app, app_args.restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_gyrokinetic_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n", gkyl_array_rio_status_msg(status.io_status));
      goto freeresources;
    }

    tfs.frame_curr = status.frame;
    tfs.t_curr = status.stime;

    // Find out what phase we are in.
    double time_count = 0.0;
    int frame_count = 0;
    int pit_curr = 0;
    for (int pit=0; pit<ctx.num_phases; pit++) {
      time_count += ctx.mr_phases[pit].duration;
      frame_count += ctx.mr_phases[pit].num_frames;
      if ((tfs.t_curr <= time_count) && (tfs.frame_curr <= frame_count)) {
        pit_curr = pit;
        break;
      }
    };
    phase_idx_init = pit_curr;

    // Change the duration and number frames so this phase reaches the expected
    // time and number of frames and not beyond.
    struct gk_multirate_phase_params *pparams = &ctx.mr_phases[phase_idx_init];
    pparams->num_frames = frame_count - tfs.frame_curr;
    pparams->duration = time_count - tfs.t_curr;

    gkyl_gyrokinetic_app_cout(app, stdout, "Restarting from frame %d", tfs.frame_curr);
    gkyl_gyrokinetic_app_cout(app, stdout, " at time = %g\n", tfs.t_curr);
  }
  else {
    gkyl_gyrokinetic_app_apply_ic(app, tfs.t_curr);

    // Write out ICs.
    reset_io_triggers(&ctx, &tfs, &trig_write_conf, &trig_write_phase, &trig_calc_intdiag);

    calc_integrated_diagnostics(&trig_calc_intdiag, app, tfs.t_curr, true, -1.0);
    write_data(&trig_write_conf, &trig_write_phase, app, tfs.t_curr, true);
  }

  if (app_args.num_steps != INT_MAX)
    phase_idx_end = 1;

  // Loop over number of number of phases;
  for (int pit=phase_idx_init; pit<phase_idx_end; pit++) {
    gkyl_gyrokinetic_app_cout(app, stdout, "\nRunning phase %d @ t = %.9e ... \n", pit, tfs.t_curr);
    struct gk_multirate_phase_params *phase_params = &ctx.mr_phases[pit];
    run_phase(app, &ctx, app_args.num_steps, &trig_write_conf, &trig_write_phase, &trig_calc_intdiag, &tfs, phase_params);
  }

  gkyl_gyrokinetic_app_stat_write(app);

  struct gkyl_gyrokinetic_stat stat = gkyl_gyrokinetic_app_stat(app); // fetch simulation statistics
  gkyl_gyrokinetic_app_cout(app, stdout, "\n");
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0)
  {
    gkyl_gyrokinetic_app_cout(app, stdout, "Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    gkyl_gyrokinetic_app_cout(app, stdout, "Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of write calls %ld\n", stat.n_io);
  gkyl_gyrokinetic_app_print_timings(app, stdout);

  freeresources:
  // simulation complete, free app
  gkyl_gyrokinetic_app_release(app);

  gkyl_gyrokinetic_comms_release(comm);
  release_ctx(&ctx);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Finalize();
#endif

  return 0;
}
