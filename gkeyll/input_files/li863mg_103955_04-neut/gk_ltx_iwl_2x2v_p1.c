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

// Define the context of the simulation. This is basically all the globals
struct gk_ltx_ctx {
  int cdim, vdim; // Dimensionality.

  char geqdsk_file[128]; // File with equilibrium.

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

  double rec_frac; // Recycling coefficient.
  double neut_gamma; // Neutral adiabatic index.
  double m_neut; // Neutral mass.
  double n0_neut, T0_neut; // Reference neutral density/temp.

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

  double t_end; // End time.
  int num_frames; // Number of output frames.
  double write_phase_freq; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

#define MAX_LTX_VV_INFO_POINTS 1000

struct ltx_lower_vessel_info {
  int num_points; // Number of points on curve.
  double t[MAX_LTX_VV_INFO_POINTS]; // index along the curve between 0 and 1.
  double R[MAX_LTX_VV_INFO_POINTS]; // Radius.
  double Z[MAX_LTX_VV_INFO_POINTS]; // Height.
};

void
get_ltx_lower_vessel_info(struct ltx_lower_vessel_info *vvi) {
  // Must go through (R(s=0),Z(s=0)) = (0.1411208988576905, 0.0).
  // May need manual modification to avoid anomalies too.
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
    0.1411209 , 0.14111931, 0.14112713, 0.14108865, 0.14367537, 0.14394567,
    0.14409344, 0.1443732 , 0.1449428 , 0.14584626, 0.14697714, 0.14817633,
    0.14932006, 0.15036927, 0.15140858, 0.15255488, 0.15389194, 0.15546736,
    0.15729399, 0.15937805, 0.16169684, 0.16415656, 0.16662234, 0.16902347,
    0.17141186, 0.17393792, 0.1767413 , 0.17989724, 0.18339635, 0.18722141,
    0.19132   , 0.19557265, 0.19987079, 0.20429364, 0.20911264, 0.21453395,
    0.22051232, 0.22689036, 0.23361348, 0.24076808, 0.24840878, 0.2564797 ,
    0.26487705, 0.27354237, 0.28245669, 0.29158629, 0.30086213, 0.31021888,
    0.31962488, 0.32906763, 0.33853077, 0.34798661, 0.35741528, 0.36681048,
    0.37616604, 0.38545935, 0.39465595, 0.40372805, 0.41267681, 0.42151796,
    0.43025304, 0.43884826, 0.44725614, 0.45544622, 0.46342614, 0.47121937,
    0.47884908, 0.48633188, 0.49368488, 0.50091932, 0.5080326 , 0.5150037 ,
    0.52181019, 0.52844457, 0.53491582, 0.5412411 , 0.54743706, 0.55352704,
    0.55953847, 0.56550496, 0.57144065, 0.5773304 , 0.58314243, 0.58884039,
    0.59439859, 0.59982752, 0.60514774, 0.6103857 , 0.615554  , 0.62063104,
    0.62557536, 0.63034712, 0.63491583, 0.63926669, 0.64340524, 0.64733751,
    0.65106447, 0.65457135, 0.65782649, 0.660798
  };
  double Z[] = {
     0.        , -0.00946427, -0.01892723, -0.02839665, -0.03750362, -0.0469502 ,
    -0.05641336, -0.06587392, -0.07532046, -0.08474134, -0.09413757, -0.10352535,
    -0.11292   , -0.12232572, -0.13173255, -0.1411269 , -0.15049601, -0.15982796,
    -0.16911406, -0.17834574, -0.1875213 , -0.19666013, -0.20579732, -0.2149517 ,
    -0.22410942, -0.23323011, -0.24226935, -0.25119166, -0.25998504, -0.26864154,
    -0.27717208, -0.28562684, -0.29405864, -0.30242558, -0.31057111, -0.31832703,
    -0.32566422, -0.33265586, -0.33931711, -0.34551163, -0.35109578, -0.35603769,
    -0.36040236, -0.36420726, -0.3673851 , -0.36987789, -0.37175562, -0.37317612,
    -0.37422303, -0.37485444, -0.37497819, -0.37458894, -0.3737709 , -0.3726319 ,
    -0.3712027 , -0.36941375, -0.36718013, -0.36448535, -0.36140513, -0.35802853,
    -0.35438623, -0.35042529, -0.34608086, -0.34133856, -0.33625055, -0.33088089,
    -0.32528134, -0.31948688, -0.31352859, -0.30742686, -0.30118433, -0.29478343,
    -0.28820777, -0.2814585 , -0.27455265, -0.26751284, -0.26035895, -0.25311462,
    -0.24580498, -0.2384586 , -0.23108732, -0.22367928, -0.21621011, -0.20865358,
    -0.20099365, -0.19324157, -0.18541447, -0.17753209, -0.16960385, -0.16161687,
    -0.15354707, -0.14537406, -0.13708584, -0.12868123, -0.12017003, -0.11156161,
    -0.10286232, -0.09407201, -0.0851854 , -0.0762    
  };
  vvi->num_points = npts;
  for (int i=0; i<npts; i++) {
    vvi->t[i] = t[i];
    vvi->R[i] = R[i];
    vvi->Z[i] = Z[i];
  }
}

// Z is constant at -0.8
// R goes from 1.5 to 1.75
void pfunc_upper(double s, double* RZ){
//  // Vertical plate.
//  RZ[0] = 0.14047;
//  RZ[1] = -(s-0.061)*0.6;
//  // Tilted plate.
//  double p0[2] = {0.14047, 0.0};
//  double p1[2] = {0.182, -0.327};
//  RZ[0] = (1-s)*p0[0]+s*p1[0];
//  RZ[1] = (1-s)*p0[1]+s*p1[1];
  // Parametrized curve.
  struct ltx_lower_vessel_info vvi;
  get_ltx_lower_vessel_info(&vvi);

  // Find indices in t that bound s.
  int idx_tlo, idx_tup;
  if (s < 1e-8) {
    idx_tlo = 0;
    idx_tup = 0;
  }
  else if (fabs(s-1.0) < 1e-8) {
    idx_tlo = vvi.num_points-1;
    idx_tup = vvi.num_points-1;
  }
  else {
    for (int i=0; i<vvi.num_points-1; i++) {
      if (vvi.t[i] <= s && s < vvi.t[i+1]) {
        idx_tlo = i;
        idx_tup = i+1;
        break;
      }
    }
  }
  // Interpolate the value of R and Z.
  double Dt = vvi.t[idx_tup]-vvi.t[idx_tlo];
  if (idx_tlo == idx_tup) {
    RZ[0] = vvi.R[idx_tlo];
    RZ[1] = vvi.Z[idx_tlo];
  }
  else {
    RZ[0] = ((s-vvi.t[idx_tlo])/Dt)*vvi.R[idx_tup] + ((vvi.t[idx_tup]-s)/Dt)*vvi.R[idx_tlo];
    RZ[1] = ((s-vvi.t[idx_tlo])/Dt)*vvi.Z[idx_tup] + ((vvi.t[idx_tup]-s)/Dt)*vvi.Z[idx_tlo];
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
  struct ltx_lower_vessel_info vvi;
  get_ltx_lower_vessel_info(&vvi);
  // Reflect across midplane.
  for (int i=0; i<vvi.num_points; i++)
    vvi.Z[i] = -vvi.Z[i];

  // Find indices in t that bound s.
  int idx_tlo, idx_tup;
  if (s < 1e-8) {
    idx_tlo = 0;
    idx_tup = 0;
  }
  else if (fabs(s-1.0) < 1e-8) {
    idx_tlo = vvi.num_points-1;
    idx_tup = vvi.num_points-1;
  }
  else {
    for (int i=0; i<vvi.num_points-1; i++) {
      if (vvi.t[i] <= s && s < vvi.t[i+1]) {
        idx_tlo = i;
        idx_tup = i+1;
        break;
      }
    }
  }
  // Interpolate the value of R and Z.
  double Dt = vvi.t[idx_tup]-vvi.t[idx_tlo];
  if (idx_tlo == idx_tup) {
    RZ[0] = vvi.R[idx_tlo];
    RZ[1] = vvi.Z[idx_tlo];
  }
  else {
    RZ[0] = ((s-vvi.t[idx_tlo])/Dt)*vvi.R[idx_tup] + ((vvi.t[idx_tup]-s)/Dt)*vvi.R[idx_tlo];
    RZ[1] = ((s-vvi.t[idx_tlo])/Dt)*vvi.Z[idx_tup] + ((vvi.t[idx_tup]-s)/Dt)*vvi.Z[idx_tlo];
  }
}

double psi_N(double psi, double psi_lcfs, double psi_axis)
{
  // Normalized psi.
  //   - psi: poloidal flux.
  //   - psi_lcfs: poloidal flux at the LCFS.
  //   - psi_axis: poloidal flux at the exis.
  return (psi - psi_axis) / (psi_lcfs - psi_axis);
}

double psi_psi_N(double psi_N, double psi_lcfs, double psi_axis)
{
  // Compute psi given the normalized psi.
  //   - psi_N: normalized poloidal flux.
  //   - psi_lcfs: poloidal flux at the LCFS.
  //   - psi_axis: poloidal flux at the exis.
  return psi_N * (psi_lcfs - psi_axis) + psi_axis;
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

void density_init_neut(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0], z = xn[1];

  struct gk_ltx_ctx *app = ctx;
  double n0_neut = app->n0_neut;
  double Lx = app->Lx;
  double Lz = app->Lz;
  double x_LCFS = app->psi_LCFS;

  double n_x = 0.5*(1+tanh(-(x-0.9*x_LCFS)/(0.4*Lx)));

  double n_z = 0.0;
  if (z <= 0.0)
    n_z = (pow(1.0/cosh(z+Lz/2.),2.0) + 1.e-6);
  else
    n_z = (pow(1.0/cosh(z-Lz/2.),2.0) + 1.e-6);

  fout[0] = 10.0*n0_neut*n_x*n_z;
}

void temp_init_neut(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[1];

  struct gk_ltx_ctx *app = ctx;
  double T0_neut = app->T0_neut;
  double Lx = app->Lx;
  double x_LCFS = app->psi_LCFS;

  fout[0] = T0_neut*(0.5*(1.0-tanh(-(x-0.8*x_LCFS)/(0.2*Lx)))+0.1);

}

void udrift_init_neut(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
}

// Profiles used in the inner radial ghost cell.

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

struct gk_ltx_ctx
create_ctx(void)
{
  int cdim = 2, vdim = 2; // Dimensionality.

  // Universal constant parameters.
  double eps0 = GKYL_EPSILON0, eV = GKYL_ELEMENTARY_CHARGE;
  double mp = GKYL_PROTON_MASS, me = GKYL_ELECTRON_MASS;
  double qi = eV; // ion charge
  double qe = -eV; // electron charge

  char geqdsk_file[128] = "../../../experiment/li_863mg_103955_04/LTX_103955_04.eqdsk";
  // Get the LCFS and axis psi.
  struct gkyl_efit_inp efit_inp = {
    // psiRZ and related inputs
    .rz_poly_order = 2,
    .flux_poly_order = 1,
  };
  // Copy eqdsk file into efit_inp.
  memcpy(efit_inp.filepath, geqdsk_file, sizeof(geqdsk_file));
  struct gkyl_efit *efit = gkyl_efit_new(&efit_inp);
  double psi_LCFS = efit->sibry;
  double psi_axis = efit->simag;
  gkyl_efit_release(efit);

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
  double psi_N_min = 1.256450306e+00; // Minimum psi_N.
  double psi_N_max = 5.725828231e-01; // Maximum psi_N.
  double psi_min = psi_psi_N(psi_N_min, psi_LCFS, psi_axis); // Minimum psi.
  double psi_max = psi_psi_N(psi_N_max, psi_LCFS, psi_axis); // Maximum psi.  
  printf("Radial grid:\n");
  printf("  Nx_core          = %3d\n",Nx_core);
  printf("  Nx_sol           = %3d\n",Nx_sol);
  printf("  Nx               = %3d\n",Nx);
  printf("  psi_LCFS         = %.9e\n",psi_LCFS        );
  printf("  psi_axis         = %.9e\n",psi_axis        );
  printf("  psi_N_min        = %.9e\n",psi_N_min       );
  printf("  psi_N_max        = %.9e\n",psi_N_max       );
  printf("  psi_min          = %.9e\n",psi_min         );
  printf("  psi_max          = %.9e\n",psi_max         );

  // Adjust psi so that psi_LCFS is at a cell boundary.
  double delta_psi = (psi_max-psi_LCFS)/Nx_core; // Cell length in psi.
  psi_min = psi_LCFS - Nx_sol*delta_psi; // Outer SOL boundary (psi increases towards core).
  psi_N_min = psi_N(psi_min, psi_LCFS, psi_axis); // Minimum psi_N.
  printf("After adjustment:\n");
  printf("  psi_N_min        = %.9e\n",psi_N_min     );
  printf("  psi_N_max        = %.9e\n",psi_N_max     );
  printf("  psi_min          = %.9e\n",psi_min       );
  printf("  psi_max          = %.9e\n",psi_max       );

  double Lz = 2.0*(M_PI-1e-14);    // Domain size along magnetic field.
  double Lx = psi_max - psi_min;
  double Lx_core = psi_LCFS - psi_min;
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

  double rec_frac = 0.6; // Recycling coefficient.
  double neut_gamma = 5.0/3.0; // Neutral adiabatic index.
  double m_neut = mi; // Neutral mass.
  double n0_neut = pow(10,17.25); // Reference neutral density.
  double T0_neut = 0.5*Ti0; // Reference neutral temperature.

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
  double vpar_max_ion_c = 1.0/sqrt(2.0);
  double mu_min_ion_c = 0.;
  double mu_max_ion_c = 1.;
  double vpar_min_elc_c = -1.0/sqrt(2.0);
  double vpar_max_elc_c = 1.0/sqrt(2.0);
  double mu_min_elc_c = 0.;
  double mu_max_elc_c = 1.;

  double t_end = 8*500.0e-6;
  int num_frames = 80;
  double write_phase_freq = 1.0; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct gk_ltx_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,
    .Lz     = Lz    ,
    .z_min = -Lz/2.,  .z_max = Lz/2.,
    .psi_min = psi_min,  .psi_max = psi_max,
    .psi_LCFS = psi_LCFS,
    .Lx = Lx,
    .Lx_core = Lx_core,
  
    .me = me,  .qe = qe,
    .mi = mi,  .qi = qi,
    .n0 = n0,  .Te0 = Te0,  .Ti0 = Ti0,
  
    .rec_frac   = rec_frac  ,
    .neut_gamma = neut_gamma,
    .m_neut     = m_neut    ,
    .n0_neut    = n0_neut   ,
    .T0_neut    = T0_neut   ,

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

    .t_end = t_end,  .num_frames = num_frames,
    .write_phase_freq = write_phase_freq,
    .int_diag_calc_num = int_diag_calc_num,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };

  // Copy eqdsk file into ctx.
  memcpy(ctx.geqdsk_file, geqdsk_file, sizeof(geqdsk_file));

  return ctx;
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

    .react_neut = {
      .num_react = 2,
      .react_type = {
        { .react_id = GKYL_REACT_IZ,
          .type_self = GKYL_SELF_ELC,
          .ion_id = GKYL_ION_H,
          .elc_nm = "elc",
          .ion_nm = "ion",
          .donor_nm = "H0",
          .charge_state = 0,
          .ion_mass = ctx.mi,
          .elc_mass = ctx.me,
        },
        { .react_id = GKYL_REACT_RECOMB,
          .type_self = GKYL_SELF_ELC,
          .ion_id = GKYL_ION_H,
          .elc_nm = "elc",
          .ion_nm = "ion",
          .recvr_nm = "H0",
          .charge_state = 0,
          .ion_mass = ctx.mi,
          .elc_mass = ctx.me,
        },
      },
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

    .react_neut = {
      .num_react = 3,
      .react_type = {
        { .react_id = GKYL_REACT_IZ,
          .type_self = GKYL_SELF_ION,
          .ion_id = GKYL_ION_H,
          .elc_nm = "elc",
          .ion_nm = "ion",
          .donor_nm = "H0",
          .charge_state = 0,
          .ion_mass = ctx.mi,
          .elc_mass = ctx.me,
        },
        { .react_id = GKYL_REACT_RECOMB,
          .type_self = GKYL_SELF_ION,
          .ion_id = GKYL_ION_H,
          .elc_nm = "elc",
          .ion_nm = "ion",
          .recvr_nm = "H0",
          .charge_state = 0,
          .ion_mass = ctx.mi,
          .elc_mass = ctx.me,
        },
        { .react_id = GKYL_REACT_CX,
          .type_self = GKYL_SELF_ION,
          .ion_id = GKYL_ION_H,
          .ion_nm = "ion",
          .partner_nm = "H0",
          .ion_mass = ctx.mi,
          .partner_mass = ctx.mi,
        },
      },
      .write_diagnostics = true,
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

  // Neutrals
  struct gkyl_gyrokinetic_neut_species H0 = {
    .name = "H0", .mass = ctx.m_neut,
    .gas_gamma = ctx.neut_gamma,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = density_init_neut,
      .ctx_upar = &ctx,
      .udrift = udrift_init_neut,
      .ctx_temp = &ctx,
      .temp = temp_init_neut,
    },

    .scaling = {
      .type = GKYL_GK_SPECIES_SCALING_RECYCLING_IZ_BALANCE,
      .impacting_ion_name =  "ion" ,
      .impacting_ion_id = GKYL_ION_H,
      .electron_name = "elc",
      .recycling_coeff = ctx.rec_frac,
      .num_boundaries = 2,
      .boundaries_dir = {1, 1,},
      .boundaries_edge = {GKYL_LOWER_EDGE, GKYL_UPPER_EDGE,},
      .write_diagnostics = true,
    },

    .num_diag_moments = 1,
    .diag_moments = { GKYL_F_MOMENT_LTE, },
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
    .rz_poly_order = 2,
    .flux_poly_order = 1,
    .reflect = true,
  };
  // Copy eqdsk file into efit_inp.
  memcpy(efit_inp.filepath, ctx.geqdsk_file, sizeof(ctx.geqdsk_file));

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

    .num_neut_species = 1,
    .neut_species = { H0 },

    .field = field,

    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0], app_args.cuts[1] },
      .comm = comm,
    },
  };

  struct gkyl_gyrokinetic_run_inp run_inp = {
    .app_inp = app_inp,
    .time_stepping = {
      .t_end = ctx.t_end,
      .num_frames = ctx.num_frames,
      .write_phase_freq = ctx.write_phase_freq,
      .int_diag_calc_num = ctx.int_diag_calc_num,
      .dt_failure_tol = ctx.dt_failure_tol,
      .num_failures_max = ctx.num_failures_max,
      .is_restart = app_args.is_restart,
      .restart_frame = app_args.restart_frame,
      .num_steps = app_args.num_steps,
    },
    .print_verbosity = {
      .enabled = true,
      .frequency = 0.1,
    },
  };

  gkyl_gyrokinetic_run_simulation(&run_inp);
  
  gkyl_gyrokinetic_comms_release(comm);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Finalize();
#endif

  return 0;
}
