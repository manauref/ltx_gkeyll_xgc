#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_efit.h>
#include <gkyl_gyrokinetic_multib.h>
#include <gkyl_gyrokinetic_run.h>
#include <gkyl_mpi_comm.h>
#include <gkyl_null_comm.h>
#include <gkyl_tok_geo.h>

#include <rt_arg_parse.h>

struct gk_ltx_ctx {
  int cdim, vdim; // Dimensionality.

  char geqdsk_file[128]; // File with equilibrium.

  // Geometry and magnetic field.
  int num_blocks; // Number of blocks.
  double psi_axis; // Psi at the magnetic axis.
  double psi_LCFS; // Radial location of the last closed flux surface.
  double psi_lower_core, psi_upper_core; // Core psi extents.
  double psi_lower_sol, psi_upper_sol; // SOL psi extents.
  double theta_lower, theta_upper; // Extents in parallel direction.
  double Lpsi_core, Lpsi_sol, Lpsi; // Box size in psi.
  double Ltheta; // Box size in theta.

  double charge_elc; // Electron charge.
  double charge_ion; // Ion charge.
  double mass_elc; // Electron mass.
  double mass_ion; // Ion mass.
  double B0; // Magnetic field.
  double n0; // Density.
  double Te0; // Electron temperature
  double Ti0; // Ion temperature.

  double nu_frac; // Fraction multiplying collision frequency.

  double anom_diff_D; // Anomalous diffusion coefficient.

  // Source parameters.
  double ndot_srcOMP;      // Amplitude of the OMP source
  double x_srcOMP;        // Radial location of the OMP source.
  double Te_srcOMP;       // Te for the OMP source.
  double Ti_srcOMP;       // Ti for the OMP source.
  double sigma_srcOMP;    // Radial spread of the OMP source.
  double floor_src;       // Source floor.

  // Grid.
  int Npsi_sol; // Number of cells in psi in the SOL.
  int Npsi_core; // Number of cells in psi in the core.
  int Ntheta; // Number of cells in theta.
  int Nvpar, Nmu; // Number of cells in vpar,mu.

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

void sol_limiter_func_upper(double s, double* RZ){
//  // Vertical plate.
//  RZ[0] = 0.14047;
//  RZ[1] = -(s-0.061)*0.6;
//  // Tilted plate.
//  double p0[2] = {0.14047, 0.0};
//  double p1[2] = {0.182, -0.327};
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

void sol_limiter_func_lower(double s, double* RZ){
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

struct gkyl_gk_block_geom*
create_ltx_iwl_gk_block_geom(void *ctx)
{
  struct gk_ltx_ctx *params = ctx;

  struct gkyl_gk_block_geom *bgeom = gkyl_gk_block_geom_new(params->cdim, params->num_blocks);

  /* Block layout and coordinates.

    z  
    ^  
    |

    |  +-----------------+-----------------+
    |  |                 |                 |
    |  |                 |                 |
    |  | b0              | b1              |
    |  + core            + sol             |
    |  |                 |                 |
    |  |                 |                 |
    |  +-----------------+-----------------+

    0 -------------------------------------------------------- -> x

  */  

  struct gkyl_efit_inp efit_inp = {
    // psiRZ and related inputs
    .rz_poly_order = 2,
    .flux_poly_order = 1,
    .reflect = true,
  };
  // Copy eqdsk file into efit_inp.
  memcpy(efit_inp.filepath, params->geqdsk_file, sizeof(params->geqdsk_file));

  // Theta limits.
  double theta_lower = params->theta_lower, theta_upper = params->theta_upper;

  // Psi limits.
  double psi_lower_core = params->psi_lower_core;
  double psi_upper_core = params->psi_upper_core;
  double psi_lower_sol  = params->psi_lower_sol ;
  double psi_upper_sol  = params->psi_upper_sol ;
  double psi_LCFS       = params->psi_LCFS ;

  // Number of cells.
  int Npsi_sol  = params->Npsi_sol ;
  int Npsi_core = params->Npsi_core;
  int Ntheta    = params->Ntheta   ;

  // Block 0: Core region.
  gkyl_gk_block_geom_set_block(bgeom, 0, &(struct gkyl_gk_block_geom_info) {
      .lower = { psi_lower_core, theta_lower },
      .upper = { psi_upper_core, theta_upper },
      .cells = { Npsi_core, Ntheta },
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_GEOMETRY_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_GEOMETRY_TOKAMAK_IWL,
          .rclose = 0.7,
          .rleft  = 0.0,
          .rright = 0.7,
          .rmin   = 0.1,
          .rmax   = 0.7,
          .zmin   = -0.3675,
          .zmax   =  0.3675,
        },
        .has_LCFS = true,
        .x_LCFS = psi_LCFS, // Location of last closed flux surface.
      },

      .connections[0] = { // x-direction.
        { .bid = 1, .dir = 0, .edge = GKYL_UPPER_POSITIVE },
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL },  // Physical boundary.
      },
      .connections[1] = { // z-direction.
        { .bid = 0, .dir = 1, .edge = GKYL_UPPER_POSITIVE }, // Periodic boundary.
        { .bid = 0, .dir = 1, .edge = GKYL_LOWER_POSITIVE }, // Periodic boundary.
      }
    }
  );

  // Block 1: SOL.
  gkyl_gk_block_geom_set_block(bgeom, 1, &(struct gkyl_gk_block_geom_info) {
      .lower = { psi_lower_sol, theta_lower },
      .upper = { psi_upper_sol, theta_upper },
      .cells = { Npsi_sol, Ntheta },
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_GEOMETRY_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_GEOMETRY_TOKAMAK_IWL,
          .rclose = 0.7,
          .rleft  = 0.1,
          .rright = 0.7,
          .rmin   = 0.1,
          .rmax   = 0.7,
          .zmin   = -0.3675,
          .zmax   =  0.3675,
          .plate_spec = true,
          .plate_func_lower = sol_limiter_func_lower,
          .plate_func_upper = sol_limiter_func_upper,
        },
        .has_LCFS = true,
        .x_LCFS = psi_LCFS, // Location of last closed flux surface.
      },
      
      .connections[0] = { // x-direction.
        { .bid = 1, .dir = 0, .edge = GKYL_PHYSICAL },  // Physical boundary.
        { .bid = 0, .dir = 0, .edge = GKYL_LOWER_POSITIVE },
      },
      .connections[1] = { // z-direction.
        { .bid = 1, .dir = 1, .edge = GKYL_PHYSICAL }, // Physical boundary.
        { .bid = 1, .dir = 1, .edge = GKYL_PHYSICAL }, // Physical boundary.
      }
    }
  );

  return bgeom;
}

// Source profiles.
void density_srcOMP(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[1];

  struct gk_ltx_ctx *app = ctx;
  double ndot_srcOMP = app->ndot_srcOMP;
  double x_srcOMP = app->x_srcOMP;
  double sigma_srcOMP = app->sigma_srcOMP;
  double floor_src = app->floor_src;

  fout[0] = ndot_srcOMP * fmax(0.5*(1.0-tanh(-(x-x_srcOMP)/(2.0*sigma_srcOMP))), floor_src);
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
  double psi_max = app->psi_upper_core;
  double Lpsi = app->Lpsi;
  double n0 = app->n0;

  fout[0] = 0.8*n0*(1.5 - tanh((9.0/Lpsi)*(psi_max-0.646767778*Lpsi-x)));
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
  double psi_max = app->psi_upper_core;
  double Lpsi = app->Lpsi;
  double Ti0 = app->Ti0;

  fout[0] = 0.8*Ti0*(1.1 - tanh((9.0/Lpsi)*(psi_max-0.646767778*Lpsi-x)));
}

void temp_init_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[1];

  struct gk_ltx_ctx *app = ctx;
  double psi_max = app->psi_upper_core;
  double Lpsi = app->Lpsi;
  double Te0 = app->Te0;

  fout[0] = 0.8*Te0*(1.1 - tanh((9.0/Lpsi)*(psi_max-0.646767778*Lpsi-x)));
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
  // Get parameters from eqdsk file.
  struct gkyl_efit_inp efit_inp = {
    // psiRZ and related inputs.
    .rz_poly_order = 2,
    .flux_poly_order = 1,
  };
  // Copy eqdsk file into efit_inp.
  memcpy(efit_inp.filepath, geqdsk_file, sizeof(geqdsk_file));
  struct gkyl_efit *efit = gkyl_efit_new(&efit_inp);
  double psi_LCFS = efit->sibry; // psi at LCFS (from eqdsk).
  double psi_axis = efit->simag; // psi at the magnetic axis.
  gkyl_efit_release(efit);

  double num_blocks = 2;

  // Grid parameters
  int Npsi = 16;
  double Npsi_core_frac = (5.0/8.0);
  int Npsi_core = floor(Npsi_core_frac*Npsi);
  int Npsi_sol = Npsi-Npsi_core;
  int Ntheta = 16;
  int Nvpar = 12;
  int Nmu = 8;
  int poly_order = 1;

  // Assume psi increases towards the core (depends on eqdsk).
  double psi_N_min = 1.256450306e+00; // Minimum psi_N.
  double psi_N_max = 5.725828231e-01; // Maximum psi_N.
  double psi_min = psi_psi_N(psi_N_min, psi_LCFS, psi_axis); // Minimum psi.
  double psi_max = psi_psi_N(psi_N_max, psi_LCFS, psi_axis); // Maximum psi.
  printf("Radial grid:\n");
  printf("  Npsi_core        = %3d\n",Npsi_core);
  printf("  Npsi_sol         = %3d\n",Npsi_sol);
  printf("  Npsi             = %3d\n",Npsi);
  printf("  psi_LCFS         = %.9e\n",psi_LCFS        );
  printf("  psi_axis         = %.9e\n",psi_axis        );
  printf("  psi_N_min        = %.9e\n",psi_N_min       );
  printf("  psi_N_max        = %.9e\n",psi_N_max       );
  printf("  psi_min          = %.9e\n",psi_min         );
  printf("  psi_max          = %.9e\n",psi_max         );

  // Adjust psi so that psi_LCFS is at a cell boundary.
  double delta_psi = (psi_max-psi_LCFS)/Npsi_core; // Cell length in psi.
  psi_min = psi_LCFS - Npsi_sol*delta_psi; // Outer SOL boundary (psi increases towards core).
  psi_N_min = psi_N(psi_min, psi_LCFS, psi_axis); // Minimum psi_N.
  printf("After adjustment:\n");
  printf("  psi_N_min        = %.9e\n",psi_N_min     );
  printf("  psi_N_max        = %.9e\n",psi_N_max     );
  printf("  psi_min          = %.9e\n",psi_min       );
  printf("  psi_max          = %.9e\n",psi_max       );
                            
  double psi_lower_core = psi_LCFS;
  double psi_upper_core = psi_max;

  double psi_lower_sol = psi_min;
  double psi_upper_sol = psi_LCFS;

  double theta_lower = -(M_PI-1e-14);
  double theta_upper =  (M_PI-1e-14);

  double Lpsi_core = psi_upper_core - psi_lower_core;
  double Lpsi_sol = psi_upper_sol - psi_lower_sol;
  double Lpsi = Lpsi_sol + Lpsi_sol;
  double Ltheta = theta_upper - theta_lower;

  printf("  psi_lower_core = %.13e\n",psi_lower_core);
  printf("  psi_upper_core = %.13e\n",psi_upper_core);
  printf("  psi_lower_sol  = %.13e\n",psi_lower_sol);
  printf("  psi_upper_sol  = %.13e\n",psi_upper_sol);
  printf("  Ntheta         = %d\n",Ntheta     );

  // Plasma parameters.
  double mi = mp;   // Hydrogen ions.
  double Te_max = 320*eV; // [J].
  double Ti_max = Te_max/2.0; // [J].
  double n_max  = 1.6e19; // [1/m^3].
  double bmag_min = 1.666901e-01; // [T].
  double bmag_max = 6.814593e-01; // [T].

  double n0  = 0.5*n_max;
  double Te0 = 0.5*Te_max;
  double Ti0 = 0.5*Ti_max;
  double B0  = 0.5*(bmag_min+bmag_max);

  double nu_frac = 1.0;
  // Electron-electron collision freq.
  double logLambda_elc = 6.6 - 0.5 * log(n0/1e20) + 1.5 * log(Ti0/eV);
  double nu_elc = nu_frac * logLambda_elc * pow(eV, 4) * n0 /
    (6*sqrt(2.) * pow(M_PI,3./2.) * pow(eps0,2) * sqrt(me) * pow(Te0,3./2.));
  // Ion-ion collision freq.
  double logLambda_ion = 6.6 - 0.5 * log(n0/1e20) + 1.5 * log(Ti0/eV);
  double nu_ion = nu_frac * logLambda_ion * pow(eV, 4) * n0 /
    (12 * pow(M_PI,3./2.) * pow(eps0,2) * sqrt(mi) * pow(Ti0,3./2.));
  double nu_elc_ion = nu_elc*sqrt(2.0);
  double nu_ion_elc = nu_elc_ion*(me/mi);

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
  double ndot_srcOMP = n_LCFS*cs_LCFS/Lc_LCFS;
  double x_srcOMP = psi_max-0.15*Lpsi;
  double Te_srcOMP = Te_max;
  double Ti_srcOMP = Ti_max;
  double sigma_srcOMP = 0.03*Lpsi;
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
  double write_phase_freq = 0.2; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.
                             //
  struct gk_ltx_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,
    // Reference parameters.
    .charge_elc = qe, 
    .charge_ion = qi, 
    .mass_elc = me, 
    .mass_ion = mi,
    .B0 = B0, 
    .n0 = n0, 
    .Te0 = Te0, 
    .Ti0 = Ti0, 

    // Geometry.
    .num_blocks = num_blocks,
    .psi_axis = psi_axis,
    .psi_LCFS = psi_LCFS,
    .psi_lower_core = psi_lower_core,
    .psi_upper_core = psi_upper_core,
    .psi_lower_sol = psi_lower_sol,
    .psi_upper_sol = psi_upper_sol,
    .theta_lower = theta_lower,
    .theta_upper = theta_upper,
    .Lpsi_core = Lpsi_core, 
    .Lpsi_sol = Lpsi_sol, 
    .Lpsi = Lpsi, 
    .Ltheta = Ltheta,
    .Npsi_sol = Npsi_sol,
    .Npsi_core = Npsi_core,
    .Ntheta = Ntheta,

    .nu_frac = nu_frac,

    .anom_diff_D = anom_diff_D,

    // Source parameters.
    .ndot_srcOMP  = ndot_srcOMP ,
    .Te_srcOMP    = Te_srcOMP   ,
    .Ti_srcOMP    = Ti_srcOMP   ,
    .x_srcOMP     = x_srcOMP    ,
    .sigma_srcOMP = sigma_srcOMP,
    .floor_src    = floor_src   ,

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
    // Velocity space cells.
    .Nvpar = Nvpar,
    .Nmu = Nmu,

    // Time stepping and I/O.
    .t_end = t_end, 
    .num_frames = num_frames, 
    .write_phase_freq = write_phase_freq,
    .int_diag_calc_num = int_diag_calc_num,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };

  // Copy eqdsk file into ctx.
  memcpy(ctx.geqdsk_file, geqdsk_file, sizeof(geqdsk_file));
  return ctx;
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Init(&argc, &argv);
  }
#endif

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  // Construct communicator for use in app.
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);

  struct gk_ltx_ctx ctx = create_ctx(); // Context for init functions.
                    
  // Construct block geometry
  struct gkyl_gk_block_geom *bgeom = create_ltx_iwl_gk_block_geom(&ctx);

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

  struct gkyl_gyrokinetic_bc elc_phys_bcs[] = {
    // block 0 BCs
    { .bidx = 0, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_FIXED_FUNC, .projection = elc_bcx_lo,},
    // block 1 BCs
    { .bidx = 1, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 1, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},
    { .bidx = 1, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},
  };

 struct gkyl_gyrokinetic_multib_species_pb mb_elc = {
    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = density_init,
      .ctx_upar = &ctx,
      .upar = upar_init,
      .ctx_temp = &ctx,
      .temp = temp_init_ion,
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

  };

  struct gkyl_gyrokinetic_multib_species_pb elc_blocks[ctx.num_blocks];
  mb_elc.block_id = 0;
  elc_blocks[0] = mb_elc;
  mb_elc.block_id = 1;
  elc_blocks[1] = mb_elc;

  struct gkyl_gyrokinetic_multib_species elc = {
    .name = "elc",
    .charge = ctx.charge_elc, .mass = ctx.mass_elc,
    .vdim = ctx.vdim,
    .lower = { ctx.vpar_min_elc_c, ctx.mu_min_elc_c},
    .upper = { ctx.vpar_max_elc_c, ctx.mu_max_elc_c},
    .cells = { ctx.Nvpar, ctx.Nmu },

    .mapc2p = {
      .mapping = mapc2p_vel_elc,
      .ctx = &ctx,
    },

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .den_ref = ctx.n0, // Density used to calculate coulomb logarithm.
      .temp_ref = ctx.Te0, // Temperature used to calculate coulomb logarithm.
      .num_cross_collisions = 1,
      .collide_with = { "ion" },
    },

    .anomalous_diffusion = {
      .anomalous_diff_id = GKYL_GK_ANOMALOUS_DIFF_D,
      .D_profile = diffusion_D_func,
      .D_profile_ctx = &ctx,
//      .write_diagnostics = true,
    },

    .num_physical_bcs = sizeof(elc_phys_bcs)/sizeof(struct gkyl_gyrokinetic_bc),
    .bcs = elc_phys_bcs,

    .blocks = elc_blocks,
    .duplicate_across_blocks = false,

    .num_diag_moments = 4,
    .diag_moments = { GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_BIMAXWELLIAN },
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    .boundary_flux_diagnostics = {
      .num_integrated_diag_moments = 1,
      .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
//      .time_integrated = true,
    },
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

  struct gkyl_gyrokinetic_bc ion_phys_bcs[] = {
    // block 0 BCs
    { .bidx = 0, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_FIXED_FUNC, .projection = ion_bcx_lo,},
    // block 1 BCs
    { .bidx = 1, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 1, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},
    { .bidx = 1, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},
  };

  struct gkyl_gyrokinetic_multib_species_pb mb_ion = {
    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = density_init,
      .ctx_upar = &ctx,
      .upar= upar_init,
      .ctx_temp = &ctx,
      .temp = temp_init_ion,
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

  };

  struct gkyl_gyrokinetic_multib_species_pb ion_blocks[ctx.num_blocks];
  mb_ion.block_id = 0;
  ion_blocks[0] = mb_ion;
  mb_ion.block_id = 1;
  ion_blocks[1] = mb_ion;

  struct gkyl_gyrokinetic_multib_species ion = {
    .name = "ion",
    .charge = ctx.charge_ion, .mass = ctx.mass_ion,
    .vdim = ctx.vdim,
    .lower = { ctx.vpar_min_ion_c, ctx.mu_min_ion_c},
    .upper = { ctx.vpar_max_ion_c, ctx.mu_max_ion_c},
    .cells = { ctx.Nvpar, ctx.Nmu },

    .mapc2p = {
      .mapping = mapc2p_vel_ion,
      .ctx = &ctx,
    },

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .den_ref = ctx.n0, // Density used to calculate coulomb logarithm.
      .temp_ref = ctx.Ti0, // Temperature used to calculate coulomb logarithm.
      .bmag_ref = ctx.B0,
      .num_cross_collisions = 1,
      .collide_with = { "elc" },
    },

    .anomalous_diffusion = {
      .anomalous_diff_id = GKYL_GK_ANOMALOUS_DIFF_D,
      .D_profile = diffusion_D_func,
      .D_profile_ctx = &ctx,
//      .write_diagnostics = true,
    },

    .num_physical_bcs = sizeof(ion_phys_bcs)/sizeof(struct gkyl_gyrokinetic_bc),
    .bcs = ion_phys_bcs,

    .blocks = ion_blocks,
    .duplicate_across_blocks = false,

    .num_diag_moments = 4,
    .diag_moments = { GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_BIMAXWELLIAN },
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    .boundary_flux_diagnostics = {
      .num_integrated_diag_moments = 1,
      .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
//      .time_integrated = true,
    },
  };

  // Field object.
  struct gkyl_gyrokinetic_multib_field_pb field_blocks[1];
  field_blocks[0] = (struct gkyl_gyrokinetic_multib_field_pb) {
    // No block specific field info for this simulation
  };

  struct gkyl_gyrokinetic_bc field_phys_bcs[] = {
    { .bidx = 0, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
    { .bidx = 1, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
  };

  struct gkyl_gyrokinetic_multib_field field = {
    .blocks = field_blocks, 
    .duplicate_across_blocks = true,

    .num_physical_bcs = sizeof(field_phys_bcs)/sizeof(struct gkyl_gyrokinetic_bc),
    .bcs = field_phys_bcs,
  };

  struct gkyl_gyrokinetic_multib app_inp = {
    .name = "gk_multib_ltx_iwl_2x2v_p1",

    .cdim = ctx.cdim,
    .poly_order = 1,
    .basis_type = app_args.basis_type,
    .cfl_frac = 1.0,

    .gk_block_geom = bgeom,
    
    .num_species = 2,
    .species = { elc, ion},

    .field = field,

    .comm = comm,
    .use_gpu = app_args.use_gpu,
  };

  struct gkyl_gyrokinetic_run_inp run_inp = {
    .app_type = GKYL_GK_MULTIB,
    .multib_app_inp = app_inp,
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
    .print_verbosity = true,
  };

  gkyl_gyrokinetic_run_simulation(&run_inp);

  gkyl_gk_block_geom_release(bgeom);
  gkyl_gyrokinetic_comms_release(comm);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Finalize();
#endif
  
  return 0;
}
