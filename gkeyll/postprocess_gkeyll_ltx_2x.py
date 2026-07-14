#[ ........................................................... ]#
#[
#[ post processing for 2x LTX simulation.
#[
#[ Manaure Francisquez.
#[ February 2026.
#[
#[ ........................................................... ]#

import numpy as np
import postgkyl as pg
import click #[ For using pgkyl commands.
from postgkyl.pgkyl import cli #[ For using pgkyl commands.
import matplotlib.pyplot as plt
import sys #[ For error exit.
import h5py #[ For saving reduced data.
import pgkylUtil as pgu #[ Some postgkyl wrappers.
#[ Append path to utilities folder.
sys.path.insert(0, '../util/')
import ltx_common_util as lcu

#[ Plotting options.
plot_grid_RZ          = False  #[ Grid on RZ plane.
plot_vs_x             = False  #[ A quantity at the outboard midplane.
plot_nT_vs_x          = False  #[ Sensity and temperature profiles vs. x.
plot_src_mom_vs_x     = False  #[ Source moments vs. x.
plot_src_int_mom_vs_t = False  #[ Source integrated moments vs. t.
plot_vs_RZ            = True  #[ A quantity at the R-Z midplane.
plot_nT_vs_x_multisim = False  #[ Density and temperature profiles vs. x for multiple sims.
plot_nu_RZ            = False  #[ Collision frequency on RZ plane.
plot_nu_vs_x          = False  #[ Collision frequency vs. x.
plot_int_mom_vs_t     = False   #[ Integrated moments vs. t.

out_data_dir  = './data/'
out_fig_dir   = './figures/'
output_prefix = 'ltx_gkeyll_'

save_data          = True    #[ Indicate whether to save data in plot to HDF5 file.
out_figure_file    = True     #[ Output a figure file?.
figure_file_format = '.png'    #[ Can be .png, .pdf, .ps, .eps, .svg.

sim_name   = 'gk_ltx_iwl_2x2v_p1'      #[ Root name of files to process.

#[ ............... End of user inputs (MAYBE) ..................... ]#

def get_equilibrium_meta(data_dir):
  #[ Return equilibrium metadata for a given shot based on the name.
  out_d = {}
  if '103955_04' in data_dir:
    out_d["R_axis"]      = 0.40045564 #[ R coord of the magnetic axis.
    out_d["Z_axis"]      = 0 #[ Z coord of the magnetic axis.
    out_d["psi_axis"]    = 1.188460310e-03 #[ psi at the magnetic axis.
    out_d["psi_lcfs"]    = -5.632462250e-03 #[ LCFS psi coordinate for 863 mg shot (LTX_103955_03.eqdsk).
    out_d["psi_conv"]    = out_d["psi_axis"] < out_d["psi_lcfs"] #[ =True psi increases outwards, =False it increases inwards.
    out_d["R_LCFS_OMP"]  = 0.605493501 #[ Major radius of the LCFS at OMP.
    out_d["Lc_LCFS_OMP"] = 6.0 #[ Connection length close to the LCFS at the OMP.
    out_d["B0"]          = 0.5*(1.666901e-01+6.814593e-01) #[ Reference magnetic field amplitude.
  elif '103955_03' in data_dir:
    out_d["R_axis"]    = 0.406052 #[ R coord of the magnetic axis.
    out_d["Z_axis"]    = 0 #[ Z coord of the magnetic axis.
    out_d["psi_axis"]  = 1.5428864200000001e-03 #[ psi at the magnetic axis.
    out_d["psi_lcfs"]  = -5.4760172700000003e-03 #[ LCFS psi coordinate for 863 mg shot (LTX_103955_03.eqdsk).
    out_d["psi_conv"]  = out_d["psi_axis"] < out_d["psi_lcfs"] #[ =True psi increases outwards, =False it increases inwards.
  elif '103795_03' in data_dir:
    out_d["R_axis"]      = 0.400392 #[ R coord of the magnetic axis.
    out_d["Z_axis"]      = 0 #[ Z coord of the magnetic axis.
    out_d["psi_axis"]    = 4.3835108700000000e-04 #[ psi at the magnetic axis.
    out_d["psi_lcfs"]    = -5.5560311699999997e-03 #[ LCFS psi coordinate for passivated Li shot (LTX_103795_03.eqdsk).
    out_d["psi_conv"]    = out_d["psi_axis"] < out_d["psi_lcfs"] #[ =True psi increases outwards, =False it increases inwards.
    out_d["R_LCFS_OMP"]  = 0.59476116093 #[ Major radius of the LCFS at OMP.
    out_d["Lc_LCFS_OMP"] = 6.0 #[ Connection length close to the LCFS at the OMP.
    out_d["B0"]          = 0.5*(1.566393e-01+6.265720e-01) #[ Reference magnetic field amplitude.
  else:
    print("get_psi_lcfs_axis: option NYI.")
    print("  data_dir: ",data_dir)
    sys.exit(1)

  return out_d

def get_plasma_meta(data_dir):
  #[ Return plasma metadata/reference parameters.
  out_d = {}
  if '103955_04' in data_dir:
    out_d["Te0"] = 0.5*320*lcu.eV #[ Reference electron temperature.
    out_d["Ti0"] = 0.5*out_d["Te0"] #[ Reference ion temperature.
    out_d["n0"]  = 0.5*1.6e19 #[ Reference density.
  elif '103795_03' in data_dir:
    out_d["Te0"] = 0.5*220*lcu.eV #[ Reference electron temperature.
    out_d["Ti0"] = 0.5*out_d["Te0"] #[ Reference ion temperature.
    out_d["n0"]  = 0.5*2.6e19 #[ Reference density.
  else:
    print("get_psi_lcfs_axis: option NYI.")
    print("  data_dir: ",data_dir)
    sys.exit(1)

  return out_d

file_fmt = '.gkyl' #[ Data file format
poly_order, basis_type = 1, 'ms' #[ Polynomial order and type of basis.

def get_interp_data_c2p(data_file, porder, basis, comp_in, c2p_file):
  #[ Get interpolated data using the c2p coordinate transformation when interpolating.
  pg_data = pg.GData(data_file, mapc2p_name=c2p_file)
  pg_interp = pg.GInterpModal(pg_data, porder, basis)
  x_out, data_out = pg_interp.interpolate(comp_in)
  for i in range(len(x_out)):
    x_out[i] = np.squeeze(x_out[i])
  # end
  data_out = np.squeeze(data_out)
  return x_out, data_out

def getInterpDataComp(file, porder, basis, comp_in, **kwargs):
  #[ Get a specific component from a multicomponent file.
  #[ Inputs
  #[   file: file name.
  #[   porder: polynomial order.
  #[   basis: basis name.
  #[   comp_in: component (string or int).
  #[   kwargs:
  #[     mapc2p: name of mapc2p file.
  #[ Available options for comp_in:
  #[   BiMaxwellianMoments in file: 'den', 'upar', 'tpar', 'tperp', 'temp' or an int.
  #[   MaxwellianMoments in file: 'den', 'upar', 'temp' or an int.
  #[   else: int
  maxwellian_comp_idx = {'den' : 0, 'upar' : 1, 'temp' : 2,}
  bimaxwellian_comp_idx = {'den' : 0, 'upar' : 1, 'tpar' : 2, 'tperp' : 3,}

  if 'mapc2p' not in kwargs:
    if isinstance(comp_in,int):
      return np.squeeze(pgu.getInterpData(file, porder, basis, comp=comp_in))

    elif isinstance(comp_in,str):
      if 'BiMaxwellianMoments' in file:
        if comp_in == 'temp':
          tpar_idx = bimaxwellian_comp_idx['tpar']
          tperp_idx = bimaxwellian_comp_idx['tperp']
          return np.squeeze( (    pgu.getInterpData(file, porder, basis, comp=tpar_idx) + \
                              2.0*pgu.getInterpData(file, porder, basis, comp=tperp_idx) )/3.0 )
        else:
          comp_idx = bimaxwellian_comp_idx[comp_in]
          return np.squeeze(pgu.getInterpData(file, porder, basis, comp=comp_idx))

      elif 'MaxwellianMoments' in file:
        comp_idx = maxwellian_comp_idx[comp_in]
        return np.squeeze(pgu.getInterpData(file, porder, basis, comp=comp_idx))
    
      else:
        print("getInterpDataComp: Component ", comp_in, " is not a valid option")
        sys.exit(1)
  else:
    #[ Transform to R-Z coordinates.
    if isinstance(comp_in,int):
      x_out, data_out = get_interp_data_c2p(file, porder, basis, comp_in, kwargs['mapc2p'])
      return x_out, data_out

    elif isinstance(comp_in,str):
      if 'BiMaxwellianMoments' in file:
        if comp_in == 'temp':
          tpar_idx = bimaxwellian_comp_idx['tpar']
          tperp_idx = bimaxwellian_comp_idx['tperp']
          x_out, tpar_out = get_interp_data_c2p(file, porder, basis, tpar_idx, kwargs['mapc2p'])
          x_out, tperp_out = get_interp_data_c2p(file, porder, basis, tperp_idx, kwargs['mapc2p'])

          return x_out, (tpar_out+2.0*tperp_out)/3.0
        else:
          comp_idx = bimaxwellian_comp_idx[comp_in]
          x_out, data_out = get_interp_data_c2p(file, porder, basis, comp_idx, kwargs['mapc2p'])
          return x_out, data_out

      elif 'MaxwellianMoments' in file:
        comp_idx = maxwellian_comp_idx[comp_in]
        x_out, data_out = get_interp_data_c2p(file, porder, basis, comp_idx, kwargs['mapc2p'])
        return x_out, data_out
    
      else:
        print("getInterpDataComp: Component ", comp_in, " is not a valid option")
        sys.exit(1)

#................................................................................#

if plot_grid_RZ:
  #[ Plot the grid on the R-Z plane.

  data_dir =    '/pscratch/sd/m/mana/gkeyll/ltx/2d/li863mg_103955_04-base/'
  wall_coords = '/global/homes/m/mana/perlmutter/gkeyll/code/gkyl-sims/ltx_gkeyll_xgc/experiment/LTXvessel.csv'  
  psi_pol =     '/pscratch/sd/m/mana/gkeyll/ltx/2d/li863mg_103955_04-base/LTX_103955_04-psi.gkyl'
  psi_contour = True
  psi_contour_num_levels = 30

  fig_file_name_root = lcu.li863_prefix+'grid_RZ'

  #[ Create a click context (needed to call pgkyl command).
  ctx = click.core.Context(cli)
  ctx.obj = {}
  ctx.obj["data"] = pg.commands.DataSpace()
  ctx.obj["verbose"] = False

  ctx.invoke(pg.commands.gk_nodes, name=sim_name, path=data_dir, wall_file=wall_coords,
    psi_file=psi_pol, contour=psi_contour, cnlevels=psi_contour_num_levels, no_show=True)

  #[ Nodes data.
  nodes = ctx.obj["data"].get_dataset(tag="nodes",idx=0)
  nodes_majorR = nodes.get_grid()
  nodes_vertZ = nodes.get_values()

  #[ Edges data.
  edges = ctx.obj["data"].get_dataset(tag="edges",idx=0)
  edges_constx = edges.get_grid()
  edges_consty = edges.get_values()

  plot_h = ctx.obj["plot_handles"]
  #[ Resize figure:
  fig_h = plot_h["figure"]
  fig_width, fig_height = fig_h.get_size_inches() #[ In inches.
  fig_h.set_size_inches(fig_width*0.85, fig_height*0.9)
  #[ Resize axis:
  ax_h = plot_h["axis"]
  ax_h.set_aspect('equal')
  ax_pos = ax_h.get_position().bounds
  ax_h.set_position([ax_pos[0]+0.06, ax_pos[1]-0.02, ax_pos[2]*1.05, ax_pos[3]*1.05])
  #[ Remove psi colorbar.
  if psi_pol:
    plot_h["psi_colorbar"].remove()
    #[ Node data.
    psi = ctx.obj["data"].get_dataset(tag="psi",idx=0)
    psi_coords = psi.get_grid()
    psi_values = psi.get_values()
  # end

  if out_figure_file:
    fig_file_name = output_prefix+fig_file_name_root

    if save_data:
      h5f = h5py.File(out_data_dir+fig_file_name+'.h5', "w")
      h5f.create_dataset('subplot00_nodes_xvalues', np.shape(nodes_majorR), dtype='f8', data=nodes_majorR)
      h5f.create_dataset('subplot00_nodes_yvalues', np.shape(nodes_vertZ ), dtype='f8', data=nodes_vertZ )
      h5f.create_dataset('subplot00_edges_constx', np.shape(edges_constx), dtype='f8', data=edges_constx)
      h5f.create_dataset('subplot00_edges_consty', np.shape(edges_consty), dtype='f8', data=edges_consty)
      if psi_pol:
        h5f.create_dataset('subplot00_psi_xvalues', np.shape(psi_coords[0]), dtype='f8', data=psi_coords[0])
        h5f.create_dataset('subplot00_psi_yvalues', np.shape(psi_coords[1]), dtype='f8', data=psi_coords[1])
        h5f.create_dataset('subplot00_psi_zvalues', np.shape(psi_values   ), dtype='f8', data=psi_values   )
      # end
      h5f.close()

    fig_h.savefig(out_fig_dir+fig_file_name+figure_file_format)
    plt.close()

  else:
    plt.show()


#................................................................................#

if plot_vs_x:
  #[ Plot a variable at the outboard midplane.

  data_dir = '/pscratch/sd/m/mana/gkeyll/ltx/2d/li863mg_103955_04-base/'
  x_axis_psi_N = True #[ Whether to put x-axis in rho_pol.

  quant      = 'elc_BiMaxwellianMoments' #[ Quantity to plot.
  quant_comp = 'temp'                    #[ Component in file (den, upar, tpar, tperp, temp, or an int).
  scale_fac  = lcu.mass_elc/lcu.eV               #[ Factor to multiply data by.
#  ylabel     = r'$n_e(\theta=0,t=0)$ (m$^{-3}$)'       #[ Label for y axis.
#  ylabel     = r'$u_{\parallel e}(\theta=0,t=0)$ (m/s)'       #[ Label for y axis.
  ylabel     = r'$T_e(\theta=0,t=0)$ (eV)'       #[ Label for y axis.
  frame      = 0                         #[ Frame number.

  fig_file_name_root = lcu.li863_prefix+'init_elc_den'

  plotz = 0.0 #[ Computational z coordinate to plot at.

  #[ File with experimental data. Set to None if you don't want to plot exp data.
#  exp_data_file = '../experiment/maan_PoP_2024/Maan_2024-fig2c-li_863mg.csv'
  exp_data_file = '../experiment/maan_PoP_2024/Maan_2024-fig2a-li_863mg.csv'
#  exp_data_file = None
  exp_scale_fac = 1.0

  file_path = data_dir+sim_name+'-'+quant+'_'+str(frame)+file_fmt

  #[ Load the grid.
  xIntC, _, nxIntC, lxIntC, dxIntC, _ = pgu.getGrid(file_path,poly_order,basis_type,location='center')
  
  #[ Get indices along z of slices we wish to plot:
  z_coord = xIntC[1]
  plotzIdx = np.argmin(np.abs(z_coord-plotz))

  #[ Load the data.
#  data = scale_fac*np.squeeze(pgu.getInterpData(file_path, poly_order, basis_type, comp=quant_comp))
  data = scale_fac*getInterpDataComp(file_path, poly_order, basis_type, quant_comp)

  data_slice = data[:,plotzIdx]

  x_coord = xIntC[0]
  xlabel = r'$\psi$ (T m$^2$)'
  if x_axis_psi_N:
    eq_meta = get_equilibrium_meta(data_dir)
    x_coord = lcu.psi_N(x_coord, eq_meta["psi_lcfs"], eq_meta["psi_axis"])
    xlabel = r'$\psi_N$'
    print(f"  psi_N_min = {x_coord[0]:.9e}")
    print(f"  psi_N_max = {x_coord[-1]:.9e}")
    if not eq_meta["psi_conv"]:
      x_coord = x_coord[::-1]
      data_slice = data_slice[::-1]

  #[ Prepare figure.
  fig_prop = (6.4, 4.)
  ax_pos   = [[0.12, 0.15, 0.86, 0.78],]
  fig_h    = plt.figure(figsize=fig_prop)
  ax_h     = [fig_h.add_axes(pos) for pos in ax_pos]

  #[ Plot data
  spl00_line0_x = x_coord
  spl00_line0_y = data_slice

  hpla = list()
  hpla.append(ax_h[0].plot(spl00_line0_x, spl00_line0_y, color=lcu.default_colors[0], linestyle=lcu.default_line_styles[0], marker=lcu.default_markers[0]))

  if exp_data_file is not None:
    #[ Plot experimental data.
    expdata = np.loadtxt(open(exp_data_file),delimiter=',')
    expdata_x, expdata_y = expdata[:,0], expdata[:,1]

    spl00_line1_x = expdata_x 
    spl00_line1_y = expdata_y

    ax_h[0].plot(spl00_line1_x, spl00_line1_y, linestyle=lcu.default_line_styles[1], color='grey')

  ax_h[0].set_xlabel(xlabel, fontsize=lcu.xy_label_font_size, labelpad=0)
  ax_h[0].set_ylabel(ylabel, fontsize=lcu.xy_label_font_size)
  ax_h[0].yaxis.get_offset_text().set_size(lcu.tick_font_size)
  lcu.set_tick_font_size(ax_h[0],lcu.tick_font_size)
  ax_h[0].set_xlim(x_coord[0], x_coord[-1])

  if out_figure_file:
    fig_file_suffix = 'z1slice'
    if abs((plotz-z_coord[0])/z_coord[0]) < 1e-5:
      fig_file_suffix = 'z1min'
    elif abs(plotz-0.5*(z_coord[0]+z_coord[-1])) < 1e-5:
      fig_file_suffix = 'z1mid'
    elif abs((plotz-z_coord[-1])/z_coord[-1]) < 1e-5:
      fig_file_suffix = 'z1max'

    fig_file_name = output_prefix+fig_file_name_root+'_'+fig_file_suffix

    if save_data:
      h5f = h5py.File(out_data_dir+fig_file_name+'.h5', "w")
      h5f.create_dataset('subplot00_line0_xvalues', np.shape(spl00_line0_x), dtype='f8', data=spl00_line0_x)
      h5f.create_dataset('subplot00_line0_yvalues', np.shape(spl00_line0_y), dtype='f8', data=spl00_line0_y)
      if exp_data_file is not None:
        h5f.create_dataset('subplot00_line1_xvalues', np.shape(spl00_line1_x), dtype='f8', data=spl00_line1_x)
        h5f.create_dataset('subplot00_line1_yvalues', np.shape(spl00_line1_y), dtype='f8', data=spl00_line1_y)
      # end
      h5f.close()

    fig_h.savefig(out_fig_dir+fig_file_name+figure_file_format)
    plt.close()

  else:
    plt.show()

#................................................................................#

if plot_nT_vs_x:
  #[ Plot density and temperature vs x:

  data_dir = '/pscratch/sd/m/mana/gkeyll/ltx/2d/li863mg_103955_04-base/'

  x_axis_psi_N = True #[ Whether to put x-axis in rho_pol.
  plot_exp_data = True #[ Whether to plot experimental data.
  
  frame = 14   #[ Frame number.
  y_labels = [
#    r'$n_e(\theta=0,t=0)$ (m$^{-3}$)',
#    r'$T_e(\theta=0,t=0)$ (eV)',
#    r'$T_i(\theta=0,t=0)$ (eV)',
    r'$n_e(\theta=0,t=0.7~\mathrm{ms})$ (m$^{-3}$)',
    r'$T_e(\theta=0,t=0.7~\mathrm{ms})$ (eV)',
    r'$T_i(\theta=0,t=0.7~\mathrm{ms})$ (eV)',
  ]

  fig_file_name_root = lcu.li863_prefix+'final_den_temp'

  plotz = 0.0 #[ Computational z coordinate to plot at.

  file_path = [
    data_dir+sim_name+'-elc_BiMaxwellianMoments_'+str(frame)+file_fmt,
    data_dir+sim_name+'-ion_BiMaxwellianMoments_'+str(frame)+file_fmt,
  ]

  #[ Load the grid.
  xIntC, _, nxIntC, lxIntC, dxIntC, _ = pgu.getGrid(file_path[0],poly_order,basis_type,location='center')
  
  #[ Get indices along z of slices we wish to plot:
  z_coord = xIntC[1]
  plotzIdx = np.argmin(np.abs(z_coord-plotz))

  #[ Load the data.
  elc_dens = getInterpDataComp(file_path[0], poly_order, basis_type, 'den')
  elc_temp = (lcu.mass_elc/lcu.eV)*getInterpDataComp(file_path[0], poly_order, basis_type, 'temp')
  ion_temp = (lcu.mass_ion/lcu.eV)*getInterpDataComp(file_path[1], poly_order, basis_type, 'temp')

  elc_dens_slice = elc_dens[:,plotzIdx]
  elc_temp_slice = elc_temp[:,plotzIdx]
  ion_temp_slice = ion_temp[:,plotzIdx]

  x_coord = xIntC[0]
  xlabel = r'$\psi$ (T m$^2$)'
  if x_axis_psi_N:
    eq_meta = get_equilibrium_meta(data_dir)
    x_coord = lcu.psi_N(x_coord, eq_meta["psi_lcfs"], eq_meta["psi_axis"])
    xlabel = r'$\psi_N$'
    print(f"  psi_N_min = {x_coord[0]:.9e}")
    print(f"  psi_N_max = {x_coord[-1]:.9e}")
    if not eq_meta["psi_conv"]:
      x_coord = x_coord[::-1]
      elc_dens_slice = elc_dens_slice[::-1]
      elc_temp_slice = elc_temp_slice[::-1]
      ion_temp_slice = ion_temp_slice[::-1]

  #[ Prepare figure.
  fig_prop = (12., 3.6)
  ax_pos   = [[0.08, 0.15, 0.25, 0.78],[0.41, 0.15, 0.25, 0.78],[0.74, 0.15, 0.25, 0.78],]
  fig_h    = plt.figure(figsize=fig_prop)
  ax_h     = [fig_h.add_axes(pos) for pos in ax_pos]

  #[ Plot data
  spl00_line0_x = x_coord
  spl01_line0_x = x_coord
  spl02_line0_x = x_coord
  spl00_line0_y = elc_dens_slice
  spl01_line0_y = elc_temp_slice
  spl02_line0_y = ion_temp_slice

  hpla, hplb, hplc = list(), list(), list()
  hpla.append(ax_h[0].plot(spl00_line0_x, spl00_line0_y, color=lcu.default_colors[0], linestyle=lcu.default_line_styles[0], marker=lcu.default_markers[0]))
  hplb.append(ax_h[1].plot(spl01_line0_x, spl01_line0_y, color=lcu.default_colors[0], linestyle=lcu.default_line_styles[0], marker=lcu.default_markers[0]))
  hplc.append(ax_h[2].plot(spl02_line0_x, spl02_line0_y, color=lcu.default_colors[0], linestyle=lcu.default_line_styles[0], marker=lcu.default_markers[0]))

  if plot_exp_data:
    #[ Load and plot experimental data.
    exp_elc_den_file = '../experiment/maan_PoP_2024/Maan_2024-fig2c-li_863mg.csv'
    exp_elc_temp_file = '../experiment/maan_PoP_2024/Maan_2024-fig2a-li_863mg.csv'

    exp_elc_den = np.loadtxt(open(exp_elc_den_file),delimiter=',')
    exp_elc_temp = np.loadtxt(open(exp_elc_temp_file),delimiter=',')

    exp_elc_den_x, exp_elc_den_y = exp_elc_den[:,0], exp_elc_den[:,1]
    exp_elc_temp_x, exp_elc_temp_y = exp_elc_temp[:,0], exp_elc_temp[:,1]

    spl00_line1_x = exp_elc_den_x 
    spl00_line1_y = exp_elc_den_y
    spl01_line1_x = exp_elc_temp_x 
    spl01_line1_y = exp_elc_temp_y

    ax_h[0].plot(spl00_line1_x, spl00_line1_y, linestyle=lcu.default_line_styles[1], color='grey')
    ax_h[1].plot(spl01_line1_x, spl01_line1_y, linestyle=lcu.default_line_styles[1], color='grey')
  #end

  for i in range(len(ax_h)):
    ax_h[i].set_xlabel(xlabel, fontsize=lcu.xy_label_font_size, labelpad=0)
    ax_h[i].set_ylabel(y_labels[i], fontsize=lcu.xy_label_font_size)
    ax_h[i].yaxis.get_offset_text().set_size(lcu.tick_font_size)
    lcu.set_tick_font_size(ax_h[i],lcu.tick_font_size)
    ax_h[i].set_xlim(x_coord[0], x_coord[-1])
  # end

  if out_figure_file:
    fig_file_suffix = 'z1slice'
    if abs((plotz-z_coord[0])/z_coord[0]) < 1e-5:
      fig_file_suffix = 'z1min'
    elif abs(plotz-0.5*(z_coord[0]+z_coord[-1])) < 1e-5:
      fig_file_suffix = 'z1mid'
    elif abs((plotz-z_coord[-1])/z_coord[-1]) < 1e-5:
      fig_file_suffix = 'z1max'

    fig_file_name = output_prefix+fig_file_name_root+'_'+fig_file_suffix

    if save_data:
      h5f = h5py.File(out_data_dir+fig_file_name+'.h5', "w")
      h5f.create_dataset('subplot00_line0_xvalues', np.shape(spl00_line0_x), dtype='f8', data=spl00_line0_x)
      h5f.create_dataset('subplot00_line0_yvalues', np.shape(spl00_line0_y), dtype='f8', data=spl00_line0_y)
      h5f.create_dataset('subplot01_line0_xvalues', np.shape(spl01_line0_x), dtype='f8', data=spl01_line0_x)
      h5f.create_dataset('subplot01_line0_yvalues', np.shape(spl01_line0_y), dtype='f8', data=spl01_line0_y)
      h5f.create_dataset('subplot02_line0_xvalues', np.shape(spl02_line0_x), dtype='f8', data=spl02_line0_x)
      h5f.create_dataset('subplot02_line0_yvalues', np.shape(spl02_line0_y), dtype='f8', data=spl02_line0_y)
      if plot_exp_data:
        h5f.create_dataset('subplot00_line1_xvalues', np.shape(spl00_line1_x), dtype='f8', data=spl00_line1_x)
        h5f.create_dataset('subplot00_line1_yvalues', np.shape(spl00_line1_y), dtype='f8', data=spl00_line1_y)
        h5f.create_dataset('subplot01_line1_xvalues', np.shape(spl01_line1_x), dtype='f8', data=spl01_line1_x)
        h5f.create_dataset('subplot01_line1_yvalues', np.shape(spl01_line1_y), dtype='f8', data=spl01_line1_y)
      # end
      h5f.close()

    fig_h.savefig(out_fig_dir+fig_file_name+figure_file_format)
    plt.close()

  else:
    plt.show()

#................................................................................#

if plot_src_mom_vs_x:
  #[ Plot source moments vs. x:

  data_dir = '/pscratch/sd/m/mana/gkeyll/ltx/2d/li863mg_103955_04-base/'

  x_axis_psi_N = True #[ Whether to put x-axis in rho_pol.
  
  frame = 0   #[ Frame number.
  y_labels = [
    r'$\dot{n}_{\mathrm{src}}$ (m$^{-3}$/s)',
    r'$T_{\mathrm{src},s}$ (eV)',
  ]

  fig_file_name_root = lcu.li863_prefix+'src_den_temp'

  plotz = 0.0 #[ Computational z coordinate to plot at.

  file_path = [
    data_dir+sim_name+'-elc_source_MaxwellianMoments_'+str(frame)+file_fmt,
    data_dir+sim_name+'-ion_source_MaxwellianMoments_'+str(frame)+file_fmt,
  ]

  #[ Load the grid.
  xIntC, _, nxIntC, lxIntC, dxIntC, _ = pgu.getGrid(file_path[0],poly_order,basis_type,location='center')
  
  #[ Get indices along z of slices we wish to plot:
  z_coord = xIntC[1]
  plotzIdx = np.argmin(np.abs(z_coord-plotz))

  #[ Load the data.
  elc_dens = getInterpDataComp(file_path[0], poly_order, basis_type, 'den')
  elc_temp = (lcu.mass_elc/lcu.eV)*getInterpDataComp(file_path[0], poly_order, basis_type, 'temp')
  ion_temp = (lcu.mass_ion/lcu.eV)*getInterpDataComp(file_path[1], poly_order, basis_type, 'temp')

  elc_dens_slice = elc_dens[:,plotzIdx]
  elc_temp_slice = elc_temp[:,plotzIdx]
  ion_temp_slice = ion_temp[:,plotzIdx]

  x_coord = xIntC[0]
  xlabel = r'$\psi$ (T m$^2$)'
  if x_axis_psi_N:
    eq_meta = get_equilibrium_meta(data_dir)
    x_coord = lcu.psi_N(x_coord, eq_meta["psi_lcfs"], eq_meta["psi_axis"])
    xlabel = r'$\psi_N$'
    print(f"  psi_N_min = {x_coord[0]:.9e}")
    print(f"  psi_N_max = {x_coord[-1]:.9e}")
    if not eq_meta["psi_conv"]:
      x_coord = x_coord[::-1]
      elc_dens_slice = elc_dens_slice[::-1]
      elc_temp_slice = elc_temp_slice[::-1]
      ion_temp_slice = ion_temp_slice[::-1]

  #[ Prepare figure.
  fig_prop = (6., 5.)
  ax_pos   = [[0.15, 0.54, 0.83, 0.40],
              [0.15, 0.11, 0.83, 0.40],]
  fig_h    = plt.figure(figsize=fig_prop)
  ax_h     = [fig_h.add_axes(pos) for pos in ax_pos]

  #[ Plot data
  spl00_line0_x = x_coord
  spl01_line0_x = x_coord
  spl01_line1_x = x_coord
  spl00_line0_y = elc_dens_slice
  spl01_line0_y = elc_temp_slice
  spl01_line1_y = ion_temp_slice

  hpla, hplb = list(), list()
  hpla.append(ax_h[0].plot(spl00_line0_x, spl00_line0_y, color=lcu.default_colors[0], linestyle=lcu.default_line_styles[0], marker=lcu.default_markers[0]))
  hplb.append(ax_h[1].plot(spl01_line0_x, spl01_line0_y, color=lcu.default_colors[0], linestyle=lcu.default_line_styles[0], marker=lcu.default_markers[0]))
  hplb.append(ax_h[1].plot(spl01_line1_x, spl01_line1_y, color=lcu.default_colors[1], linestyle=lcu.default_line_styles[1], marker=lcu.default_markers[0]))

  for i in range(len(ax_h)):
    ax_h[i].set_ylabel(y_labels[i], fontsize=lcu.xy_label_font_size)
    ax_h[i].yaxis.get_offset_text().set_size(lcu.tick_font_size)
    lcu.set_tick_font_size(ax_h[i],lcu.tick_font_size)
    ax_h[i].set_xlim(x_coord[0], x_coord[-1])
  # end

  ax_h[1].legend([hplb[0][0],hplb[1][0]],['e$^-$','H$^+$'],fontsize=lcu.legend_font_size, frameon=False, loc='upper right')
  plt.setp( ax_h[0].get_xticklabels(), visible=False)
  ax_h[1].set_xlabel(xlabel, fontsize=lcu.xy_label_font_size, labelpad=0)
  ax_h[0].set_ylim(0.0, 7e22)
  ax_h[1].set_ylim(0.0, 250.0)

  if out_figure_file:
    fig_file_suffix = 'z1slice'
    if abs((plotz-z_coord[0])/z_coord[0]) < 1e-5:
      fig_file_suffix = 'z1min'
    elif abs(plotz-0.5*(z_coord[0]+z_coord[-1])) < 1e-5:
      fig_file_suffix = 'z1mid'
    elif abs((plotz-z_coord[-1])/z_coord[-1]) < 1e-5:
      fig_file_suffix = 'z1max'

    fig_file_name = output_prefix+fig_file_name_root+'_'+fig_file_suffix

    if save_data:
      h5f = h5py.File(out_data_dir+fig_file_name+'.h5', "w")
      h5f.create_dataset('subplot00_line0_xvalues', np.shape(spl00_line0_x), dtype='f8', data=spl00_line0_x)
      h5f.create_dataset('subplot00_line0_yvalues', np.shape(spl00_line0_y), dtype='f8', data=spl00_line0_y)
      h5f.create_dataset('subplot01_line0_xvalues', np.shape(spl01_line0_x), dtype='f8', data=spl01_line0_x)
      h5f.create_dataset('subplot01_line0_yvalues', np.shape(spl01_line0_y), dtype='f8', data=spl01_line0_y)
      h5f.create_dataset('subplot01_line1_xvalues', np.shape(spl01_line1_x), dtype='f8', data=spl01_line1_x)
      h5f.create_dataset('subplot01_line1_yvalues', np.shape(spl01_line1_y), dtype='f8', data=spl01_line1_y)
      h5f.close()

    fig_h.savefig(out_fig_dir+fig_file_name+figure_file_format)
    plt.close()

  else:
    plt.show()

#................................................................................#

if plot_src_int_mom_vs_t:
  #[ Plot source integrated moments vs. t:

  data_dir = '/pscratch/sd/m/mana/gkeyll/ltx/2d/li863mg_103955_04-base/'

  xlabel = r'Time (ms)'
  y_labels = [
    r'$\dot{N}_{\mathrm{src}}$ (s$^{-1}$)',
    r'$P_{\mathrm{src},s}$ (kW)',
  ]

  fig_file_name_root = lcu.li863_prefix+'src_Ndot_Power'

  file_path = [
    data_dir+sim_name+'-elc_source_integrated_moms'+file_fmt,
    data_dir+sim_name+'-ion_source_integrated_moms'+file_fmt,
  ]

  #[ Load the data.
  time, int_mom_elc = pgu.readDynVector(file_path[0])
  time, int_mom_ion = pgu.readDynVector(file_path[1])

  vol_fac = 2.0*np.pi
  Ndot    = vol_fac*int_mom_elc[:,0]
  pow_elc = vol_fac*0.5*lcu.mass_elc*int_mom_elc[:,2]
  pow_ion = vol_fac*0.5*lcu.mass_ion*int_mom_ion[:,2]

  #[ Conver time to ms and power to kW.
  time *= 1e3
  pow_elc *= 1e-3
  pow_ion *= 1e-3

  #[ Prepare figure.
  fig_prop = (6., 5.)
  ax_pos   = [[0.15, 0.54, 0.83, 0.40],
              [0.15, 0.11, 0.83, 0.40],]
  fig_h    = plt.figure(figsize=fig_prop)
  ax_h     = [fig_h.add_axes(pos) for pos in ax_pos]

  #[ Plot data
  spl00_line0_x = time
  spl01_line0_x = time
  spl01_line1_x = time
  spl00_line0_y = Ndot
  spl01_line0_y = pow_elc
  spl01_line1_y = pow_ion

  hpla, hplb = list(), list()
  hpla.append(ax_h[0].plot(spl00_line0_x, spl00_line0_y, color=lcu.default_colors[0], linestyle=lcu.default_line_styles[0], marker=lcu.default_markers[0]))
  hplb.append(ax_h[1].plot(spl01_line0_x, spl01_line0_y, color=lcu.default_colors[0], linestyle=lcu.default_line_styles[0], marker=lcu.default_markers[0]))
  hplb.append(ax_h[1].plot(spl01_line1_x, spl01_line1_y, color=lcu.default_colors[1], linestyle=lcu.default_line_styles[1], marker=lcu.default_markers[0]))

  for i in range(len(ax_h)):
    ax_h[i].set_ylabel(y_labels[i], fontsize=lcu.xy_label_font_size)
    ax_h[i].yaxis.get_offset_text().set_size(lcu.tick_font_size)
    lcu.set_tick_font_size(ax_h[i],lcu.tick_font_size)
    ax_h[i].set_xlim(time[0], time[-1])
  # end

  ax_h[1].legend([hplb[0][0],hplb[1][0]],['e$^-$','H$^+$'],fontsize=lcu.legend_font_size, frameon=False, loc='upper right')
  plt.setp( ax_h[0].get_xticklabels(), visible=False)
  ax_h[1].set_xlabel(xlabel, fontsize=lcu.xy_label_font_size, labelpad=0)
  ax_h[0].set_ylim(3.5e21, 3.8e21)
  ax_h[1].set_ylim(100.0, 250.0)

  if out_figure_file:
    fig_file_name = output_prefix+fig_file_name_root

    if save_data:
      h5f = h5py.File(out_data_dir+fig_file_name+'.h5', "w")
      h5f.create_dataset('subplot00_line0_xvalues', np.shape(spl00_line0_x), dtype='f8', data=spl00_line0_x)
      h5f.create_dataset('subplot00_line0_yvalues', np.shape(spl00_line0_y), dtype='f8', data=spl00_line0_y)
      h5f.create_dataset('subplot01_line0_xvalues', np.shape(spl01_line0_x), dtype='f8', data=spl01_line0_x)
      h5f.create_dataset('subplot01_line0_yvalues', np.shape(spl01_line0_y), dtype='f8', data=spl01_line0_y)
      h5f.create_dataset('subplot01_line1_xvalues', np.shape(spl01_line1_x), dtype='f8', data=spl01_line1_x)
      h5f.create_dataset('subplot01_line1_yvalues', np.shape(spl01_line1_y), dtype='f8', data=spl01_line1_y)
      h5f.close()

    fig_h.savefig(out_fig_dir+fig_file_name+figure_file_format)
    plt.close()

  else:
    plt.show()

#................................................................................#

if plot_vs_RZ:
  #[ Plot a variable on the R-Z midplane.

  data_dir = '/pscratch/sd/m/mana/gkeyll/ltx/2d/li863mg_103955_04-base/'

  quant      = 'elc_BiMaxwellianMoments' #[ Quantity to plot.
  quant_comp = 'den'                    #[ Component in file (den, upar, tpar, tperp, temp, or an int).
  scale_fac  = 1.0               #[ Factor to multiply data by.
#  scale_fac  = lcu.mass_ion/lcu.eV               #[ Factor to multiply data by.
  zlabel     = r'$n_e(t=0.7~\mathrm{ms})$ (m$^{-3}$)'       #[ Label for y axis.
#  zlabel     = r'$u_{\parallel i}(t=0.7~\mathrm{ms})$ (km/s)'       #[ Label for y axis.
#  zlabel     = r'$T_{\parallel e}(t=0.7~\mathrm{ms})$ (eV)'       #[ Label for y axis.
#  zlabel     = r'$T_{\perp e}(t=0.7~\mathrm{ms})$ (eV)'       #[ Label for y axis.
  frame      = 14                         #[ Frame number.

  fig_file_name_root = lcu.li863_prefix+'final_elc_den_RZ'

  wall_file = '/global/homes/m/mana/perlmutter/gkeyll/code/gkyl-sims/ltx_gkeyll_xgc/experiment/LTXvessel.csv'

  xlabel = r'$R$ (m)'
  ylabel = r'$Z$ (m)'

  file_path = data_dir+sim_name+'-'+quant+'_'+str(frame)+file_fmt
  c2p_path  = data_dir+sim_name+'-mapc2p_deflated'+file_fmt

  #[ Load the data.
  xInt, data = getInterpDataComp(file_path, poly_order, basis_type, quant_comp, mapc2p=c2p_path)
  data *= scale_fac

  #[ Prepare figure.
  fig_prop = (5.5, 5.8)
  ax_pos   = [[0.15, 0.1, 0.65, 0.85],]
  cbax_pos = [0.82, 0.1, 0.02, 0.85]
  fig_h    = plt.figure(figsize=fig_prop)
  ax_h     = [fig_h.add_axes(pos) for pos in ax_pos]
  cbax_h   = fig_h.add_axes(cbax_pos)

  ax_h[0].set_aspect('equal')

  #[ Plot data
  spl00_x = xInt[0]
  spl00_y = xInt[1]
  spl00_z = data

  hpla = list()
  hpla.append(ax_h[0].pcolormesh(spl00_x, spl00_y, spl00_z, cmap='inferno'))

  #[ Plot wall
  wall_data = np.loadtxt(open(wall_file),delimiter=',')
  wall_h = ax_h[0].plot(wall_data[:,0],wall_data[:,1],color="grey")

  hcba = plt.colorbar(hpla[0], ax=ax_h[0], cax=cbax_h)
  hcba.ax.tick_params(labelsize=lcu.tick_font_size)
  hcba.set_label(zlabel, rotation=90, labelpad=0, fontsize=lcu.colorbar_label_font_size)
  hcba.ax.yaxis.get_offset_text().set_fontsize(lcu.tick_font_size)

  ax_h[0].set_xlabel(xlabel, fontsize=lcu.xy_label_font_size)
  ax_h[0].xaxis.get_offset_text().set_size(lcu.tick_font_size)
  ax_h[0].set_ylabel(ylabel, fontsize=lcu.xy_label_font_size, labelpad=-3)
  ax_h[0].yaxis.get_offset_text().set_size(lcu.tick_font_size)
  lcu.set_tick_font_size(ax_h[0],lcu.tick_font_size)

  if out_figure_file:
    fig_file_name = output_prefix+fig_file_name_root

    if save_data:
      h5f = h5py.File(out_data_dir+fig_file_name+'.h5', "w")
      h5f.create_dataset('subplot00_xvalues', np.shape(spl00_x), dtype='f8', data=spl00_x)
      h5f.create_dataset('subplot00_yvalues', np.shape(spl00_y), dtype='f8', data=spl00_y)
      h5f.create_dataset('subplot00_zvalues', np.shape(spl00_z), dtype='f8', data=spl00_z)
      h5f.close()

    fig_h.savefig(out_fig_dir+fig_file_name+figure_file_format)
    plt.close()

  else:
    plt.show()

#................................................................................#

if plot_nT_vs_x_multisim:
  #[ Plot density and temperature vs x for multimple sims:

  data_dir = [
    '/pscratch/sd/m/mana/gkeyll/ltx/2d/li863mg_103955_04-base/',
    '/pscratch/sd/m/mana/gkeyll/ltx/2d/lipass_103795_03-base/',
  ]
  legend_strs = [
    lcu.li863_legend,
    lcu.lipass_legend,
  ]

  x_axis_psi_N = True #[ Whether to put x-axis in rho_pol.
  plot_exp_data = False #[ Whether to plot experimental data.
  
  frame = 20   #[ Frame number.
  y_labels = [
#    r'$n_e(\theta=0,t=0)$ (m$^{-3}$)',
#    r'$T_e(\theta=0,t=0)$ (eV)',
#    r'$T_i(\theta=0,t=0)$ (eV)',
    r'$n_e(\theta=0,t=1~\mathrm{ms})$ (m$^{-3}$)',
    r'$T_e(\theta=0,t=1~\mathrm{ms})$ (eV)',
    r'$T_i(\theta=0,t=1~\mathrm{ms})$ (eV)',
  ]

  fig_file_name_root = lcu.liComp863pass_prefix+'final_den_temp'

  plotz = 0.0 #[ Computational z coordinate to plot at.

  #[ Prepare figure.
  fig_prop = (12., 3.6)
  ax_pos   = [[0.08, 0.15, 0.25, 0.78],[0.41, 0.15, 0.25, 0.78],[0.74, 0.15, 0.25, 0.78],]
  fig_h    = plt.figure(figsize=fig_prop)
  ax_h     = [fig_h.add_axes(pos) for pos in ax_pos]
  
  file_path = list()
  spl00_line0_x = list()
  spl01_line0_x = list()
  spl02_line0_x = list()
  spl00_line0_y = list()
  spl01_line0_y = list()
  spl02_line0_y = list()
  hpla, hplb, hplc = list(), list(), list()
  for dI in range(len(data_dir)):

    file_path.append([
        data_dir[dI]+sim_name+'-elc_BiMaxwellianMoments_'+str(frame)+file_fmt,
        data_dir[dI]+sim_name+'-ion_BiMaxwellianMoments_'+str(frame)+file_fmt,
    ])
  
    #[ Load the grid.
    xIntC, _, nxIntC, lxIntC, dxIntC, _ = pgu.getGrid(file_path[dI][0],poly_order,basis_type,location='center')
    
    #[ Get indices along z of slices we wish to plot:
    z_coord = xIntC[1]
    plotzIdx = np.argmin(np.abs(z_coord-plotz))
  
    #[ Load the data.
    elc_dens = getInterpDataComp(file_path[dI][0], poly_order, basis_type, 'den')
    elc_temp = (lcu.mass_elc/lcu.eV)*getInterpDataComp(file_path[dI][0], poly_order, basis_type, 'temp')
    ion_temp = (lcu.mass_ion/lcu.eV)*getInterpDataComp(file_path[dI][1], poly_order, basis_type, 'temp')
  
    elc_dens_slice = elc_dens[:,plotzIdx]
    elc_temp_slice = elc_temp[:,plotzIdx]
    ion_temp_slice = ion_temp[:,plotzIdx]
  
    x_coord = xIntC[0]
    xlabel = r'$\psi$ (T m$^2$)'
    if x_axis_psi_N:
      eq_meta = get_equilibrium_meta(data_dir[dI])
      x_coord = lcu.psi_N(x_coord, eq_meta["psi_lcfs"], eq_meta["psi_axis"])
      xlabel = r'$\psi_N$'
      if not eq_meta["psi_conv"]:
        x_coord = x_coord[::-1]
        elc_dens_slice = elc_dens_slice[::-1]
        elc_temp_slice = elc_temp_slice[::-1]
        ion_temp_slice = ion_temp_slice[::-1]
  
    #[ Plot data
    spl00_line0_x.append(x_coord)
    spl01_line0_x.append(x_coord)
    spl02_line0_x.append(x_coord)
    spl00_line0_y.append(elc_dens_slice)
    spl01_line0_y.append(elc_temp_slice)
    spl02_line0_y.append(ion_temp_slice)
  
    hpla.append(ax_h[0].plot(spl00_line0_x[dI], spl00_line0_y[dI], color=lcu.default_colors[dI], linestyle=lcu.default_line_styles[dI], marker=lcu.default_markers[dI]))
    hplb.append(ax_h[1].plot(spl01_line0_x[dI], spl01_line0_y[dI], color=lcu.default_colors[dI], linestyle=lcu.default_line_styles[dI], marker=lcu.default_markers[dI]))
    hplc.append(ax_h[2].plot(spl02_line0_x[dI], spl02_line0_y[dI], color=lcu.default_colors[dI], linestyle=lcu.default_line_styles[dI], marker=lcu.default_markers[dI]))
  
    if plot_exp_data:
      #[ Load and plot experimental data.
      exp_elc_den_file = '../experiment/maan_PoP_2024/Maan_2024-fig2c-li_863mg.csv'
      exp_elc_temp_file = '../experiment/maan_PoP_2024/Maan_2024-fig2a-li_863mg.csv'
  
      exp_elc_den = np.loadtxt(open(exp_elc_den_file),delimiter=',')
      exp_elc_temp = np.loadtxt(open(exp_elc_temp_file),delimiter=',')
  
      exp_elc_den_x, exp_elc_den_y = exp_elc_den[:,0], exp_elc_den[:,1]
      exp_elc_temp_x, exp_elc_temp_y = exp_elc_temp[:,0], exp_elc_temp[:,1]
  
      spl00_line1_x = exp_elc_den_x 
      spl00_line1_y = exp_elc_den_y
      spl01_line1_x = exp_elc_temp_x 
      spl01_line1_y = exp_elc_temp_y
  
      ax_h[0].plot(spl00_line1_x, spl00_line1_y, linestyle=lcu.default_line_styles[1], color='grey')
      ax_h[1].plot(spl01_line1_x, spl01_line1_y, linestyle=lcu.default_line_styles[1], color='grey')
    #end

  ax_h[0].legend([hpla[0][0],hpla[1][0]],legend_strs,fontsize=lcu.legend_font_size, frameon=False, loc='upper right')
  for i in range(len(ax_h)):
    ax_h[i].set_xlabel(xlabel, fontsize=lcu.xy_label_font_size, labelpad=0)
    ax_h[i].set_ylabel(y_labels[i], fontsize=lcu.xy_label_font_size)
    ax_h[i].yaxis.get_offset_text().set_size(lcu.tick_font_size)
    lcu.set_tick_font_size(ax_h[i],lcu.tick_font_size)
    ax_h[i].set_xlim(x_coord[0], x_coord[-1])
  # end

  if out_figure_file:
    fig_file_suffix = 'z1slice'
    if abs((plotz-z_coord[0])/z_coord[0]) < 1e-5:
      fig_file_suffix = 'z1min'
    elif abs(plotz-0.5*(z_coord[0]+z_coord[-1])) < 1e-5:
      fig_file_suffix = 'z1mid'
    elif abs((plotz-z_coord[-1])/z_coord[-1]) < 1e-5:
      fig_file_suffix = 'z1max'

    fig_file_name = output_prefix+fig_file_name_root+'_'+fig_file_suffix

    if save_data:
      h5f = h5py.File(out_data_dir+fig_file_name+'.h5', "w")
      h5f.create_dataset('subplot00_line0_xvalues', np.shape(spl00_line0_x), dtype='f8', data=spl00_line0_x)
      h5f.create_dataset('subplot00_line0_yvalues', np.shape(spl00_line0_y), dtype='f8', data=spl00_line0_y)
      h5f.create_dataset('subplot01_line0_xvalues', np.shape(spl01_line0_x), dtype='f8', data=spl01_line0_x)
      h5f.create_dataset('subplot01_line0_yvalues', np.shape(spl01_line0_y), dtype='f8', data=spl01_line0_y)
      h5f.create_dataset('subplot02_line0_xvalues', np.shape(spl02_line0_x), dtype='f8', data=spl02_line0_x)
      h5f.create_dataset('subplot02_line0_yvalues', np.shape(spl02_line0_y), dtype='f8', data=spl02_line0_y)
      if plot_exp_data:
        h5f.create_dataset('subplot00_line1_xvalues', np.shape(spl00_line1_x), dtype='f8', data=spl00_line1_x)
        h5f.create_dataset('subplot00_line1_yvalues', np.shape(spl00_line1_y), dtype='f8', data=spl00_line1_y)
        h5f.create_dataset('subplot01_line1_xvalues', np.shape(spl01_line1_x), dtype='f8', data=spl01_line1_x)
        h5f.create_dataset('subplot01_line1_yvalues', np.shape(spl01_line1_y), dtype='f8', data=spl01_line1_y)
      # end
      h5f.close()

    fig_h.savefig(out_fig_dir+fig_file_name+figure_file_format)
    plt.close()

  else:
    plt.show()

#................................................................................#

def calc_plasma_frequency(n, m, eps0, eV):
  #[ Calculate the plasma frequency
  return np.sqrt(n*eV*eV/m/eps0)

def calc_coulomb_log(ns, nr, ms, mr, Ts, Tr, qs, qr, bmag_mid, eps0, hbar, eV):
  #[ Calculate the Coulomb Logarithm
  vts = np.sqrt(Ts/ms) #[ Thermal velocity for species s
  vtr = np.sqrt(Tr/mr)  #[ Thermal velocity for species r
  wps = calc_plasma_frequency(ns,ms, eps0, eV) #[ Plasma Frequency for species s
  wpr = calc_plasma_frequency(nr,mr, eps0, eV) #[ Plasma frequency for species r
  wcs = qs*bmag_mid/ms #[ Cyclotron frequency for species s
  wcr = qr*bmag_mid/mr #[ Cyclotron frequency for species r
  inner1 = (wps*wps + wcs*wcs)/(Ts/ms + 3*Ts/ms) + (wpr*wpr + wcr*wcr)/(Tr/mr + 3*Ts/ms)
  u = 3*(vts*vts + vtr*vtr) #[ Relative velocity
  msr = ms*mr/(ms+mr) #[ Reduced mass
  inner2 = max(np.abs(qs*qr)/(4*np.pi*eps0*msr*u*u), hbar/(2*np.sqrt(eV)*msr*u))
  inner = (1/inner1)*(1/inner2/inner2) + 1
  return 0.5*np.log(inner)

def calc_normNu(ns, nr, ms, mr, qs, qr, Ts, Tr, bmag_mid, eps0, hbar, eV):
  #[ Calculate the normNu
  clog = calc_coulomb_log(ns,nr,ms,mr,Ts, Tr, qs, qr, bmag_mid, eps0, hbar, eV)
  cross_fac2 = 0.0
  if (np.abs(ms - mr)/mr < 1e-16):
    cross_fac2 = 1.0
  else:
    cross_fac2 = 2.0
  return cross_fac2 * (1.0/ms)*(1.0/mr+1.0/ms)*np.power(qs*qr,2)*clog/(3.0*np.power(2.0*np.pi,1.5)*np.power(eps0,2))


if plot_nu_vs_x:
  #[ Plot the collision frequency vs. x.

  plot_nuStar = False
#  data_dir = '/pscratch/sd/m/mana/gkeyll/ltx/2d/li863mg_103955_04-base/'
#  fig_file_name_root = lcu.li863_prefix+'final_'
  data_dir = '/pscratch/sd/m/mana/gkeyll/ltx/2d/lipass_103795_03-base/'
  fig_file_name_root = lcu.lipass_prefix+'final_'
  frame    = 20
  species  = ['elc','elc'] #[ A combination of 'ion' and 'elc', e.g. ['elc','ion'] for nu_ei.
  x_axis_psi_N = True #[ Whether to put x-axis in rho_pol.
  
  #[ Get indices along z of slices we wish to plot:
  plotz = 0.0
  title_str = r'$\theta=0,~t=1~\mathrm{ms}$'

  eq_meta = get_equilibrium_meta(data_dir)
  Lc = eq_meta["Lc_LCFS_OMP"]
  R0 = eq_meta["R_LCFS_OMP"]
  r0 = eq_meta["R_LCFS_OMP"] - eq_meta["R_axis"]
  B0 = eq_meta["B0"]
#  nuStar_Fac = q0*R0/np.power(r0/R0, 3./2.)
  nuStar_Fac = Lc/(2.0*np.pi*np.power(r0/R0, 3./2.))

  plasma_meta = get_plasma_meta(data_dir)
  n0 = plasma_meta["n0"]
  Te0 = plasma_meta["Te0"]
  Ti0 = plasma_meta["Ti0"]

  me = lcu.mass_elc
  mi = lcu.mass_ion
  qe = lcu.charge_elc
  qi = lcu.charge_ion

  ylabel = ''
  nu_file_str = ''
  normNu = -1.0
  if species[0] == 'elc':
    if species[1] == 'elc':
      ylabel = r'$\nu_{ee}$'
      nu_file_str = 'nu_ee'
      normNu = calc_normNu(n0, n0, me, me, qe, qe, Te0, Te0, B0, lcu.eps0, lcu.hbar, lcu.eV)
    elif species[1] == 'ion':                                                            
      ylabel = r'$\nu_{ei}$'                                                       
      nu_file_str = 'nu_ei'                                                           
      normNu = calc_normNu(n0, n0, me, mi, qe, qi, Te0, Ti0, B0, lcu.eps0, lcu.hbar, lcu.eV)
  elif species[0] == 'ion':                                                              
    if species[1] == 'elc':                                                              
      ylabel = r'$\nu_{ie}$'                                                       
      nu_file_str = 'nu_ie'                                                           
      normNu = calc_normNu(n0, n0, mi, me, qi, qe, Ti0, Te0, B0, lcu.eps0, lcu.hbar, lcu.eV)
    elif species[1] == 'ion':                                                            
      ylabel = r'$\nu_{ii}$'                                                       
      nu_file_str = 'nu_ii'                                                           
      normNu = calc_normNu(n0, n0, mi, mi, qi, qi, Ti0, Ti0, B0, lcu.eps0, lcu.hbar, lcu.eV)

  if plot_nuStar:
    assert species[0] == species[1], 'Species must be the same is plot_nuStar == True'
    if species[0] == 'elc':
      ylabel = r'$\nu_{e}^*=\nu_{ee}L_c/(2\pi\epsilon^{3/2}v_{te})$'
      nu_file_str = 'nu_eStar'
    else:
      ylabel = r'$\nu_{i}^*=\nu_{ii}L_c/(2\pi\epsilon^{3/2}v_{ti})$'
      nu_file_str = 'nu_iStar'

  file_path = data_dir+sim_name+'-%s_BiMaxwellianMoments_'+str(frame)+file_fmt

  #[ Load the grid.
  xIntC, _, nxIntC, lxIntC, dxIntC, _ = pgu.getGrid(file_path % 'elc',poly_order,basis_type,location='center')
  
  x_coord = xIntC[0]

  #[ Get indices along z of slices we wish to plot:
  z_coord = xIntC[1]
  plotzIdx = np.argmin(np.abs(z_coord-plotz))

  #[ Load the data.
  dens_r = getInterpDataComp(file_path % species[1], poly_order, basis_type, 'den')
  vtsq_r = getInterpDataComp(file_path % species[1], poly_order, basis_type, 'temp')
  vtsq_s = getInterpDataComp(file_path % species[0], poly_order, basis_type, 'temp')

  dens_r_slice = dens_r[:,plotzIdx]
  vtsq_r_slice = vtsq_r[:,plotzIdx]
  vtsq_s_slice = vtsq_s[:,plotzIdx]

  #[ Collision frequency.
  nu_sr = normNu * dens_r_slice / np.power(np.sqrt(vtsq_s_slice+vtsq_r_slice),3)
  if plot_nuStar:
    nu_sr = nu_sr * nuStar_Fac / np.sqrt(vtsq_s_slice)

  xlabel = r'$\psi$ (T m$^2$)'
  if x_axis_psi_N:
    eq_meta = get_equilibrium_meta(data_dir)
    x_coord = lcu.psi_N(x_coord, eq_meta["psi_lcfs"], eq_meta["psi_axis"])
    xlabel = r'$\psi_N$'
    if not eq_meta["psi_conv"]:
      x_coord = x_coord[::-1]
      data_slice = nu_sr[::-1]

  #[ Prepare figure.
  fig_prop = (6.5, 3.5)
  ax_pos   = [[0.18, 0.16, 0.8, 0.76],]
  fig_h    = plt.figure(figsize=fig_prop)
  ax_h     = [fig_h.add_axes(pos) for pos in ax_pos]

  #[ Plot data
  spl00_line0_x = x_coord
  spl00_line0_y = data_slice

  hpla = list()
  hpla.append(ax_h[0].plot(spl00_line0_x, spl00_line0_y, color=lcu.default_colors[0],
                           linestyle=lcu.default_line_styles[0], marker=lcu.default_markers[0]))

  ax_h[0].set_xlabel(xlabel, fontsize=lcu.xy_label_font_size, labelpad=0)
  ax_h[0].set_ylabel(ylabel, fontsize=lcu.xy_label_font_size, labelpad=0)
  ax_h[0].set_title(title_str, fontsize=lcu.xy_label_font_size)
  lcu.set_tick_font_size(ax_h[0],lcu.tick_font_size)
  ax_h[0].set_xlim(x_coord[0], x_coord[-1])

  if out_figure_file:
    fig_file_suffix = 'z1slice'
    if abs((plotz-z_coord[0])/z_coord[0]) < 1e-5:
      fig_file_suffix = 'z1min'
    elif abs(plotz-0.5*(z_coord[0]+z_coord[-1])) < 1e-5:
      fig_file_suffix = 'z1mid'
    elif abs((plotz-z_coord[-1])/z_coord[-1]) < 1e-5:
      fig_file_suffix = 'z1max'

    fig_file_name = output_prefix+fig_file_name_root+nu_file_str+'_'+fig_file_suffix

    if save_data:
      h5f = h5py.File(out_data_dir+fig_file_name+'.h5', "w")
      h5f.create_dataset('subplot00_nodes_xvalues', np.shape(spl00_line0_x), dtype='f8', data=spl00_line0_x)
      h5f.create_dataset('subplot00_nodes_yvalues', np.shape(spl00_line0_y), dtype='f8', data=spl00_line0_y)
      h5f.close()

    fig_h.savefig(out_fig_dir+fig_file_name+figure_file_format)
    plt.close()

  else:
    plt.show()


#................................................................................#

if plot_int_mom_vs_t:
  #[ Plot integrated moments vs. t:

#  data_dir = '/pscratch/sd/m/mana/gkeyll/ltx/2d/li863mg_103955_04-base/'
#  fig_file_name_root = lcu.li863_prefix+'int_moms'
  data_dir = '/pscratch/sd/m/mana/gkeyll/ltx/2d/lipass_103795_03-base/'
  fig_file_name_root = lcu.lipass_prefix+'int_moms'

  xlabel = r'Time (ms)'
  y_labels = [
    r'$\left\langle n_s\right\rangle$',
    r'$\left\langle n_s u_{\parallel s}\right\rangle$',
    r'$\left\langle n_s\left(u_{\parallel s}^2+3v_{ts}^2\right)\right\rangle$',
  ]

  file_path = [
    data_dir+sim_name+'-elc_integrated_moms'+file_fmt,
    data_dir+sim_name+'-ion_integrated_moms'+file_fmt,
  ]

  #[ Load the data.
  time, int_mom_elc = pgu.readDynVector(file_path[0])
  time, int_mom_ion = pgu.readDynVector(file_path[1])

  m0_elc = int_mom_elc[:,0]
  m1_elc = int_mom_elc[:,1]
  m2_elc = int_mom_elc[:,2]+int_mom_elc[:,3]
  m0_ion = int_mom_ion[:,0]
  m1_ion = int_mom_ion[:,1]
  m2_ion = int_mom_ion[:,2]+int_mom_ion[:,3]

  #[ Conver time to ms.
  time *= 1e3

  #[ Normalize to the final value.
  m0_elc *= 1.0/np.abs(m0_elc[-1])
  m1_elc *= 1.0/np.abs(m1_elc[-1])
  m2_elc *= 1.0/np.abs(m2_elc[-1])
  m0_ion *= 1.0/np.abs(m0_ion[-1])
  m1_ion *= 1.0/np.abs(m1_ion[-1])
  m2_ion *= 1.0/np.abs(m2_ion[-1])

  #[ Prepare figure.
  fig_prop = (6., 6.5)
  ax_pos   = [[0.15, 0.69, 0.83, 0.27],
              [0.15, 0.40, 0.83, 0.27],
              [0.15, 0.11, 0.83, 0.27],]
  fig_h    = plt.figure(figsize=fig_prop)
  ax_h     = [fig_h.add_axes(pos) for pos in ax_pos]

  #[ Plot data
  spl00_line0_x = time
  spl10_line0_x = time
  spl20_line0_x = time
  spl00_line0_y = m0_elc
  spl10_line0_y = m1_elc
  spl20_line0_y = m2_elc

  spl00_line1_x = time
  spl10_line1_x = time
  spl20_line1_x = time
  spl00_line1_y = m0_ion
  spl10_line1_y = m1_ion
  spl20_line1_y = m2_ion

  hpla, hplb, hplc = list(), list(), list()
  hpla.append(ax_h[0].plot(spl00_line0_x, spl00_line0_y, color=lcu.default_colors[0], linestyle=lcu.default_line_styles[0], marker=lcu.default_markers[0]))
  hpla.append(ax_h[0].plot(spl00_line1_x, spl00_line1_y, color=lcu.default_colors[1], linestyle=lcu.default_line_styles[1], marker=lcu.default_markers[1]))
  hpla.append(ax_h[1].plot(spl10_line0_x, spl10_line0_y, color=lcu.default_colors[0], linestyle=lcu.default_line_styles[0], marker=lcu.default_markers[0]))
  hpla.append(ax_h[1].plot(spl10_line1_x, spl10_line1_y, color=lcu.default_colors[1], linestyle=lcu.default_line_styles[1], marker=lcu.default_markers[1]))
  hpla.append(ax_h[2].plot(spl20_line0_x, spl20_line0_y, color=lcu.default_colors[0], linestyle=lcu.default_line_styles[0], marker=lcu.default_markers[0]))
  hpla.append(ax_h[2].plot(spl20_line1_x, spl20_line1_y, color=lcu.default_colors[1], linestyle=lcu.default_line_styles[1], marker=lcu.default_markers[1]))

  for i in range(len(ax_h)):
    ax_h[i].set_ylabel(y_labels[i], fontsize=lcu.xy_label_font_size)
    ax_h[i].yaxis.get_offset_text().set_size(lcu.tick_font_size)
    lcu.set_tick_font_size(ax_h[i],lcu.tick_font_size)
    ax_h[i].set_xlim(time[0], 1.0)
#    ax_h[i].set_xlim(time[0], time[-1])
  # end

  ax_h[0].legend([hpla[0][0],hpla[1][0]],['e$^-$','H$^+$'],fontsize=lcu.legend_font_size, frameon=False, loc='upper right')
  ax_h[-1].set_xlabel(xlabel, fontsize=lcu.xy_label_font_size, labelpad=0)
  for i in range(2):
    plt.setp( ax_h[i].get_xticklabels(), visible=False)

  if out_figure_file:
    fig_file_name = output_prefix+fig_file_name_root

    if save_data:
      h5f = h5py.File(out_data_dir+fig_file_name+'.h5', "w")
      h5f.create_dataset('subplot00_line0_xvalues', np.shape(spl00_line0_x), dtype='f8', data=spl00_line0_x)
      h5f.create_dataset('subplot00_line0_yvalues', np.shape(spl00_line0_y), dtype='f8', data=spl00_line0_y)
      h5f.create_dataset('subplot00_line1_xvalues', np.shape(spl00_line1_x), dtype='f8', data=spl00_line1_x)
      h5f.create_dataset('subplot00_line1_yvalues', np.shape(spl00_line1_y), dtype='f8', data=spl00_line1_y)
      h5f.create_dataset('subplot10_line0_xvalues', np.shape(spl10_line0_x), dtype='f8', data=spl10_line0_x)
      h5f.create_dataset('subplot10_line0_yvalues', np.shape(spl10_line0_y), dtype='f8', data=spl10_line0_y)
      h5f.create_dataset('subplot10_line1_xvalues', np.shape(spl10_line1_x), dtype='f8', data=spl10_line1_x)
      h5f.create_dataset('subplot10_line1_yvalues', np.shape(spl10_line1_y), dtype='f8', data=spl10_line1_y)
      h5f.create_dataset('subplot20_line0_xvalues', np.shape(spl20_line0_x), dtype='f8', data=spl20_line0_x)
      h5f.create_dataset('subplot20_line0_yvalues', np.shape(spl20_line0_y), dtype='f8', data=spl20_line0_y)
      h5f.create_dataset('subplot20_line1_xvalues', np.shape(spl20_line1_x), dtype='f8', data=spl20_line1_x)
      h5f.create_dataset('subplot20_line1_yvalues', np.shape(spl20_line1_y), dtype='f8', data=spl20_line1_y)
      h5f.close()

    fig_h.savefig(out_fig_dir+fig_file_name+figure_file_format)
    plt.close()

  else:
    plt.show()
#................................................................................#
