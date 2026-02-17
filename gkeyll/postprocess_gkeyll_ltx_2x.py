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
import matplotlib.pyplot as plt
import sys #[ For error exit.
import h5py #[ For saving reduced data.
import pgkylUtil as pgu #[ Some postgkyl wrappers.

#[ Plotting options.
plot_vs_x = True  #[ Plot a quantity at the outboard midplane.

out_data_dir  = './data/'
out_fig_dir   = './figures/'
output_prefix = 'ltxb_gkeyll_'

save_data          = True    #[ Indicate whether to save data in plot to HDF5 file.
out_figure_file    = True     #[ Output a figure file?.
figure_file_format = '.png'    #[ Can be .png, .pdf, .ps, .eps, .svg.

sim_name   = 'gk_ltx_iwl_2x2v_p1'      #[ Root name of files to process.

#[ ............... End of user inputs (MAYBE) ..................... ]#

eV        = 1.602176487e-19
me, mp    = 9.10938215e-31, 1.672621637e-27
charge_elc, charge_ion = -eV, eV
mass_elc, mass_ion = me, mp

def get_equilibrium_meta(data_dir):
  #[ Return the axis and LCFS psi for a given shot based on the name.
  out_d = {}
  if '103955_03' in data_dir:
    out_d["R_axis"]    = 0.406052 #[ R coord of the magnetic axis.
    out_d["Z_axis"]    = 0 #[ Z coord of the magnetic axis.
    out_d["psi_axis"]  = 1.5428864200000001e-03 #[ psi at the magnetic axis.
    out_d["psi_lcfs"]  = -5.4760172700000003e-03 #[ LCFS psi coordinate for 863 mg shot (LTX_103955_03.eqdsk).
    out_d["psi_conv"]  = out_d["psi_axis"] < out_d["psi_lcfs"] #[ =True psi increases outwards, =False it increases inwards.
  elif '103795_03' in data_dir:
    out_d["R_axis"]    = 0.400392 #[ R coord of the magnetic axis.
    out_d["Z_axis"]    = 0 #[ Z coord of the magnetic axis.
    out_d["psi_axis"]  = 4.3835108700000000e-04 #[ psi at the magnetic axis.
    out_d["psi_lcfs"]  = -5.5560311699999997e-03 #[ LCFS psi coordinate for 863 mg shot (LTX_103955_03.eqdsk).
    out_d["psi_conv"]  = out_d["psi_axis"] < out_d["psi_lcfs"] #[ =True psi increases outwards, =False it increases inwards.
  else:
    print("get_psi_lcfs_axis: option NYI.")
    sys.exit(1)

  return out_d

file_fmt = '.gkyl' #[ Data file format
poly_order, basis_type = 1, 'ms' #[ Polynomial order and type of basis.

#[ Some RGB colors. These are MATLAB-like.
default_blue     = [0, 0.4470, 0.7410]
default_orange   = [0.8500, 0.3250, 0.0980]
default_green    = [0.4660, 0.6740, 0.1880]
default_purple   = [0.4940, 0.1840, 0.5560]
default_red      = [0.6350, 0.0780, 0.1840]
default_sky_blue = [0.3010, 0.7450, 0.9330]
grey             = [0.5, 0.5, 0.5]
#[ Colors in a single array.
default_colors = [default_blue,default_orange,default_green,default_purple,default_red,default_sky_blue,grey,'black']

#[ LineStyles in a single array.
default_line_styles = ['-','--',':','-.','None','None','None','None']
default_markers     = ['None','None','None','None','o','d','s','+']

#[ Some fontsizes used in plots.
xy_label_font_size       = 17
title_font_size          = 17
colorbar_label_font_size = 17
tick_font_size           = 14
legend_font_size         = 14
text_font_size           = 16

#[ Set the font size of the ticks to a given size.
def set_tick_font_size(axIn,fontSizeIn):
  axIn.tick_params(axis='both', which='major', labelsize=fontSizeIn)

#[ Normalized poloidal flux coordinate rho.
def rho_psi(psi, psisep, psi_ax):
  return np.sqrt((psi - psi_ax) / (psisep - psi_ax))

#[ Flux in terms of normalized poloidal flux rho.
def psi_rho(rho, psisep, psi_ax):
  return np.power(rho,2)*(psisep - psi_ax) + psi_ax

#[ Normalized psi.
def psi_N(psi, psi_lcfs, psi_axis):
  return (psi - psi_axis) / (psi_lcfs - psi_axis)

#def getInterpDataComp(file, porder, basis, comp_in):
#  #[ Get a specific component from a multicomponent file.
#  #[ Inputs
#  #[   file: file name.
#  #[   porder: polynomial order.
#  #[   basis: basis name.
#  #[   comp_in: component (string or int).
#  #[ Available options for comp_in:
#  #[   BiMaxwellianMoments in file: 'den', 'upar', 'tpar', 'tperp', 'temp' or an int.
#  #[   MaxwellianMoments in file: 'den', 'upar', 'temp' or an int.
#  #[   else: int
#  maxwellian_comp_idx = {'den' : 0, 'upar' : 1, 'temp' : 2,}
#  bimaxwellian_comp_idx = {'den' : 0, 'upar' : 1, 'tpar' : 2, 'tperp' : 2,}
#
#  if isinstance(comp_in,int):
#    return np.squeeze(pgu.getInterpData(file, porder, basis, comp=comp_in))
#
#  elif isinstance(comp_in,str):
#    if 'MaxwellianMoments' in file:
#      comp_idx = maxwellian_comp_idx[comp_in]
#      return np.squeeze(pgu.getInterpData(file, porder, basis, comp=comp_idx))
#  
#    elif 'BiMaxwellianMoments' in file:
#      if comp_in == 'temp':
#
#      comp_idx = bimaxwellian_comp_idx[comp_in]
#      return np.squeeze(pgu.getInterpData(file, porder, basis, comp=comp_idx))
#  
#    else:

#................................................................................#

if plot_vs_x:
  #[ Plot a variable at the outboard midplane.

  data_dir = '/Users/mfrancis/Documents/gkeyll/code/gkyl-sims/ltx/sim_data/numerical_eq/2x/li_863mg_103955_03/gn0/'
  x_axis_psi_N = True #[ Whether to put x-axis in rho_pol.

  quant      = 'elc_BiMaxwellianMoments' #[ Quantity to plot.
  quant_comp = 0                    #[ Component in file (den, upar, tpar, tperp, or an int).
  scale_fac  = 1.0               #[ Factor to multiply data by.
  ylabel     = r'$n_e(\theta=0,t=0)$ (m$^{-3}$)'       #[ Label for y axis.
  frame      = 0                         #[ Frame number.

  fig_file_name_root = 'elc_den'

  plotz = 0.0 #[ Computational z coordinate to plot at.

  #[ File with experimental data. Set to None if you don't want to plot exp data.
  exp_data_file = '/Users/mfrancis/Documents/gkeyll/code/gkyl-sims/ltx/maan_2024_data/Maan_2024-fig2c-li_863mg.csv'
  exp_scale_fac = 1.0

  file_path = data_dir+sim_name+'-'+quant+'_'+str(frame)+file_fmt

  #[ Load the grid.
  xIntC, _, nxIntC, lxIntC, dxIntC, _ = pgu.getGrid(file_path,poly_order,basis_type,location='center')
  
  #[ Get indices along z of slices we wish to plot:
  z_coord = xIntC[1]
  plotzIdx = np.argmin(np.abs(z_coord-plotz))

  #[ Load the data.
  data = scale_fac*np.squeeze(pgu.getInterpData(file_path, poly_order, basis_type, comp=quant_comp))

  data_slice = data[:,plotzIdx]

  x_coord = xIntC[0]
  xlabel = r'$\psi$ (T m$^2$)'
  if x_axis_psi_N:
    eq_meta = get_equilibrium_meta(data_dir)
    x_coord = psi_N(x_coord, eq_meta["psi_lcfs"], eq_meta["psi_axis"])
    xlabel = r'$\psi_N$'
    print(f"  rho_min = {x_coord[0]:.9e}")
    print(f"  rho_max = {x_coord[-1]:.9e}")
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
  hpla.append(ax_h[0].plot(spl00_line0_x, spl00_line0_y, color=default_colors[0], linestyle=default_line_styles[0], marker=default_markers[0]))

  if exp_data_file is not None:
    #[ Plot experimental data.
    expdata = np.loadtxt(open(exp_data_file),delimiter=',')
    expdata_x, expdata_y = expdata[:,0], expdata[:,1]

    spl00_line1_x = expdata_x 
    spl00_line1_y = expdata_y

    ax_h[0].plot(spl00_line1_x, spl00_line1_y, linestyle=default_line_styles[1], color='grey')

  ax_h[0].set_xlabel(xlabel, fontsize=xy_label_font_size, labelpad=0)
  ax_h[0].set_ylabel(ylabel, fontsize=xy_label_font_size)
  ax_h[0].yaxis.get_offset_text().set_size(tick_font_size)
  set_tick_font_size(ax_h[0],tick_font_size)
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
      h5f.close()

    fig_h.savefig(out_fig_dir+fig_file_name+figure_file_format)
    plt.close()

  else:
    plt.show()
