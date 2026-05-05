#[ ........................................................... ]#
#[
#[ Comparison of 2D Gkeyll and XGC LTX-beta simulations.
#[
#[ This script plots reduced data and should be able to run
#[ without specific dependencies needed to process raw Gkeyll
#[ or XGC data, using only:
#[   - numpy
#[   - matplotlib
#[   - h5py
#[
#[ Manaure Francisquez & George Wilkie.
#[ February 2026.
#[
#[ ........................................................... ]#

import numpy as np
import matplotlib.pyplot as plt
import h5py
from matplotlib.collections import LineCollection #[ For plotting edges.
import matplotlib.tri as mtri
import sys #[ For error exit.
#[ Append path to utilities folder.
sys.path.insert(0, '../util/')
import ltx_common_util as lcu

plot_grids_RZ_li863mg           = False  #[ Grids on the R-Z plane of 863 mg simulation.
plot_den_omp_init_li863mg       = False  #[ Initial n at outboard midplane (OMP) of 863 mg simulation.
plot_den_temp_omp_init_li863mg  = False  #[ Initial n and T profiles at OMP.
plot_den_temp_omp_final_li863mg = True  #[ Final n and T profiles at OMP.
plot_moms_RZ_final_li863mg      = False  #[ Final n, upar, Tpar and Tperp profiles on R-Z plane.

gke_data_dir = '../gkeyll/data/' #[ Location of reduced Gkeyll data.
xgc_data_dir = '../xgc/data/' #[ Location of reduced XGC data.

out_data_dir  = './data/'
out_fig_dir   = './figures/'
output_prefix = 'ltx_sims_'

output_figure_file = True     #[ Output a figure file?.
figure_file_format = '.png'    #[ Can be .png, .pdf, .ps, .eps, .svg.
save_data          = True    #[ Indicate whether to save data in plot to HDF5 file.

#[ Vacuum vessel wall coordinates.
ltx_vv_file = '../experiment/LTXvessel.csv'  

#[ ............... End of user inputs (MAYBE) ..................... ]#

#[ Default colors for Gkeyll and XGC.
gke_color = lcu.default_colors[0]
xgc_color = lcu.default_colors[1]
#[ Default line styles for Gkeyll and XGC.
gke_linestyle = lcu.default_line_styles[0]
xgc_linestyle = lcu.default_line_styles[1]

#................................................................................#

if plot_grids_RZ_li863mg:
  #[ Plot R-Z grids for the 863 mg Li simulations.
  fig_name = lcu.li863_prefix+'grids_RZ'

  nodes_color = "black" #[ Color for cell nodes.
  edges_color = lcu.default_blue #[ Color for cell edges.
  wall_color = "grey" #[ Color for the wall

  gke_data_file = gke_data_dir+'ltx_gkeyll_li863mg_grid_RZ.h5'
  xgc_data_file = xgc_data_dir+'ltx_xgc_li863mg_grid_RZ.h5'

  #[ Load Gkeyll nodes, edges and psi.
  gke_data = h5py.File(gke_data_file, "r")
  gke_spl00_nodes_x      = lcu.h5data_to_numpy_array(gke_data, 'subplot00_nodes_xvalues')
  gke_spl00_nodes_y      = lcu.h5data_to_numpy_array(gke_data, 'subplot00_nodes_yvalues')
  gke_spl00_edges_constx = lcu.h5data_to_numpy_array(gke_data, 'subplot00_edges_constx')
  gke_spl00_edges_consty = lcu.h5data_to_numpy_array(gke_data, 'subplot00_edges_consty')
  gke_spl00_psi_x        = lcu.h5data_to_numpy_array(gke_data, 'subplot00_psi_xvalues')
  gke_spl00_psi_y        = lcu.h5data_to_numpy_array(gke_data, 'subplot00_psi_yvalues')
  gke_spl00_psi_z        = lcu.h5data_to_numpy_array(gke_data, 'subplot00_psi_zvalues')
  gke_data.close()

  #[ Load XGC nodes, edges and psi.
  xgc_data = h5py.File(xgc_data_file, "r")
  xgc_spl01_nodes_x      = lcu.h5data_to_numpy_array(xgc_data, 'subplot00_nodes_xvalues')
  xgc_spl01_nodes_y      = lcu.h5data_to_numpy_array(xgc_data, 'subplot00_nodes_yvalues')
  xgc_spl01_conn         = lcu.h5data_to_numpy_array(xgc_data, 'subplot00_line_connections')
#  xgc_spl01_edges_constx = lcu.h5data_to_numpy_array(xgc_data, 'subplot00_edges_constx')
#  xgc_spl01_edges_consty = lcu.h5data_to_numpy_array(xgc_data, 'subplot00_edges_consty')
#  xgc_spl01_psi_x        = lcu.h5data_to_numpy_array(xgc_data, 'subplot00_psi_xvalues')
#  xgc_spl01_psi_y        = lcu.h5data_to_numpy_array(xgc_data, 'subplot00_psi_yvalues')
#  xgc_spl01_psi_z        = lcu.h5data_to_numpy_array(xgc_data, 'subplot00_psi_zvalues')
  xgc_data.close()

  triang = mtri.Triangulation(xgc_spl01_nodes_x,xgc_spl01_nodes_y,xgc_spl01_conn)

  #[ Load wall coordinates.
  wall_data = np.loadtxt(open(ltx_vv_file),delimiter=',')
  wall_x = wall_data[:,0]
  wall_y = wall_data[:,1]

  #[ Prepare figure.
  fig_prop = (10.0, 8.)
  ax_pos   = [[0.1, 0.1, 0.4, 0.82],[0.58, 0.1, 0.4, 0.82]]
  fig_h    = plt.figure(figsize=fig_prop)
  ax_h     = [fig_h.add_axes(pos) for pos in ax_pos]

  #[ Plot Gkeyll nodes, edges and psi.
  hpla = list()
  hpla.append(ax_h[0].plot(gke_spl00_nodes_x,gke_spl00_nodes_y,marker=".", color=nodes_color, linestyle="none"))
  hpla.append(ax_h[0].add_collection(LineCollection(gke_spl00_edges_constx, color=edges_color)))
  hpla.append(ax_h[0].add_collection(LineCollection(gke_spl00_edges_consty, color=edges_color)))
  hpla.append(ax_h[0].contour(gke_spl00_psi_x, gke_spl00_psi_y, gke_spl00_psi_z))
  plt.text(0.06, 0.93, r'(a) Gkeyll', fontsize=16, color='black', transform=ax_h[0].transAxes)

  #[ Plot XGC nodes, edges and psi.
  hplb = list()
#  hplb.append(ax_h[1].plot(xgc_spl01_nodes_x,xgc_spl01_nodes_y,marker=".", color=nodes_color, linestyle="none"))
  hplb.append(ax_h[1].triplot(triang, color=edges_color,linewidth=0.2))
#  hplb.append(ax_h[1].add_collection(LineCollection(xgc_spl01_edges_consty, color=edges_color)))
#  hplb.append(ax_h[1].add_collection(LineCollection(xgc_spl01_edges_constx, color=edges_color)))
#  hplb.append(ax_h[1].add_collection(LineCollection(xgc_spl01_edges_consty, color=edges_color)))
#  hplb.append(ax_h[1].contour(xgc_spl01_psi_x, xgc_spl01_psi_y, xgc_spl01_psi_z))
  plt.text(0.06, 0.93, r'(b) XGC', fontsize=16, color='black', transform=ax_h[1].transAxes)

  ax_h[0].set_ylabel(r'Z (m)', fontsize=lcu.xy_label_font_size)
  wall_h = list()
  for i in range(len(ax_h)):
    ax_h[i].set_aspect('equal')
    ax_h[i].set_xlabel(r'R (m)', fontsize=lcu.xy_label_font_size, labelpad=-2)
    lcu.set_tick_font_size(ax_h[i],lcu.tick_font_size)
    wall_h.append(ax_h[i].plot(wall_x,wall_y,color=wall_color))

  if output_figure_file:
    fig_file_name = output_prefix+fig_name

    if save_data:
      h5f = h5py.File(out_data_dir+fig_file_name+'.h5', "w")
      h5f.create_dataset('subplot00_gke_nodes_xvalues', np.shape(gke_spl00_nodes_x     ), dtype='f8', data=gke_spl00_nodes_x     )
      h5f.create_dataset('subplot00_gke_nodes_yvalues', np.shape(gke_spl00_nodes_y     ), dtype='f8', data=gke_spl00_nodes_y     )
      h5f.create_dataset('subplot00_gke_edges_constx' , np.shape(gke_spl00_edges_constx), dtype='f8', data=gke_spl00_edges_constx)
      h5f.create_dataset('subplot00_gke_edges_consty' , np.shape(gke_spl00_edges_consty), dtype='f8', data=gke_spl00_edges_consty)
      h5f.create_dataset('subplot00_gke_psi_xvalues'  , np.shape(gke_spl00_psi_x       ), dtype='f8', data=gke_spl00_psi_x       )
      h5f.create_dataset('subplot00_gke_psi_yvalues'  , np.shape(gke_spl00_psi_y       ), dtype='f8', data=gke_spl00_psi_y       )
      h5f.create_dataset('subplot00_gke_psi_zvalues'  , np.shape(gke_spl00_psi_z       ), dtype='f8', data=gke_spl00_psi_z       )
      h5f.create_dataset('subplot00_wall_xvalues'     , np.shape(wall_x), dtype='f8', data=wall_x)
      h5f.create_dataset('subplot00_wall_yvalues'     , np.shape(wall_y), dtype='f8', data=wall_y)
      h5f.create_dataset('subplot01_xgc_nodes_xvalues', np.shape(xgc_spl01_nodes_x     ), dtype='f8', data=xgc_spl01_nodes_x     )
      h5f.create_dataset('subplot01_xgc_nodes_yvalues', np.shape(xgc_spl01_nodes_y     ), dtype='f8', data=xgc_spl01_nodes_y     )
#      h5f.create_dataset('subplot01_xgc_edges_constx' , np.shape(xgc_spl01_edges_constx), dtype='f8', data=xgc_spl01_edges_constx)
#      h5f.create_dataset('subplot01_xgc_edges_consty' , np.shape(xgc_spl01_edges_consty), dtype='f8', data=xgc_spl01_edges_consty)
#      h5f.create_dataset('subplot01_xgc_psi_xvalues'  , np.shape(xgc_spl01_psi_x       ), dtype='f8', data=xgc_spl01_psi_x       )
#      h5f.create_dataset('subplot01_xgc_psi_yvalues'  , np.shape(xgc_spl01_psi_y       ), dtype='f8', data=xgc_spl01_psi_y       )
#      h5f.create_dataset('subplot01_xgc_psi_zvalues'  , np.shape(xgc_spl01_psi_z       ), dtype='f8', data=xgc_spl01_psi_z       )
      h5f.create_dataset('subplot01_wall_xvalues'     , np.shape(wall_x), dtype='f8', data=wall_x)
      h5f.create_dataset('subplot01_wall_yvalues'     , np.shape(wall_y), dtype='f8', data=wall_y)
      h5f.close()

    fig_h.savefig(out_fig_dir+fig_file_name+figure_file_format)
    plt.close()
  else:
    plt.show()

#................................................................................#

if plot_den_omp_init_li863mg:
  #[ Plot radial profiles of density and temperature.
  fig_name = lcu.li863_prefix+'init_elc_den_omp'

  gke_data_file = gke_data_dir+'ltx_gkeyll_li863mg_init_elc_den_z1mid.h5'
  xgc_data_file = xgc_data_dir+'ltx_xgc_li863mg_init_elc_den_omp.h5'

  #[ Prepare figure.
  fig_prop = (6.0, 3.6)
  ax_pos   = [[0.15, 0.15, 0.83, 0.78],]
  fig_h    = plt.figure(figsize=fig_prop)
  ax_h     = [fig_h.add_axes(pos) for pos in ax_pos]

  gke_data = h5py.File(gke_data_file, "r")
  gke_spl00_line0_x = lcu.h5data_to_numpy_array(gke_data, 'subplot00_line0_xvalues')
  gke_spl00_line0_y = lcu.h5data_to_numpy_array(gke_data, 'subplot00_line0_yvalues')
  gke_data.close()

  xgc_data = h5py.File(xgc_data_file, "r")
  xgc_spl00_line0_x = lcu.h5data_to_numpy_array(xgc_data, 'subplot00_line0_xvalues')
  xgc_spl00_line0_y = lcu.h5data_to_numpy_array(xgc_data, 'subplot00_line0_yvalues')
  xgc_data.close()

  hpla = list()
  hpla.append(ax_h[0].plot(gke_spl00_line0_x, gke_spl00_line0_y, linestyle=gke_linestyle, color=gke_color))
  hpla.append(ax_h[0].plot(xgc_spl00_line0_x, xgc_spl00_line0_y, linestyle=xgc_linestyle, color=xgc_color))

  xmin = [min(np.amin(gke_spl00_line0_x),np.amin(xgc_spl00_line0_x)),]
  xmax = [max(np.amax(gke_spl00_line0_x),np.amax(xgc_spl00_line0_x)),]
  ymin = [min(np.amin(gke_spl00_line0_y),np.amin(xgc_spl00_line0_y)),]
  ymax = [max(np.amax(gke_spl00_line0_y),np.amax(xgc_spl00_line0_y)),]

  ax_h[0].set_ylabel(r'$n_e(\theta=0,t=0)~(\mathrm{m}^{-3})$', fontsize=lcu.xy_label_font_size)
  for i in range(len(ax_h)):
    ax_h[i].set_xlabel(r'$\psi_N$', fontsize=lcu.xy_label_font_size, labelpad=-2)
    ax_h[i].set_xlim(xmin[i], xmax[i])
    ax_h[i].set_ylim(0., ax_h[i].get_ylim()[1])
    lcu.set_tick_font_size(ax_h[0+i],lcu.tick_font_size)
    hmagx = ax_h[i].yaxis.get_offset_text().set_size(lcu.tick_font_size)
  # end
  ax_h[0].legend([hpla[0][0],hpla[1][0]],['Gkeyll','XGC'],fontsize=lcu.legend_font_size, frameon=False, loc='upper right')

  if output_figure_file:
    fig_file_name = output_prefix+fig_name

    if save_data:
      h5f = h5py.File(out_data_dir+fig_file_name+'.h5', "w")
      h5f.create_dataset('subplot00_gke_line0_xvalues', np.shape(gke_spl00_line0_x), dtype='f8', data=gke_spl00_line0_x)
      h5f.create_dataset('subplot00_gke_line0_yvalues', np.shape(gke_spl00_line0_y), dtype='f8', data=gke_spl00_line0_y)
      h5f.create_dataset('subplot00_xgc_line0_xvalues', np.shape(xgc_spl00_line0_x), dtype='f8', data=xgc_spl00_line0_x)
      h5f.create_dataset('subplot00_xgc_line0_yvalues', np.shape(xgc_spl00_line0_y), dtype='f8', data=xgc_spl00_line0_y)
      h5f.close()

    fig_h.savefig(out_fig_dir+fig_file_name+figure_file_format)
    plt.close()
  else:
    plt.show()

#................................................................................#

def plot_den_temp_omp(gke_data_file, xgc_data_file, fig_name, y_labels):
  #[ Plot density and temperature at the OMP.

  #[ Prepare figure.
  fig_prop = (12., 3.6)
  ax_pos   = [[0.08, 0.15, 0.25, 0.78],[0.41, 0.15, 0.25, 0.78],[0.74, 0.15, 0.25, 0.78],]
  fig_h     = plt.figure(figsize=fig_prop)
  ax_h      = [fig_h.add_axes(pos) for pos in ax_pos]

  gke_data = h5py.File(gke_data_file, "r")
  gke_spl00_line0_x = lcu.h5data_to_numpy_array(gke_data, 'subplot00_line0_xvalues')
  gke_spl00_line0_y = lcu.h5data_to_numpy_array(gke_data, 'subplot00_line0_yvalues')
  gke_spl01_line0_x = lcu.h5data_to_numpy_array(gke_data, 'subplot01_line0_xvalues')
  gke_spl01_line0_y = lcu.h5data_to_numpy_array(gke_data, 'subplot01_line0_yvalues')
  gke_spl02_line0_x = lcu.h5data_to_numpy_array(gke_data, 'subplot02_line0_xvalues')
  gke_spl02_line0_y = lcu.h5data_to_numpy_array(gke_data, 'subplot02_line0_yvalues')
  gke_data.close()

  xgc_data = h5py.File(xgc_data_file, "r")
  xgc_spl00_line0_x = lcu.h5data_to_numpy_array(xgc_data, 'subplot00_line0_xvalues')
  xgc_spl00_line0_y = lcu.h5data_to_numpy_array(xgc_data, 'subplot00_line0_yvalues')
  xgc_spl01_line0_x = lcu.h5data_to_numpy_array(xgc_data, 'subplot01_line0_xvalues')
  xgc_spl01_line0_y = lcu.h5data_to_numpy_array(xgc_data, 'subplot01_line0_yvalues')
  xgc_spl02_line0_x = lcu.h5data_to_numpy_array(xgc_data, 'subplot02_line0_xvalues')
  xgc_spl02_line0_y = lcu.h5data_to_numpy_array(xgc_data, 'subplot02_line0_yvalues')
  xgc_data.close()

  hpla, hplb, hplc = list(), list(), list()
  hpla.append(ax_h[0].plot(gke_spl00_line0_x, gke_spl00_line0_y, linestyle=gke_linestyle, color=gke_color))
  hplb.append(ax_h[1].plot(gke_spl01_line0_x, gke_spl01_line0_y, linestyle=gke_linestyle, color=gke_color))
  hplc.append(ax_h[2].plot(gke_spl02_line0_x, gke_spl02_line0_y, linestyle=gke_linestyle, color=gke_color))
  hpla.append(ax_h[0].plot(xgc_spl00_line0_x, xgc_spl00_line0_y, linestyle=xgc_linestyle, color=xgc_color))
  hplb.append(ax_h[1].plot(xgc_spl01_line0_x, xgc_spl01_line0_y, linestyle=xgc_linestyle, color=xgc_color))
  hplc.append(ax_h[2].plot(xgc_spl02_line0_x, xgc_spl02_line0_y, linestyle=xgc_linestyle, color=xgc_color))

  xmin = [
    min(np.amin(gke_spl00_line0_x),np.amin(xgc_spl00_line0_x)),
    min(np.amin(gke_spl01_line0_x),np.amin(xgc_spl01_line0_x)),
    min(np.amin(gke_spl02_line0_x),np.amin(xgc_spl02_line0_x)),
  ]
  xmax = [
    max(np.amax(gke_spl00_line0_x),np.amax(xgc_spl00_line0_x)),
    max(np.amax(gke_spl01_line0_x),np.amax(xgc_spl01_line0_x)),
    max(np.amax(gke_spl02_line0_x),np.amax(xgc_spl02_line0_x)),
  ]
  ymin = [
    min(np.amin(gke_spl00_line0_y),np.amin(xgc_spl00_line0_y)),
    min(np.amin(gke_spl01_line0_y),np.amin(xgc_spl01_line0_y)),
    min(np.amin(gke_spl02_line0_y),np.amin(xgc_spl02_line0_y)),
  ]
  ymax = [
    max(np.amax(gke_spl00_line0_y),np.amax(xgc_spl00_line0_y)),
    max(np.amax(gke_spl01_line0_y),np.amax(xgc_spl01_line0_y)),
    max(np.amax(gke_spl02_line0_y),np.amax(xgc_spl02_line0_y)),
  ]

  for i in range(len(ax_h)):
    ax_h[i].set_xlabel(r'$\psi_N$', fontsize=lcu.xy_label_font_size, labelpad=0)
    ax_h[i].set_ylabel(y_labels[i], fontsize=lcu.xy_label_font_size)
    ax_h[i].set_xlim(xmin[i], xmax[i])
    ax_h[i].set_ylim(0., ax_h[i].get_ylim()[1])
    lcu.set_tick_font_size(ax_h[i],lcu.tick_font_size)
    hmagx = ax_h[i].yaxis.get_offset_text().set_size(lcu.tick_font_size)
  plt.text(0.82, 0.88, r'(a)', fontsize=lcu.text_font_size, color='black', fontweight='regular', transform=ax_h[0].transAxes)
  plt.text(0.82, 0.88, r'(b)', fontsize=lcu.text_font_size, color='black', fontweight='regular', transform=ax_h[1].transAxes)
  plt.text(0.82, 0.88, r'(c)', fontsize=lcu.text_font_size, color='black', fontweight='regular', transform=ax_h[2].transAxes)
  ax_h[0].legend([hpla[0][0],hpla[1][0]],['Gkeyll','XGC'],fontsize=lcu.legend_font_size, frameon=False, loc='center left')

  if output_figure_file:
    fig_file_name = output_prefix+fig_name

    if save_data:
      h5f = h5py.File(out_data_dir+fig_file_name+'.h5', "w")
      h5f.create_dataset('subplot00_gke_line0_xvalues', np.shape(gke_spl00_line0_x), dtype='f8', data=gke_spl00_line0_x)
      h5f.create_dataset('subplot00_gke_line0_yvalues', np.shape(gke_spl00_line0_y), dtype='f8', data=gke_spl00_line0_y)
      h5f.create_dataset('subplot00_xgc_line0_xvalues', np.shape(xgc_spl00_line0_x), dtype='f8', data=xgc_spl00_line0_x)
      h5f.create_dataset('subplot00_xgc_line0_yvalues', np.shape(xgc_spl00_line0_y), dtype='f8', data=xgc_spl00_line0_y)
                                        
      h5f.create_dataset('subplot01_gke_line0_xvalues', np.shape(gke_spl01_line0_x), dtype='f8', data=gke_spl01_line0_x)
      h5f.create_dataset('subplot01_gke_line0_yvalues', np.shape(gke_spl01_line0_y), dtype='f8', data=gke_spl01_line0_y)
      h5f.create_dataset('subplot01_xgc_line0_xvalues', np.shape(xgc_spl01_line0_x), dtype='f8', data=xgc_spl01_line0_x)
      h5f.create_dataset('subplot01_xgc_line0_yvalues', np.shape(xgc_spl01_line0_y), dtype='f8', data=xgc_spl01_line0_y)
                                        
      h5f.create_dataset('subplot02_gke_line0_xvalues', np.shape(gke_spl02_line0_x), dtype='f8', data=gke_spl02_line0_x)
      h5f.create_dataset('subplot02_gke_line0_yvalues', np.shape(gke_spl02_line0_y), dtype='f8', data=gke_spl02_line0_y)
      h5f.create_dataset('subplot02_xgc_line0_xvalues', np.shape(xgc_spl02_line0_x), dtype='f8', data=xgc_spl02_line0_x)
      h5f.create_dataset('subplot02_xgc_line0_yvalues', np.shape(xgc_spl02_line0_y), dtype='f8', data=xgc_spl02_line0_y)
      h5f.close()

    fig_h.savefig(out_fig_dir+fig_file_name+figure_file_format)
    plt.close()
  else:
    plt.show()

if plot_den_temp_omp_init_li863mg:
  #[ Plot initial density and temperature at the OMP.
  fig_name = lcu.li863_prefix+'init_den_temp_omp'
  y_labels = [r'$n_e(\theta=0,t=0)$ (m$^{-3}$)', r'$T_e(\theta=0,t=0)$ (eV)', r'$T_i(\theta=0,t=0)$ (eV)',]

  gke_data_file = gke_data_dir+'ltx_gkeyll_li863mg_init_den_temp_z1mid.h5'
  #xgc_data_file = xgc_data_dir+'ltx_xgc_li863mg_init_den_temp_omp.h5'
  xgc_data_file = xgc_data_dir+'ltx103955-04_v0_xgc_init_moments_omp.h5'

  plot_den_temp_omp(gke_data_file, xgc_data_file, fig_name, y_labels)

#................................................................................#

if plot_den_temp_omp_final_li863mg:
  #[ Plot final density and temperature at the OMP.
  fig_name = lcu.li863_prefix+'final_den_temp_omp'
  y_labels = [
    r'$n_e(\theta=0)$ (m$^{-3}$)',
    r'$T_e(\theta=0)$ (eV)',
    r'$T_i(\theta=0)$ (eV)',
  ]

  gke_data_file = gke_data_dir+'ltx_gkeyll_li863mg_final_den_temp_z1mid.h5'
#  xgc_data_file = xgc_data_dir+'ltx_xgc_li863mg_final_den_temp_omp.h5'
  xgc_data_file = xgc_data_dir+'ltx103955-04_v0_xgc_final_moments_omp.h5'

  plot_den_temp_omp(gke_data_file, xgc_data_file, fig_name, y_labels)

#................................................................................#

if plot_moms_RZ_final_li863mg:

  xgc_mesh_file= xgc_data_dir+'ltx_xgc_li863mg_grid_RZ.h5'
  xgc_mesh_data = h5py.File(xgc_mesh_file, "r")
  rnodes      = lcu.h5data_to_numpy_array(xgc_mesh_data, 'subplot00_nodes_xvalues')
  znodes     = lcu.h5data_to_numpy_array(xgc_mesh_data, 'subplot00_nodes_yvalues')
  conn = lcu.h5data_to_numpy_array(xgc_mesh_data, 'subplot00_line_connections')
#  xgc_spl00_edges_constx = lcu.h5data_to_numpy_array(xgc_data, 'subplot00_edges_constx')
#  xgc_spl00_edges_consty = lcu.h5data_to_numpy_array(xgc_data, 'subplot00_edges_consty')
#  xgc_spl00_psi_x        = lcu.h5data_to_numpy_array(xgc_data, 'subplot00_psi_xvalues')
#  xgc_spl00_psi_y        = lcu.h5data_to_numpy_array(xgc_data, 'subplot00_psi_yvalues')
#  xgc_spl00_psi_z        = lcu.h5data_to_numpy_array(xgc_data, 'subplot00_psi_zvalues')
  xgc_mesh_data.close()

  triang = mtri.Triangulation(rnodes,znodes,conn)


  #[ Final n, upar, Tpar and Tperp profiles on R-Z plane.
  fig_name = lcu.li863_prefix+'final_moms_RZ'

  species = 'ion' #[ Species name, 'elc' or 'ion'

  gke_data_file = [
    gke_data_dir+'ltx_gkeyll_li863mg_final_'+species+'_den_RZ.h5',
    gke_data_dir+'ltx_gkeyll_li863mg_final_'+species+'_upar_RZ.h5',
    gke_data_dir+'ltx_gkeyll_li863mg_final_'+species+'_tpar_RZ.h5',
    gke_data_dir+'ltx_gkeyll_li863mg_final_'+species+'_tperp_RZ.h5',
  ]
  xgc_data_file = [
    xgc_data_dir+'ltx_xgc_li863mg_final_'+species+'_den_RZ.h5',
    xgc_data_dir+'ltx_xgc_li863mg_final_'+species+'_upar_RZ.h5',
    xgc_data_dir+'ltx_xgc_li863mg_final_'+species+'_tpar_RZ.h5',
    xgc_data_dir+'ltx_xgc_li863mg_final_'+species+'_tperp_RZ.h5',
  ]

  wall_file = '../experiment/LTXvessel.csv'

  xlabel = r'$R$ (m)'
  ylabel = r'$Z$ (m)'

  #[ Load data.
  gke_data = [h5py.File(gke_data_file[i], "r") for i in range(len(gke_data_file))]
  spl00_x = lcu.h5data_to_numpy_array(gke_data[0], 'subplot00_xvalues')
  spl00_y = lcu.h5data_to_numpy_array(gke_data[0], 'subplot00_yvalues')
  spl00_z = lcu.h5data_to_numpy_array(gke_data[0], 'subplot00_zvalues')
  spl01_x = lcu.h5data_to_numpy_array(gke_data[1], 'subplot00_xvalues')
  spl01_y = lcu.h5data_to_numpy_array(gke_data[1], 'subplot00_yvalues')
  spl01_z = lcu.h5data_to_numpy_array(gke_data[1], 'subplot00_zvalues')
  spl02_x = lcu.h5data_to_numpy_array(gke_data[2], 'subplot00_xvalues')
  spl02_y = lcu.h5data_to_numpy_array(gke_data[2], 'subplot00_yvalues')
  spl02_z = lcu.h5data_to_numpy_array(gke_data[2], 'subplot00_zvalues')
  spl03_x = lcu.h5data_to_numpy_array(gke_data[3], 'subplot00_xvalues')
  spl03_y = lcu.h5data_to_numpy_array(gke_data[3], 'subplot00_yvalues')
  spl03_z = lcu.h5data_to_numpy_array(gke_data[3], 'subplot00_zvalues')
  for i in range(len(gke_data_file)):
    gke_data[i].close()

  xgc_data = [h5py.File(xgc_data_file[i], "r") for i in range(len(xgc_data_file))]
  spl10_x = lcu.h5data_to_numpy_array(xgc_data[0], 'subplot00_xvalues')
  spl10_y = lcu.h5data_to_numpy_array(xgc_data[0], 'subplot00_yvalues')
  spl10_z = lcu.h5data_to_numpy_array(xgc_data[0], 'subplot00_zvalues')
  spl11_x = lcu.h5data_to_numpy_array(xgc_data[1], 'subplot00_xvalues')
  spl11_y = lcu.h5data_to_numpy_array(xgc_data[1], 'subplot00_yvalues')
  spl11_z = lcu.h5data_to_numpy_array(xgc_data[1], 'subplot00_zvalues')
  spl12_x = lcu.h5data_to_numpy_array(xgc_data[2], 'subplot00_xvalues')
  spl12_y = lcu.h5data_to_numpy_array(xgc_data[2], 'subplot00_yvalues')
  spl12_z = lcu.h5data_to_numpy_array(xgc_data[2], 'subplot00_zvalues')
  spl13_x = lcu.h5data_to_numpy_array(xgc_data[3], 'subplot00_xvalues')
  spl13_y = lcu.h5data_to_numpy_array(xgc_data[3], 'subplot00_yvalues')
  spl13_z = lcu.h5data_to_numpy_array(xgc_data[3], 'subplot00_zvalues')
  for i in range(len(xgc_data_file)):
    xgc_data[i].close()

  #[ Prepare figure.
  fig_prop  = (16., 8.)
  ax_pos    = [[0.08, 0.52, 0.20, 0.40],[0.29, 0.52, 0.20, 0.40],[0.50, 0.52, 0.20, 0.40],[0.71, 0.52, 0.20, 0.40],
               [0.08, 0.10, 0.20, 0.40],[0.29, 0.10, 0.20, 0.40],[0.50, 0.10, 0.20, 0.40],[0.71, 0.10, 0.20, 0.40], ]
  cbax_pos  = [[0.26, 0.10, 0.01, 0.82],[0.47, 0.10, 0.01, 0.82],[0.68, 0.10, 0.01, 0.82],[0.89, 0.10, 0.01, 0.82], ]
  fig_h     = plt.figure(figsize=fig_prop)
  ax_h      = [fig_h.add_axes(pos) for pos in ax_pos]
  cbax_h    = [fig_h.add_axes(pos) for pos in cbax_pos]

  for i in range(len(ax_h)):
    ax_h[i].set_aspect('equal')

  zmin, zmax = np.zeros(4), np.zeros(4)
  zmin[0], zmax[0] = min(np.amin(spl00_z),np.amin(spl10_z)), max(np.amax(spl00_z),np.amax(spl10_z))
  zmin[1], zmax[1] = min(np.amin(spl01_z),np.amin(spl11_z)), max(np.amax(spl01_z),np.amax(spl11_z))
  zmin[2], zmax[2] = min(np.amin(spl02_z),np.amin(spl12_z)), max(np.amax(spl02_z),np.amax(spl12_z))
  zmin[3], zmax[3] = min(np.amin(spl03_z),np.amin(spl13_z)), max(np.amax(spl03_z),np.amax(spl13_z))

  #[ Plot data.
  hpla, hplb, hplc, hpld = list(), list(), list(), list(),
  hpla.append(ax_h[0].pcolormesh(spl00_x, spl00_y, spl00_z, cmap='inferno', vmin=zmin[0], vmax=zmax[0]))
  hplb.append(ax_h[1].pcolormesh(spl01_x, spl01_y, spl01_z, cmap='inferno', vmin=zmin[1], vmax=zmax[1]))
  hplc.append(ax_h[2].pcolormesh(spl02_x, spl02_y, spl02_z, cmap='inferno', vmin=zmin[2], vmax=zmax[2]))
  hpld.append(ax_h[3].pcolormesh(spl03_x, spl03_y, spl03_z, cmap='inferno', vmin=zmin[3], vmax=zmax[3]))

  hple, hplf, hplg, hplh = list(), list(), list(), list(),
#  hple.append(ax_h[4+0].pcolormesh(spl10_x, spl10_y, spl10_z, cmap='inferno', vmin=zmin[0], vmax=zmax[0]))
#  hplf.append(ax_h[4+1].pcolormesh(spl11_x, spl11_y, spl11_z, cmap='inferno', vmin=zmin[1], vmax=zmax[1]))
#  hplg.append(ax_h[4+2].pcolormesh(spl12_x, spl12_y, spl12_z, cmap='inferno', vmin=zmin[2], vmax=zmax[2]))
#  hplh.append(ax_h[4+3].pcolormesh(spl13_x, spl13_y, spl13_z, cmap='inferno', vmin=zmin[3], vmax=zmax[3]))
  hple.append(ax_h[4+0].tripcolor(triang, spl10_z, cmap='inferno', vmin=zmin[0], vmax=zmax[0]))
  hplf.append(ax_h[4+1].tripcolor(triang, spl11_z, cmap='inferno', vmin=zmin[1], vmax=zmax[1]))
  hplg.append(ax_h[4+2].tripcolor(triang, spl12_z, cmap='inferno', vmin=zmin[2], vmax=zmax[2]))
  hplh.append(ax_h[4+3].tripcolor(triang, spl13_z, cmap='inferno', vmin=zmin[3], vmax=zmax[3]))

  #[ Plot wall
  wall_data = np.loadtxt(open(wall_file),delimiter=',')
  wall_h = list()
  for i in range(len(ax_h)):
    wall_h.append(ax_h[i].plot(wall_data[:,0],wall_data[:,1],color="grey"))

  #[ Add colorbar.
  hpl_gke = [hpla, hplb, hplc, hpld]
  hcb_all = [plt.colorbar(hpl_gke[i][0], ax=ax_h[i], cax=cbax_h[i]) for i in range(len(hpl_gke))]
  for i in range(len(hcb_all)):
    hcb_all[i].ax.tick_params(labelsize=lcu.tick_font_size)
    hcb_all[i].ax.yaxis.get_offset_text().set_fontsize(lcu.tick_font_size)

  for i in range(4):
    ax_h[4+i].set_xlabel(xlabel, fontsize=lcu.xy_label_font_size)
    ax_h[4+i].xaxis.get_offset_text().set_size(lcu.tick_font_size)
    plt.setp( ax_h[i].get_xticklabels(), visible=False)

  for i in range(3):
    plt.setp( ax_h[1+i].get_yticklabels(), visible=False)
    plt.setp( ax_h[5+i].get_yticklabels(), visible=False)

  ax_h[0].set_ylabel(ylabel, fontsize=lcu.xy_label_font_size, labelpad=-3)
  ax_h[0].yaxis.get_offset_text().set_size(lcu.tick_font_size)
  ax_h[4].set_ylabel(ylabel, fontsize=lcu.xy_label_font_size, labelpad=-3)
  ax_h[4].yaxis.get_offset_text().set_size(lcu.tick_font_size)
  txt_labels = ['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)',]
  for i in range(len(ax_h)):
    lcu.set_tick_font_size(ax_h[i],lcu.tick_font_size)
    plt.text(0.8, 0.9, txt_labels[i], fontsize=lcu.text_font_size, color='black', transform=ax_h[i].transAxes,)

  plt.text(-0.6, 0.34, r'Gkeyll', fontsize=32, color='black', transform=ax_h[0].transAxes, rotation=90,)
  plt.text(-0.6, 0.39, r'XGC', fontsize=32, color='black', transform=ax_h[4].transAxes, rotation=90,)
  plt.text(0.3, 1.05, r'$n_{:s}$ (m${:s}$)'.format(species[0],'^{-3}')       , fontsize=20, color='black', transform=ax_h[0].transAxes,)
  plt.text(0.3, 1.05, r'$u_{:s}$ (km/s)'.format('{\parallel '+species[0]+'}'), fontsize=20, color='black', transform=ax_h[1].transAxes,)
  plt.text(0.3, 1.05, r'$T_{:s}$ (eV)'.format('{\parallel '+species[0]+'}')  , fontsize=20, color='black', transform=ax_h[2].transAxes,)
  plt.text(0.3, 1.05, r'$T_{:s}$ (eV)'.format('{\perp '+species[0]+'}')      , fontsize=20, color='black', transform=ax_h[3].transAxes,)

  if output_figure_file:
    fig_file_name = output_prefix+fig_name

    if save_data:
      h5f = h5py.File(out_data_dir+fig_file_name+'.h5', "w")
      h5f.create_dataset('subplot00_xvalues', np.shape(spl00_x), dtype='f8', data=spl00_x)
      h5f.create_dataset('subplot00_yvalues', np.shape(spl00_y), dtype='f8', data=spl00_y)
      h5f.create_dataset('subplot00_zvalues', np.shape(spl00_z), dtype='f8', data=spl00_z)
      h5f.create_dataset('subplot01_xvalues', np.shape(spl01_x), dtype='f8', data=spl01_x)
      h5f.create_dataset('subplot01_yvalues', np.shape(spl01_y), dtype='f8', data=spl01_y)
      h5f.create_dataset('subplot01_zvalues', np.shape(spl01_z), dtype='f8', data=spl01_z)
      h5f.create_dataset('subplot02_xvalues', np.shape(spl02_x), dtype='f8', data=spl02_x)
      h5f.create_dataset('subplot02_yvalues', np.shape(spl02_y), dtype='f8', data=spl02_y)
      h5f.create_dataset('subplot02_zvalues', np.shape(spl02_z), dtype='f8', data=spl02_z)
      h5f.create_dataset('subplot03_xvalues', np.shape(spl03_x), dtype='f8', data=spl03_x)
      h5f.create_dataset('subplot03_yvalues', np.shape(spl03_y), dtype='f8', data=spl03_y)
      h5f.create_dataset('subplot03_zvalues', np.shape(spl03_z), dtype='f8', data=spl03_z)
      h5f.create_dataset('subplot10_xvalues', np.shape(spl10_x), dtype='f8', data=spl10_x)
      h5f.create_dataset('subplot10_yvalues', np.shape(spl10_y), dtype='f8', data=spl10_y)
      h5f.create_dataset('subplot10_zvalues', np.shape(spl10_z), dtype='f8', data=spl10_z)
      h5f.create_dataset('subplot11_xvalues', np.shape(spl11_x), dtype='f8', data=spl11_x)
      h5f.create_dataset('subplot11_yvalues', np.shape(spl11_y), dtype='f8', data=spl11_y)
      h5f.create_dataset('subplot11_zvalues', np.shape(spl11_z), dtype='f8', data=spl11_z)
      h5f.create_dataset('subplot12_xvalues', np.shape(spl12_x), dtype='f8', data=spl12_x)
      h5f.create_dataset('subplot12_yvalues', np.shape(spl12_y), dtype='f8', data=spl12_y)
      h5f.create_dataset('subplot12_zvalues', np.shape(spl12_z), dtype='f8', data=spl12_z)
      h5f.create_dataset('subplot13_xvalues', np.shape(spl13_x), dtype='f8', data=spl13_x)
      h5f.create_dataset('subplot13_yvalues', np.shape(spl13_y), dtype='f8', data=spl13_y)
      h5f.create_dataset('subplot13_zvalues', np.shape(spl13_z), dtype='f8', data=spl13_z)
      h5f.close()

    fig_h.savefig(out_fig_dir+fig_file_name+figure_file_format)
    plt.close()
  else:
    plt.show()

#................................................................................#
