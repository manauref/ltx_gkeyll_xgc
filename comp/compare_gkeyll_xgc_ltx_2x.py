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
import sys #[ For error exit.
#[ Append path to utilities folder.
sys.path.insert(0, '../util/')
import ltx_common_util as lcu

plot_den_omp_init_li863mg = True  #[ Initial n at outboard midplane (OMP) of 863 mg simulation.
plot_den_temp_omp_final = False  #[ Final n and T profiles at OMP.

gke_data_dir = '../gkeyll/data/' #[ Location of reduced Gkeyll data.
xgc_data_dir = '../xgc/data/' #[ Location of reduced XGC data.

out_data_dir  = './data/'
out_fig_dir   = './figures/'
output_prefix = 'ltx_sims_'

output_figure_file = True     #[ Output a figure file?.
figure_file_format = '.png'    #[ Can be .png, .pdf, .ps, .eps, .svg.
save_data          = True    #[ Indicate whether to save data in plot to HDF5 file.

#[ ............... End of user inputs (MAYBE) ..................... ]#

#[ Default colors for Gkeyll and XGC.
gke_color = lcu.default_colors[0]
xgc_color = lcu.default_colors[1]
#[ Default line styles for Gkeyll and XGC.
gke_linestyle = lcu.default_line_styles[0]
xgc_linestyle = lcu.default_line_styles[1]

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
      h5f.create_dataset('gke_subplot00_line0_xvalues', np.shape(gke_spl00_line0_x), dtype='f8', data=gke_spl00_line0_x)
      h5f.create_dataset('gke_subplot00_line0_yvalues', np.shape(gke_spl00_line0_y), dtype='f8', data=gke_spl00_line0_y)
      h5f.create_dataset('xgc_subplot00_line0_xvalues', np.shape(xgc_spl00_line0_x), dtype='f8', data=xgc_spl00_line0_x)
      h5f.create_dataset('xgc_subplot00_line0_yvalues', np.shape(xgc_spl00_line0_y), dtype='f8', data=xgc_spl00_line0_y)
      h5f.close()

    fig_h.savefig(out_fig_dir+fig_file_name+figure_file_format)
    plt.close()
  else:
    plt.show()

#................................................................................#


if plot_den_temp_omp_final:
  #[ Plot radial profiles of density and temperature.
  fig_name = 'nT_omp'

  gke_data_file = gke_data_dir+'ltx_gkeyll_elc_den_z1mid.h5'
  xgc_data_file = xgc_data_dir+'ltx_xgc_elc_den_omp.h5'

  #[ Prepare figure.
  figProp1 = (12., 3.6)
  ax1Pos   = [[0.095, 0.15, 0.25, 0.78],[0.4075, 0.15, 0.25, 0.78],[0.74, 0.15, 0.25, 0.78],]
  fig1     = plt.figure(figsize=figProp1)
  ax1      = [fig1.add_axes(pos) for pos in ax1Pos]

  gke_data = h5py.File(gke_data_file, "r")
  gke_spl00_line0_x = h5data_to_numpy_array(gke_data, 'subplot00_line0_xvalues')
  gke_spl00_line0_y = h5data_to_numpy_array(gke_data, 'subplot00_line0_yvalues')
  gke_spl01_line0_x = h5data_to_numpy_array(gke_data, 'subplot01_line0_xvalues')
  gke_spl01_line0_y = h5data_to_numpy_array(gke_data, 'subplot01_line0_yvalues')
  gke_spl02_line0_x = h5data_to_numpy_array(gke_data, 'subplot02_line0_xvalues')
  gke_spl02_line0_y = h5data_to_numpy_array(gke_data, 'subplot02_line0_yvalues')
  gke_data.close()

  xgc_data = h5py.File(xgc_data_file, "r")
  xgc_spl00_line0_x = h5data_to_numpy_array(xgc_data, 'subplot00_line0_xvalues')
  xgc_spl00_line0_y = h5data_to_numpy_array(xgc_data, 'subplot00_line0_yvalues')
  xgc_spl01_line0_x = h5data_to_numpy_array(xgc_data, 'subplot01_line0_xvalues')
  xgc_spl01_line0_y = h5data_to_numpy_array(xgc_data, 'subplot01_line0_yvalues')
  xgc_spl02_line0_x = h5data_to_numpy_array(xgc_data, 'subplot02_line0_xvalues')
  xgc_spl02_line0_y = h5data_to_numpy_array(xgc_data, 'subplot02_line0_yvalues')
  xgc_data.close()

  hpl1a, hpl1b, hpl1c = list(), list(), list()
  hpl1a.append(ax1[0].plot(gke_spl00_line0_x, gke_spl00_line0_y, linestyle=gke_linestyle, color=gke_color))
  hpl1b.append(ax1[1].plot(gke_spl01_line0_x, gke_spl01_line0_y, linestyle=gke_linestyle, color=gke_color))
  hpl1c.append(ax1[2].plot(gke_spl02_line0_x, gke_spl02_line0_y, linestyle=gke_linestyle, color=gke_color))
  hpl1a.append(ax1[0].plot(xgc_spl00_line0_x, xgc_spl00_line0_y, linestyle=xgc_linestyle, color=xgc_color))
  hpl1b.append(ax1[1].plot(xgc_spl01_line0_x, xgc_spl01_line0_y, linestyle=xgc_linestyle, color=xgc_color))
  hpl1c.append(ax1[2].plot(xgc_spl02_line0_x, xgc_spl02_line0_y, linestyle=xgc_linestyle, color=xgc_color))

  ymin = [
    np.amin([gke_spl00_line0_y,xgc_spl00_line0_y]),
    np.amin([gke_spl01_line0_y,xgc_spl01_line0_y]),
    np.amin([gke_spl02_line0_y,xgc_spl02_line0_y]),
  ]
  ymax = [
    np.amax([gke_spl00_line0_y,xgc_spl00_line0_y]),
    np.amax([gke_spl01_line0_y,xgc_spl01_line0_y]),
    np.amax([gke_spl02_line0_y,xgc_spl02_line0_y]),
  ]

  ax1[0].set_ylabel(r'$\left\langle n_e\right\rangle_{\theta,z,t}~(\mathrm{m}^{-3})$', fontsize=xyLabelFontSize)
  ax1[1].set_ylabel(r'$\left\langle T_e\right\rangle_{\theta,z,t}$ (eV)', fontsize=xyLabelFontSize)
  ax1[2].set_ylabel(r'$\left\langle T_i\right\rangle_{\theta,z,t}$ (eV)', fontsize=xyLabelFontSize)
  for i in range(3):
    ax1[i].set_xlabel(r'$\psi_N$', fontsize=xyLabelFontSize, labelpad=-2)
    ax1[i].set_xlim(ymin[i], ymax[i])
    ax1[i].set_ylim(0., ax1[i].get_ylim()[1])
    lcu.set_tick_font_size(ax1[0+i],tickFontSize)
    hmagx = ax1[i].yaxis.get_offset_text().set_size(tickFontSize)
  plt.text(0.82, 0.88, r'(a)', fontsize=textFontSize, color='black', fontweight='regular', transform=ax1[0].transAxes)
  plt.text(0.82, 0.88, r'(b)', fontsize=textFontSize, color='black', fontweight='regular', transform=ax1[1].transAxes)
  plt.text(0.82, 0.88, r'(c)', fontsize=textFontSize, color='black', fontweight='regular', transform=ax1[2].transAxes)
  ax1[0].legend([hpl1a[0][0],hpl1a[1][0]],['Gkeyll','XGC'],fontsize=legendFontSize, frameon=False, loc='center right')

  if output_figure_file:
    fig1.savefig(out_fig_dir+output_prefix+fig_name+figure_file_format)
    plt.close()
  else:
    plt.show()

#................................................................................#

