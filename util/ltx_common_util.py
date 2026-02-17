#[ ........................................................... ]#
#[
#[ Some common utilities used for postprocessing.
#[
#[ Manaure Francisquez & George Wilkie.
#[ February 2026.
#[
#[ ........................................................... ]#

import numpy as np

#[ Prefixes in filenames for 863 mg and passivated Li shots.
li863_prefix = 'li863mg_'
lipass_prefix = 'liPass_'

eV        = 1.602176487e-19
me, mp    = 9.10938215e-31, 1.672621637e-27
charge_elc, charge_ion = -eV, eV
mass_elc, mass_ion = me, mp

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

def set_tick_font_size(axIn, fontSizeIn):
  #[ Set the font size of the ticks to a given size.
  #[   - axIn: axis handle to modify the font size in.
  #[   - fontSizeIn: desired font size.
  axIn.tick_params(axis='both', which='major', labelsize=fontSizeIn)

def rho_psi(psi, psi_lcfs, psi_axis):
  #[ Normalized poloidal flux coordinate rho.
  #[   - psi: poloidal flux.
  #[   - psi_lcfs: poloidal flux at the LCFS.
  #[   - psi_axis: poloidal flux at the exis.
  return np.sqrt((psi - psi_axis) / (psi_lcfs - psi_axis))

def psi_rho(rho, psi_lcfs, psi_axis):
  #[ Poloidal flux psi in terms of normalized poloidal flux rho.
  #[   - rhp: Normalized poloidal flux coordinate rho.
  #[   - psi_lcfs: poloidal flux at the LCFS.
  #[   - psi_axis: poloidal flux at the exis.
  return np.power(rho,2)*(psi_lcfs - psi_axis) + psi_axis

def psi_N(psi, psi_lcfs, psi_axis):
  #[ Normalized psi.
  #[   - psi: poloidal flux.
  #[   - psi_lcfs: poloidal flux at the LCFS.
  #[   - psi_axis: poloidal flux at the exis.
  return (psi - psi_axis) / (psi_lcfs - psi_axis)

def h5data_to_numpy_array(h5file_handle, dataset_name):
  #[ Convert an HDF5 data set to a numpy array.
  #[   - h5file_handle: handle to open HDF5 file.
  #[   - dataset_name: name of dataset in HDF5 file to conver to numpy array.
  h5dat = h5file_handle[dataset_name] 
  nparr = np.zeros(h5dat.shape, dtype=h5dat.dtype)
  h5dat.read_direct(nparr)
  return nparr
