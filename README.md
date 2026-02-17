# Gkeyll & XGC simulations of LTX-&beta;

We compare simulations of the Lithium Tokamak eXperiment &beta; (LTX-&beta;) at the Princeton Plasma Physics Laboratory (PPPL) using the Gkeyll and XGC gyrokinetic solvers.

## Dependencies

At minimum one needs a Python environment with
- numpy
- matplotlib
- h5py

in order to manipulate reduced (plotted) data stored in HDF5.

### Gkeyll dependencies

If you wish to process raw Gkeyll data you need to install [postgkyl](https://github.com/ammarhakim/postgkyl).

### XGC dependencies

Put XGC dendencies here.

## How to use this repository

Raw data from Gkeyll and XGC code is too large to store in GitHub. So to compare results between these two codes we will follow these 2 steps:
1. In each of the `gkeyll/` and `xgc/` folders we postprocess the raw data. This produces code-specific plots, and HDF5 files with the reduced (plotted) data.
2. In `comp/` we read the Gkeyll and XGC HDF5 files and make plots of the results from both codes together.

Since one person may have access to Gkeyll data but not XGC data, or viceversa, this approach allows us to share the reduced data so that both users can plot Gkeyll and XGC (reduced) data together.

### Directory organization

In each of `gkeyll/` and `xgc/` place the script that post-processes raw data. Then place the plots and reduced rata into two folders:
- figures
- data

Then in `comp/` we place one or more scripts that plot reduced data from both the XGC and Gkeyll codes, and places figures in the `comp/figures/` folder.

### File name convention

The name of data and figure files in `gkeyll/` and `xgc/` should begin with `ltx_gkeyll_` and `ltx_xgc_`, respecitively.

Also, it is advised that HDF5 files with reduced data for 1D plots name the datasets with the format
- `subplot#$_line%_xvalues`
- `subplot#$_line%_yvalues`

where `#` and `$` are the indices indicating the subplot number (in case a plot as multiple subplots), and `%` is an index indicating which line of the plot the dataset is for. In the case of 2D plots, the associated HDF5 files should use the following formats for the names of the datasets:
- `subplot#$_xvalues`
- `subplot#$_yvalues`
- `subplot#$_zvalues`
