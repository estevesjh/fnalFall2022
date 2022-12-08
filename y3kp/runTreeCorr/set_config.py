""" Setup files and bins

- Set the lambda/redshift bins
- Set the treecorr param file
- Set the output file systems
"""
import os
import numpy as np

Rmin_phys_mpc = 3.
Rmax_phys_mpc = 20.
nbins_phys_mpc = 4

lnrp_bins_phys_mpc = np.linspace(np.log(Rmin_phys_mpc), np.log(Rmax_phys_mpc), nbins_phys_mpc+1)
rp_bins_phys_mpc = np.exp(lnrp_bins_phys_mpc)
rp_phys_mpc = np.sqrt(rp_bins_phys_mpc[:-1]*rp_bins_phys_mpc[1:])

# setup treecorr
config ={
         'dec_col': 'dec',
         'dec_units': 'degrees',
         'ra_col': 'ra',
         'ra_units': 'degrees',
         'sep_units': 'arcmin',
         #'max_sep': 125,
         #'min_sep': 2.5,
         'nbins': nbins_phys_mpc,
         'verbose': 0
        }

# 3d metric
config3 = {
         'nbins': 31,
         'max_sep': 200,
         'min_sep': 5,
         'min_rpar': -40,
         'max_rpar': nbins_phys_mpc,
         'verbose': 0
        }
