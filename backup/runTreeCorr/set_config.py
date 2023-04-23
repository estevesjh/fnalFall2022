""" Setup files and bins

- Set the lambda/redshift bins
- Set the treecorr param file
- Set the output file systems
"""
import os
import numpy as np

Rmin_phys_mpc = 5.
Rmax_phys_mpc = 35.
nbins_phys_mpc = 6

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
         'nbins': nbins_phys_mpc - 2,
         'verbose': 0
        }

# 3d metric
config3 = {
         'nbins': nbins_phys_mpc,
         'min_sep': Rmin_phys_mpc,
         'max_sep': Rmax_phys_mpc,
         'min_rpar': -100,
         'max_rpar': 100,
         'verbose': 0
        }
