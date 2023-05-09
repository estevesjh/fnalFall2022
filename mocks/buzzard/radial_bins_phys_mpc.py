#!/usr/bin/env python
import numpy as np
Rmin_phys_mpc = 0.0323
Rmax_phys_mpc = 30
nbins_phys_mpc = 15

lnrp_bins_phys_mpc = np.linspace(np.log(Rmin_phys_mpc), np.log(Rmax_phys_mpc), nbins_phys_mpc+1)
rp_bins_phys_mpc = np.exp(lnrp_bins_phys_mpc)
rp_phys_mpc = np.sqrt(rp_bins_phys_mpc[:-1]*rp_bins_phys_mpc[1:])

# Shear Profile Binning
# theta - angle separation arcmin
theta_min_arcmin = 20/60.
theta_max_arcmin = 120.
nbins_theta = 15

lnrp_bins_theta = np.linspace(np.log(theta_min_arcmin), np.log(theta_max_arcmin), nbins_theta+1)
_theta_bins_arcmin = np.exp(lnrp_bins_theta)
theta_arcmin = np.sqrt(_theta_bins_arcmin[:-1]*_theta_bins_arcmin[1:])
