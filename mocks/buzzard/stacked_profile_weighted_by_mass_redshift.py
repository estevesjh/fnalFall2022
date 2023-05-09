#!/usr/bin/env python
import numpy as np
import sys
__author__ = 'Heidi Wu'
# hywu@boisestate.edu

# Johnny's version copied from Heidi

def stacked_profile_weighted_by_mass_redshift(lnM_select, z_select, profile_select, lnM_all, z_all, profile_all, dm=0.1, dz=0.05):
    #### set up the bins for mass and redshift
    min_m = min(lnM_all)-dm
    max_m = max(lnM_all)+dm
    min_z = min(z_all)-dz
    max_z = max(z_all)+dz

    m_bins = np.arange(min_m, max_m+dm, dm)
    z_bins = np.arange(min_z, max_z+dz, dz)
    nM = len(m_bins)-1
    nz = len(z_bins)-1
    #print('nM, nz', nM, nz)

    nr = np.shape(profile_select)[1]#rbp.nbins_phys_mpc
    profile_weighted = np.zeros(nr)
    cov_weighted = np.zeros((nr,nr))
    weight_norm = 0

    for iz in range(nz):
        z_lo = z_bins[iz]
        z_hi = z_bins[iz+1]
        for iM in range(nM):
            m_lo = m_bins[iM]
            m_hi = m_bins[iM+1]
            select_bin = (lnM_select >= m_lo)&(lnM_select < m_hi)&(z_select>=z_lo)&(z_select<z_hi)
            weight = len(lnM_select[select_bin]) * 1.
            weight_norm += weight
            select_all = (lnM_all >= m_lo)&(lnM_all < m_hi)&(z_all>=z_lo)&(z_all<z_hi)

            #pdf1_list[iz, iM] = weight
            #pdf2_list[iz, iM] = len(lnM_all[select_all])


            if weight > 0 and len(lnM_all[select_all]) > 0:
                profile_weighted += (np.mean(profile_all[select_all, :], axis=0)*weight)
                
    profile_weighted /= weight_norm

    #pdf1_list = np.concatenate(pdf1_list)
    #pdf2_list = np.concatenate(pdf2_list)
    #diff_list = pdf2_list - pdf1_list

    return profile_weighted

def cov_stacked_profile_weighted_by_mass_redshift(lnM_select, z_select, profile_select, lnM_all, z_all, profile_all, dm=0.1, dz=0.05):
    #### set up the bins for mass and redshift
    min_m = min(lnM_all)-dm
    max_m = max(lnM_all)+dm
    min_z = min(z_all)-dz
    max_z = max(z_all)+dz

    m_bins = np.arange(min_m, max_m+dm, dm)
    z_bins = np.arange(min_z, max_z+dz, dz)
    nM = len(m_bins)-1
    nz = len(z_bins)-1
    #print('nM, nz', nM, nz)

    nr = np.shape(profile_select)[1]#rbp.nbins_phys_mpc
    profile_weighted = np.zeros(nr)
    cov_weighted = np.zeros((nr,nr))
    weight_norm = 0

    for iz in range(nz):
        z_lo = z_bins[iz]
        z_hi = z_bins[iz+1]
        for iM in range(nM):
            m_lo = m_bins[iM]
            m_hi = m_bins[iM+1]
            select_bin = (lnM_select >= m_lo)&(lnM_select < m_hi)&(z_select>=z_lo)&(z_select<z_hi)
            weight = len(lnM_select[select_bin]) * 1.
            weight_norm += weight
            select_all = (lnM_all >= m_lo)&(lnM_all < m_hi)&(z_all>=z_lo)&(z_all<z_hi)

            #pdf1_list[iz, iM] = weight
            #pdf2_list[iz, iM] = len(lnM_all[select_all])


            if weight > 0 and len(lnM_all[select_all]) > 0:
                profile_weighted += (np.mean(profile_all[select_all, :], axis=0)*weight)
                
                # compute covariance directly
                # cov(x_i, x_j) = \frac{1}{nhalos-1} \sum_{nhalos} \left[x_i -\bar{x}\right]\left[x_j -\bar{x}\right]
                nhalos = np.count_nonzero(select_all)
                diff = profile_all[select_all, :] - np.mean(profile_all[select_all, :], axis=0)[np.newaxis,:]
                sum_cov = np.sum(np.array([np.outer(diff[ix], diff[ix]) for ix in range(nhalos)]),axis=0)
                cov_weighted += weight*sum_cov/(nhalos-1)

    profile_weighted /= weight_norm
    cov_weighted /= weight_norm
    return cov_weighted, profile_weighted