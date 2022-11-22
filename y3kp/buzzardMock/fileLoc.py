#!/usr/bin/env python
import os

class FileLocs(object):
    def __init__(self, machine='nersc'):
        if machine == 'fnal':
            # self.halo_run_loc = '/bsuhome/hwu/scratch/Buzzard/buzzard-1.9.8/'
            self.halo_run_loc = '/data/des61.a/data/johnny/Buzzard/Buzzard_v2.0.0/y3_rm/'
            self.halo_run_fname = self.halo_run_loc + 'buzzard-1.9.8_3y3a_run_halos_lambda_chisq_mu_star.fit'
            self.profile_output_fname = self.halo_run_loc + 'buzzard-1.9.8_3y3a_run_halos_profiles.fit'
            
            ## mu-star catalog
            self.data_loc   = '/data/des61.a/data/johnny/Buzzard/Buzzard_v2.0.0/y3_rm/output/to_heidi/'
            self.data_fname = self.data_loc + 'buzzard_v2_lambda_gt20_mu-star_rhod.fits'
            
        if machine == 'nersc':
            self.halo_run_loc = '/global/cfs/cdirs/des/jesteves/data/buzzard/v1.9.8/y3_rm/'
            self.halo_run_fname = self.halo_run_loc + 'buzzard-1.9.8_3y3a_run_halos_lambda_chisq_mu_star.fit'
            self.profile_output_fname = self.halo_run_loc + 'buzzard-1.9.8_3y3a_run_halos_profiles.fit'
            
            ## mu-star catalog
            self.data_fname = self.halo_run_loc + 'output/to_heidi/buzzard_v2_lambda_gt20_mu-star_rhod.fits'
            
            ## boost-factor data y1
            self.boost_dir = '/global/cfs/cdirs/des/jesteves/data/boost_factor/'
            self.boost_tamas_y1 = self.boost_dir+'y1/tamas_profiles/'
        
        # mock catalogs
        self.mock_fname = self.halo_run_loc + 'mock_buzzard_like_y3_v0.fits'
        self.mock_boost_factor_1d = self.halo_run_loc + 'beta_table_zl_y1_like.npz'
        self.mock_boost_factor_2d = self.halo_run_loc + 'beta_table_zl_radii_y1_like.npz'