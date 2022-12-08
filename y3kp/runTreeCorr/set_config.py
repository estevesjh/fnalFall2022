""" Setup files and bins

- Set the lambda/redshift bins
- Set the treecorr param file
- Set the output file systems
"""
import os
import numpy as np

# setup treecorr
config ={
         'dec_col': 'dec',
         'dec_units': 'degrees',
         'ra_col': 'ra',
         'ra_units': 'degrees',
         'sep_units': 'arcmin',
         'max_sep': 125,
         'min_sep': 2.5,
         'nbins': 20,
         'verbose': 0
        }

# 3d metric
config3 = {
         'nbins': 31,
         'max_sep': 200,
         'min_sep': 5,
         'min_rpar': -40,
         'max_rpar': 40,
         'verbose': 0
        }

# setup output file system
class SetupFiles(object):
    root = './'
    
    def __init__(self, tag, Npatches):
        self.Npatches = Npatches
        self.root_out = SetupFiles.root+'output/'
        self.root_dir = self.root_out+'kmeans%i/'%(self.Npatches)
        self.root_tmp = self.root_out+'tmp/'
        
        # initialize variables
        self.setup_files(tag)
        
        # make dirs
        self.make_roots()
        
    def setup_files(self, tag):
        self.label = tag
        self.outfile = self.root_out+self.label+'.npz'
        self.tmp_fname_base = self.root_out+self.label+'_jk_patches_{:05d}.npy'
        self.fname_kmeans_centers= self.root_tmp+self.label+'_centers.npy'
        
        # create root directories
        self.make_roots()
        
    def get_outfile(self, lbd_bin, z_bin):
        fname = '%s_l%i-z%i.npz'%(self.label, lbd_bin, z_bin)
        return self.root_dir+fname
    
    def make_roots(self):
        self.make_dir(SetupFiles.root)
        self.make_dir(self.root_out)
        self.make_dir(self.root_dir)
        self.make_dir(self.root_tmp)
        
    def make_dir(self, dirname):
        if not os.path.isdir(dirname):
            os.mkdir(dirname)