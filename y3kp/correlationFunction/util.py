import os
import numpy as np
import h5py
import healpy as hp
import matplotlib.pyplot as plt
from astropy.table import Table

from jackEstimatorKmeans import JackKniferKmeans

# Fname Setup
path = '/project/projectdirs/des/www/y3_cats/'
rm_fname = path+'y3_gold_2.2.1_wide_sofcol_run_redmapper_v0.5.1_redmagic_12_3_19.h5'

columns = ['mem_match_id','ra','dec','z_lambda','lambda_chisq']
path_rm = 'catalog/redmapper/lgt5'
path_ran = 'randoms/redmapper/lgt5'

# Tree Corr Setup
import treecorr
config ={
         'dec_col': 'dec',
         'dec_units': 'degrees',
         'ra_col': 'ra',
         'ra_units': 'degrees',
         'sep_units': 'arcmin',
         'max_sep': 100,
         'min_sep': 1,
         'nbins': 20,
         'verbose': 0
        }

def load_data():
    rm = Table(read_hdf5(rm_fname, 'catalog/redmapper/lgt5', columns=columns))
    ran = Table(read_hdf5(rm_fname, 'randoms/redmapper/lgt5', columns=None))

    # 20 times the rm catalog
    Nran = int(20*len(rm))
    ix = np.random.randint(len(ran), size=Nran)
    ran = ran[ix]
    return rm, ran

def read_hdf5(fname, path, columns=None):
    """Read the hdf5 files for a given path
    if columns is None read all columns
    """
    h5  = h5py.File(fname,'r+')
    h5group = h5[path]
    if columns is None: columns = list(h5group.keys())
    
    out = dict()
    for col in columns:
        out[col] = h5group[col][:]
    
    return out

class SetupFiles(object):
    root = './'
    root_out = root+'output/'
    root_tmp = root+'tmp/'
    
    def __init__(self, tag, Npatches):
        self.label = tag+'_kmeans%i'%Npatches
        self.outfile = SetupFiles.root_out+self.label+'.npz'
        self.tmp_fname_base = SetupFiles.root_tmp+self.label+'_jk_patches_{:05d}.npy'
        self.fname_kmeans_centers= SetupFiles.root_tmp+self.label+'_centers.npy'
        
        # create root directories
        self.make_roots()
        
    def make_roots(self):
        self.make_dir(SetupFiles.root)
        self.make_dir(SetupFiles.root_out)
        self.make_dir(SetupFiles.root_tmp)
        
    def make_dir(self, dirname):
        if not os.path.isdir(dirname):
            os.mkdir(dirname)