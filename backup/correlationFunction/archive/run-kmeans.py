import numpy as np
import h5py
import healpy as hp
import matplotlib.pyplot as plt
from astropy.table import Table

from jackEstimatorKmeans import JackKniferKmeans

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
# JK Setup
Npatches = 300
fname_base = './data/xi_kmeans{}_jk_patches_{:05d}.npy'
fname_kmeans_centers = './data/kmeans%i_centers.npy'%(Npatches)

# Fname Setup
path = '/project/projectdirs/des/www/y3_cats/'
rm_fname = path+'y3_gold_2.2.1_wide_sofcol_run_redmapper_v0.5.1_redmagic_12_3_19.h5'

columns = ['mem_match_id','ra','dec','z_lambda','lambda_chisq']
path_rm = 'catalog/redmapper/lgt5'
path_ran = 'randoms/redmapper/lgt5'

############## START CODE ################
# load data
print('load data')
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

rm = Table(read_hdf5(rm_fname, 'catalog/redmapper/lgt5', columns=columns))
ran = Table(read_hdf5(rm_fname, 'randoms/redmapper/lgt5', columns=None))

# 20 times the rm catalog
Nran = int(20*len(rm))
ix = np.random.randint(len(ran), size=Nran)
ran = ran[ix]

print('Start JK Estimator')
jk =  JackKniferKmeans(rm['ra'], rm['dec'], Npatches, fname=fname_kmeans_centers)
jk.show_stats()
jk.write(fname_kmeans_centers)

# random group labels
ran_labels = jk.add_randoms(ran['ra'],ran['dec'])

print('Run Treecorr')
def get_angular_correlation(data, randoms, config=config):
    # load data and random catalogs
    cat = treecorr.Catalog(ra=data['ra'], dec=data['dec'], ra_units='degrees', dec_units='degrees')
    rcat = treecorr.Catalog(ra=randoms['ra'], dec=randoms['dec'], ra_units='degrees', dec_units='degrees')
    
    # process data catalog
    dd = treecorr.NNCorrelation(config)
    dd.process(cat)
    
    # process random catalog (it takes longer)
    rr = treecorr.NNCorrelation(config)
    rr.process(rcat)
    
    # return Landy & Szalay estimator
    xi, varxi = dd.calculateXi(rr)
    r = np.exp(dd.meanlogr)
    sig = np.sqrt(varxi)
    
    # clean memory
    rcat = cat = 0
    rr = dd = 0
    return r, xi

print('Running trecorr on patches')
for kk in range(Npatches):
    
    print('Patch %i'%kk)
    fname = fname_base.format(Npatches, kk)
    
    # find indices
    idx = jk.get_mask(kk)
    idx1= jk.get_mask(kk, ran_labels)
    
    # run correlation functions
    _r, _xi = get_angular_correlation(rm[idx], ran[idx1])
    
    # append results
    np.save(fname, _xi)
    print('Saved results: %s'%fname)
    print('')
    
# nohup python run.py >log 2>&1 &