import numpy as np
import h5py
import healpy as hp
import matplotlib.pyplot as plt
from astropy.table import Table

from jackEstimator import JackKnifer

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
Nside = 8 
Nside_res = 64 # footprint
fname_base = './data/xi_{}_jk_patches_{:05d}.npy'

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

print('build footprint')
def radec_to_pix(ra,dec,nside=4):
    thetas,phis = np.radians(90-dec),np.radians(ra)
    return hp.ang2pix(nside, thetas, phis,nest=False)

npix = hp.nside2npix(Nside_res)
hpx_rm = np.array(radec_to_pix(rm['ra'],rm['dec'],nside=Nside_res))
hpx_ran = np.array(radec_to_pix(ran['ra'],ran['dec'],nside=Nside_res))

hpxmap = np.zeros(npix, dtype=np.int)
w, values = np.unique(hpx_ran,return_counts=True)
hpxmap[w] = 1

print('Start JK Estimator')
nside_jk = 8 #this is the low-resolution nside to define the JK masks. 
#total de patches a serem removidos
jk =  JackKnifer(nside_jk, hpxmap, frac_thr= 0.7)
npatches = jk.npatches
print('total Jk patches are', npatches)

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
for kk in range(npatches):
    print('Patch %i'%kk)
    fname = fname_base.format(Nside, kk)
    
    # find indices
    idx = jk.get_cat_indices(hpx_rm, kk)
    idx1= jk.get_cat_indices(hpx_ran, kk)
    
    # run correlation functions
    _r, _xi = get_angular_correlation(rm[idx], ran[idx1])
    
    # append results
    np.save(fname, _xi)
    print('Saved results: %s'%fname)
    print('')
    
# nohup python run.py 2> &1 > log &
