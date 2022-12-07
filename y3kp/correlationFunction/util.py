""" Utils

- Load the redMaPPer Y3 catalog on NERSC
"""
import h5py
import numpy as np
from astropy.table import Table
from astropy.io.fits import getdata
import treecorr

# import bins
from set_bins_files import lbd_bins, z_bins

# avoid sbatch error
import os
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
# export HDF5_USE_FILE_LOCKING='FALSE'

# interpolate comoving distance
from astropy.cosmology import Planck18 as cosmo
from scipy.interpolate import interp1d

h = cosmo.H(0).value/100.
zgrid = np.linspace(0., 10., 10000)
rcomov = cosmo.comoving_distance(np.array(zgrid)).value/h
get_rcomov = interp1d(zgrid, rcomov)

# Fname Setup
path = '/project/projectdirs/des/www/y3_cats/'
rm_fname = path+'y3_gold_2.2.1_wide_sofcol_run_redmapper_v0.5.1_redmagic_12_3_19.h5'

columns = ['mem_match_id','ra','dec','z_lambda','lambda_chisq']
path_rm = 'catalog/redmapper/lgt5'
path_ran = 'randoms/redmapper/lgt5'

from fileLoc import FileLocs
floc = FileLocs(machine='nersc')
path = floc.halo_run_loc
def load_mock(nz=3, n_patches=20):
    datas = []
    for i in range(nz):
        fname = floc.mock_nbody_fname
        fname2 = floc.mock_random_fname
        #fname = path+'data_%i.fits'%i
        #fname2 = path+'rnd_%i.fits'%i
        
        data = Table(getdata(fname))
        rnd = Table(getdata(fname2))
        data['rcomov'] = get_rcomov(np.array(data['z']))
        rnd['rcomov'] = get_rcomov(np.array(rnd['z']))
        
        c1 = treecorr.Catalog(ra=data['ra'], dec=data['dec'], r=data['rcomov'],
                               npatch=n_patches, ra_units='degrees', dec_units='degrees')
        c2 = treecorr.Catalog(ra=rnd['ra'], dec=rnd['dec'], r=rnd['rcomov'],
                               npatch=n_patches, ra_units='degrees', dec_units='degrees')
        
        datas.append([c1, c2])
    return datas

def load_data():
    rm = Table(read_hdf5(rm_fname, 'catalog/redmapper/lgt5', columns=columns))
    ran = Table(read_hdf5(rm_fname, 'randoms/redmapper/lgt5', columns=None))
    
    # compute rcomov
    rm['rcomov'] = get_rcomov(np.array(rm['z_lambda']))
    ran['rcomov'] = get_rcomov(np.array(ran['ztrue']))
    
    return rm, ran

def apply_bin_cut(lbd_bin_id, z_bin_id, rm, ran, n_patches=20, 
                  ran_factor=30, is_3d=False, is_all=False):
    # select lambda and redshift bin
    zmin, zmax = z_bins[z_bin_id], z_bins[z_bin_id+1]
    lmin, lmax = lbd_bins[lbd_bin_id], lbd_bins[lbd_bin_id+1]
    
    if is_all:
        lmin, lmax = lbd_bins[1], lbd_bins[-1]
    
    # redMaPPer mask
    lbd = np.array(rm['lambda_chisq'])
    zcls = np.array(rm['z_lambda'])
    
    rm_mask = (lbd>=lmin)&(lbd<lmax)
    rm_mask&= (zcls>=zmin)&(zcls<zmax)
    
    # Random Mask
    lbd_ran = np.array(ran['avg_lambdaout'])
    z_ran = np.array(ran['ztrue'])

    ran_mask = (lbd_ran>=lmin)&(lbd_ran<lmax)
    ran_mask&= (z_ran>=zmin)&(z_ran<zmax)
        
    data = Table(rm[["ra","dec","rcomov"]][rm_mask])
    random = Table(ran[["ra","dec","rcomov"]][ran_mask])
    
    # select ran_factor times the redMaPPer catalog; default 30 times 
    Nran_size = int(ran_factor*len(data))
    idx = np.random.randint(len(random), size=Nran_size)
    random = random[idx]

    # initialize treecorr catalogs
    if is_3d: # redshift space distance
        cat = treecorr.Catalog(ra=data['ra'], dec=data['dec'], r=data['rcomov'],
                               npatch=n_patches, ra_units='degrees', dec_units='degrees')

        rcat = treecorr.Catalog(ra=random['ra'], dec=random['dec'], r=random['rcomov'],
                                patch_centers=cat.patch_centers, ra_units='degrees', dec_units='degrees')
    else: # angular distance
        cat = treecorr.Catalog(ra=data['ra'], dec=data['dec'], npatch=n_patches, 
                           ra_units='degrees', dec_units='degrees')

        rcat = treecorr.Catalog(ra=random['ra'], dec=random['dec'], patch_centers=cat.patch_centers,
                            ra_units='degrees', dec_units='degrees')
        
    return cat, rcat

def read_hdf5(fname, path, columns=None):
    """Read the hdf5 files for a given path
    if columns is None read all columns
    """
    h5  = h5py.File(fname,'r')
    h5group = h5[path]
    if columns is None: columns = list(h5group.keys())
    
    out = dict()
    for col in columns:
        out[col] = h5group[col][:]
    h5.close()
    return out

# TBD: If needed
# Load Files
# class LoadClusterData(object):
#     """ Load the redMaPPer Y3 cluster sample
#     """
#     columns = ['mem_match_id','ra','dec','z_lambda','lambda_chisq']
#     path_rm = 'catalog/redmapper/lgt5'
#     path_ran = 'randoms/redmapper/lgt5'

#     def __init__(self):
#         pass
    
#     def set_fnames(self, env='NERSC'):
#         if env=='NERSC':
#             self.path = '/project/projectdirs/des/www/y3_cats/'
#             self.fname = self.path+'y3_gold_2.2.1_wide_sofcol_run_redmapper_v0.5.1_redmagic_12_3_19.h5'
#         else:
#             print('Env Error: please go to nersc to load this file')
