import numpy as np
from nbodykit.lab import *
from nbodykit import setup_logging, style
from astropy.io.fits import getdata
from astropy.table import Table
from scipy.interpolate import InterpolatedUnivariateSpline
from fileLoc import FileLocs
floc = FileLocs(machine='nersc')
cosmo = cosmology.Cosmology(h=0.7).match(Omega0_m=0.3)

print('halo file name: %s'%floc.mock_nbody_fname)
data = FITSCatalog(floc.mock_nbody_fname)

print('Random file name: %s'%floc.mock_random_fname)
randoms = FITSCatalog(floc.mock_random_fname)

n_patches = 10

# bins
zbins = np.array([0.2, 0.32, 0.373, 0.51, 0.65])
zmin_list = np.array([0.2, 0.373, 0.51])
zmax_list = np.array([0.32, 0.51, 0.64])
zmeans = np.array([0.25, 0.44, 0.575])

def redshift_slice(jointCat, zmin, zmax):
    datas = []
    for cat in jointCat:
        mask = (cat['z']>=zmin)&(cat['z']<zmax)
        datas.append(cat[mask].copy())
    return datas

# add Cartesian position column
data['Position'] = transform.SkyToCartesian(data['ra'], data['dec'], data['z'], cosmo=cosmo)
randoms['Position'] = transform.SkyToCartesian(randoms['ra'], randoms['dec'], randoms['z'], cosmo=cosmo)

# Make FKP weighting
# the sky fraction, used to compute volume in n(z)
FSKY = 1. # a made-up value
zhist = RedshiftHistogram(randoms, FSKY, cosmo, redshift='z')
alpha = 1.0 * data.csize / randoms.csize
nofz = InterpolatedUnivariateSpline(zhist.bin_centers, alpha*zhist.nbar)

# add the n(z) columns to the FKPCatalog
randoms['NZ'] = nofz(randoms['z'])
data['NZ'] = nofz(data['z'])

data['FKPWeight'] = 1.0 / (1 + data['NZ'] * 1e4)
randoms['FKPWeight'] = 1.0 / (1 + randoms['NZ'] * 1e4)

def dask_to_table(cat):
    out = Table()
    columns = list(cat.columns)
    for col in columns:
        out[col] = np.array(cat[col])
    return out

def get_labels(data, random, n_patches=10):
    import treecorr
    cat = treecorr.Catalog(ra=data['ra'], dec=data['dec'], ra_units='deg', dec_units='deg', npatch=n_patches)
    labels = cat.patch

    rcat = treecorr.Catalog(ra=random['ra'], dec=random['ra'], ra_units='deg', dec_units='deg', patch_centers=cat.patch_centers)
    labels_rnd = rcat.patch
    
    data['labels'] = labels
    random['labels'] = labels_rnd
    return data, random

i = 0
fname_base = floc.halo_run_loc+'nbody_output/cat/%s_z%i_k%i.fits'
data = dask_to_table(data)
randoms = dask_to_table(randoms)
for zmin, zmax in zip(zmin_list, zmax_list):
    print('redshift bin: %.2f < z < %.2f'%(zmin, zmax))
    d, r = redshift_slice([data, randoms], zmin, zmax)
    d, r = get_labels(d, r, n_patches=n_patches)
    for k in n_patches:
        fname = fname_base%("data", i, k)
        fname2 = fname_base%("random", i, k)
        dk = d[d['labels']!=k] 
        rk = r[r['labels']!=k] 
        dk.write(fname, overwrite=True)
        rk.write(fname2, overwrite=True)

    i+=1