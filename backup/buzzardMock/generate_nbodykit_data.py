import numpy as np
from nbodykit.lab import *
from nbodykit import setup_logging, style
from astropy.io.fits import getdata
from astropy.table import Table
from scipy.interpolate import InterpolatedUnivariateSpline

from fileLoc import FileLocs
floc = FileLocs(machine='nersc')

print('halo file name: %s'%floc.mock_nbody_fname)
data = FITSCatalog(floc.mock_nbody_fname)

print('Random file name: %s'%floc.mock_random_fname)
randoms = FITSCatalog(floc.mock_random_fname)

cosmo = cosmology.Cosmology(h=0.7).match(Omega0_m=0.3)

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

print(data['z'])
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

# def dask_to_table(cat):
#     out = Table()
#     columns = list(cat.columns)
#     for col in columns:
#         out[col] = np.array(cat[col])
#     return out
# i = 0
# for zmin, zmax in zip(zmin_list, zmax_list):
#     d, r = redshift_slice([data, randoms], zmin, zmax)
#     tb1 = dask_to_table(d)
#     tb2 = dask_to_table(r)
#     tb1.write(floc.halo_run_loc+'data_%i.fits'%i, overwrite=True)
#     tb2.write(floc.halo_run_loc+'rnd_%i.fits'%i, overwrite=True)
#     i+=1