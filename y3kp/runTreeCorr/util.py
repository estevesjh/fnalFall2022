""" Utils

- Load the redMaPPer Y3 catalog on NERSC
"""
import os
import h5py
import numpy as np
from astropy.table import Table
from astropy.io.fits import getdata
import treecorr

# specify our cosmology
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)

# interpolate comoving distance
from scipy.interpolate import interp1d

h = cosmo.H(0).value/100.
zgrid = np.linspace(0., 10., 10000)
rcomov = cosmo.comoving_distance(np.array(zgrid)).value/h
get_rcomov = interp1d(zgrid, rcomov)

# Setup file loc
from fileLoc import FileLocs
floc = FileLocs(machine='nersc')
path = floc.halo_run_loc
indir = path+'/wp/input/'
outdir = path+'/wp/output/'

fname_base = indir+'mock_%s_box%i.fits'
fname_out = outdir+'mock_%s_box%03d.npz'

if not os.path.isdir(outdir):
    os.makedirs(outdir)

def load_mock(nsize=14, n_patches=20, is_3d=False):
    datas = []
    outfiles = []
    for i in range(nsize):
        fname1 = fname_base%('data', i)
        fname2 = fname_base%('rnd', i)
        
        data = Table(getdata(fname1))
        rnd = Table(getdata(fname2))
        if is_3d:
            data['rcomov'] = get_rcomov(np.array(data['z']))
            rnd['rcomov'] = get_rcomov(np.array(rnd['z']))
            
            c1 = treecorr.Catalog(ra=data['ra'], dec=data['dec'], r=data['rcomov'],
                                    npatch=n_patches, ra_units='degrees', dec_units='degrees')
            c2 = treecorr.Catalog(ra=rnd['ra'], dec=rnd['dec'], r=rnd['rcomov'],
                                    npatch=n_patches, ra_units='degrees', dec_units='degrees')

        else:
            c1 = treecorr.Catalog(ra=data['ra'], dec=data['dec'],
                                    npatch=n_patches, ra_units='degrees', dec_units='degrees')
            c2 = treecorr.Catalog(ra=rnd['ra'], dec=rnd['dec'],
                                    npatch=n_patches, ra_units='degrees', dec_units='degrees')


        datas.append([c1, c2])
        outfiles.append(fname_out%('treeCorr',i))

    return datas, outfiles