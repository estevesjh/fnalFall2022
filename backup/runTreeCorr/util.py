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

# from make_data import zmeans
zmeans = np.array([0.25, 0.44, 0.575])
alpha = 50

# angular conversion factor
rad2deg = 180/np.pi
conv_factor = {"radians": 1.,
               "degrees": 1.*rad2deg, 
               "arcmin" : 60.*rad2deg, 
               "arcsec" : 360.*rad2deg}

h = cosmo.H(0).value/100.
zgrid = np.linspace(0., 1., 1000)
rcomov = cosmo.comoving_distance(zgrid).value
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
    zvalues = []
    for i in range(nsize):
        fname1 = fname_base%('data', i)
        fname2 = fname_base%('rnd', i)
        fname = fname_out%('kmeans%i'%n_patches,i)

        data = Table(getdata(fname1))
        rnd = Table(getdata(fname2))

        cut = np.random.randint(len(rnd), size=int(alpha*len(data)))
        rnd = rnd[cut]
        
        if is_3d:
            fname = fname_out%('kmeans%i_3d'%n_patches,i)
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

        ij = np.argmin(np.abs(zmeans-np.mean(data['z'])))

        datas.append([c1, c2])
        outfiles.append(fname)
        zvalues.append(np.mean(zmeans[ij]))

    return datas, outfiles, zvalues

def theta_converter(r,z):
    # radii or z needs to be a float number 
    theta_vec = r/cosmo.angular_diameter_distance(z).value
    return theta_vec*conv_factor["arcmin"]

def convert_profile(z, profile):
    theta_vec = theta_converter(radii, z)
    interp_gamma_t = interp1d(theta_vec, profile, bounds_error=False)(theta)
    return interp_gamma_t