import os
import numpy as np
from astropy.io.fits import getdata
from astropy.table import Table
from scipy.interpolate import interp1d

# specify our cosmology
from astropy.cosmology import FlatLambdaCDM
cosmo0 = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)

# Setup file loc
from fileLoc import FileLocs
floc = FileLocs(machine='nersc')
path = floc.halo_run_loc

# load data
mock = Table(getdata(floc.mock_fname))
random = Table(getdata(floc.mock_random_fname))

data = Table(mock[['RA', 'DEC', 'Z']])
data.rename_columns(['RA','DEC','Z'], ['ra','dec','z'])

sigma_z0 = 0.01
def add_noise(z):
    return z+sigma_z0*(1+z)*np.random.normal(size=z.size)

data['z_noise'] = add_noise(data['z'])
random['z_noise'] = add_noise(random['z'])

# compute FKPWeight
### Nz
# redshit vector
fsky = 1/8. # DES footprint
z = np.where(data['z']<0.2, np.nan, data['z'])
z = np.where(z>0.65, np.nan, z)

zbins = np.linspace(0.001, 1.5, 100)
zm   = zbins[:-1] + np.diff(zbins)
vbins = fsky*(4/3.)*np.pi*(cosmo0.comoving_distance(zbins[1:])**3-cosmo0.comoving_distance(zbins[:-1])**3)
nz  = np.histogram(z, bins=zbins)[0]/vbins
nofz = interp1d(zm, nz)

# add the n(z) columns to the FKPCatalog
data['NZ'] = nofz(data['z'])
random['NZ'] = nofz(random['z'])

data['FKPWeight'] = 1.0 / (1 + data['NZ'] * 1e4)
random['FKPWeight'] = 1.0 / (1 + random['NZ'] * 1e4)

# do box separation
# bins
zmin_list = np.array([0.2, 0.51])
zmax_list = np.array([0.32, 0.64])
zmeans = np.array([0.25, 0.575])

def get_mask(x, xmin, xmax):
    return (x<=xmax)&(x>=xmin)

datas = []
randoms = []
for zmin, zmax in zip(zmin_list, zmax_list):
    mask = get_mask(data['z'], zmin, zmax)
    maskr = get_mask(random['z'], zmin, zmax)
    rnd = random[maskr]
    #rnd = rnd[np.random.randint(len(rnd), size=int(40*len(data[mask])))]
    
    datas.append(data[mask])
    randoms.append(rnd)

# save data
outdir = path+'/nbody_output/input/'

if not os.path.isdir(outdir):
    os.makedirs(outdir)

fname_base = outdir+'mock_%s_box%i.fits'
for i in range(len(datas)):
    fname = fname_base%('data', i)
    fnamer = fname_base%('rnd', i)
    datas[i].write(fname, overwrite=True)
    randoms[i].write(fnamer, overwrite=True)
    print('Saved file: %s'%fname)
    print('Saved file: %s'%fnamer)
    print()