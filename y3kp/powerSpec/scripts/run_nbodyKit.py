#!/usr/bin/env python
from nbodykit.lab import *
from nbodykit import setup_logging, style
import os, fitsio

import sys
import numpy as np
from astropy.io.fits import getdata
from astropy.table import Table
import argparse

from nbodykit import CurrentMPIComm
comm = CurrentMPIComm.get()

# load file
from fileLoc import FileLocs
floc = FileLocs(machine='nersc')
path = floc.halo_run_loc
inputdir = path+'/nbody_output/input/'

setup_logging("info")
cosmo = cosmology.Cosmology(h=0.7).match(Omega0_m=0.3)
parser = argparse.ArgumentParser(description='Run Nbodykit in Npatches and estimate a JK covariance')
parser.add_argument('tag', default='xi', type=str, help='kind: correlation function or power spectrum')

# kind
# z_col = 'z'
z_col = 'z_noise'

# bin setup
# radii = np.logspace(np.log10(1.), np.log10(200), 200+1)
# radii = np.logspace(np.log10(0.3), np.log10(150), 50+1)
radii = np.arange(1.0, 251, 1.)

def runPower(infile1, infile2, outfile, col=z_col):
    data = FITSCatalog(infile1)
    rnd = FITSCatalog(infile2)
        
    # do projection
    data['Position'] = transform.SkyToCartesian(data['ra'], data['dec'], data[col], cosmo=cosmo)
    rnd['Position'] = transform.SkyToCartesian(rnd['ra'], rnd['dec'], rnd[col], cosmo=cosmo)
    
    fkp = FKPCatalog(data, rnd)
    
    # clear memory
    rnd = data = 0
    
    mesh = fkp.to_mesh(Nmesh=1800, nbar='NZ', fkp_weight='FKPWeight', window='tsc')
    result = ConvolvedFFTPower(mesh, poles=list(np.arange(0, 32+2, 2)), dk=0.0025, kmin=0.001, kmax=2.)
    result.save(outfile)
    pass

def run2PCF(infile1, infile2, outfile, col=z_col):
    data = FITSCatalog(infile1)
    rnd = FITSCatalog(infile2)

    # do projection
    data['Position'] = transform.SkyToCartesian(data['ra'], data['dec'], data[col], cosmo=cosmo)
    rnd['Position'] = transform.SkyToCartesian(rnd['ra'], rnd['dec'], rnd[col], cosmo=cosmo)

    cr = SurveyData2PCF('projected', data, rnd, radii, cosmo=cosmo,
                        pimax=200, ra='ra', dec='dec', redshift=col)
    cr.save(outfile)
    # save2PCF(outfile, cr)

def save2PCF(fname, pcf):
    r, corr = pcf.corr['r'], pcf.corr['corr']
    print('Save file: %s'%fname)
    np.savez(fname, r=r, mean=corr)

function_dict = {'xi': run2PCF, 'power':runPower}
# status, files1, files2, outfiles = get_file_list(z_bin, tag=tag)
def main(tag='xi'):
    name = '%s_%s'%(tag, z_col)
    # add all
    files1 = [inputdir+'mock_data_box%i.fits'%jj for jj in range(2)]
    files2 = [inputdir+'mock_rnd_box%i.fits'%jj for jj in range(2)]
    outfiles = [path+'nbody_output/%s_box%i_all.json'%(name,jj) for jj in range(2)]

    func = function_dict[tag]
    with TaskManager(cpus_per_task=640, use_all_cpus=True) as tm:
        print('Correlation Function')
        n_files = len(files1)
        for kk in tm.iterate(range(0, n_files)):
            print('%i - runing file: %s'%(kk, files1[kk]))
            func(files1[kk], files2[kk], outfiles[kk])
        
if __name__ == "__main__":
    args = parser.parse_args()
    main(tag=args.tag)


#cr = SurveyData2PCF('1d', data, rnd, radii, cosmo=cosmo,  ra='ra', dec='dec', redshift='z', show_progress=True)
#cr.save("2pcf_%i.json"%i)

# srun -n 4 python example.py
# print('Power Spectrum')
# fkp = FKPCatalog(data, rnd)
# mesh = fkp.to_mesh(Nmesh=1024, nbar='NZ', fkp_weight='FKPWeight', window='sym20')
# result = FFTPower(mesh, '1d', poles=[0,2,4], dk=0.005, kmin=0.005, kmax=2.)

# compute the multipoles
#r = ConvolvedFFTPower(mesh, poles=[0,2,4], dk=0.005, kmin=0.005, kmax=2.)
#r.save("Cpower_%i.json"%i)
