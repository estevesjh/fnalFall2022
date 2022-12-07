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

# from set_bins_files import SetupFiles
# load file
from fileLoc import FileLocs
floc = FileLocs(machine='nersc')
path = floc.halo_run_loc
fname_base = path+'nbody_output/tmp/%s_mock_z%i'

setup_logging("info")
cosmo = cosmology.Cosmology(h=0.7).match(Omega0_m=0.3)
parser = argparse.ArgumentParser(description='Run Nbodykit in Npatches and estimate a JK covariance')
parser.add_argument('z_bin', default=1, type=int, help='redshift bin id')
parser.add_argument('tag', default='xi', type=str, help='kind: correlation function or power spectrum')

# bin setup
# radii = np.logspace(np.log10(1.), np.log10(200), 200+1)
radii = np.linspace(0.1, 200, 200+1)

def runPower(infile1, infile2, outfile):
    data = FITSCatalog(infile1)
    rnd = FITSCatalog(infile2)
    fkp = FKPCatalog(data, rnd)
    mesh = fkp.to_mesh(Nmesh=1024, nbar='NZ', fkp_weight='FKPWeight', window='sym20')
    result = ConvolvedFFTPower(mesh, poles=[0,2,4], dk=0.005, kmin=0.005, kmax=2.)
    result.save(outfile)

def run2PCF(infile1, infile2, outfile):
    data = FITSCatalog(infile1)
    rnd = FITSCatalog(infile2)
    cr = SurveyData2PCF('projected', data, rnd, radii, cosmo=cosmo,  
                        pimax=200, ra='ra', dec='dec', redshift='z')
    cr.save(outfile)
    # save2PCF(outfile, cr)

def save2PCF(fname, pcf):
    r, corr = pcf.corr['r'], pcf.corr['corr']
    print('Save file: %s'%fname)
    np.savez(fname, r=r, mean=corr)

def get_file_list(z_bin, tag='xi'):
    import glob
    files1 = glob.glob(fname_base%('data', z_bin)+'*')
    files2 = glob.glob(fname_base%('rnd', z_bin)+'*')
    n_patches = len(files1)
    if n_patches>0:
        outfiles = [fname_base%(tag, z_bin)+"_k%i.json"%kk for kk in range(n_patches)]
        return True, files1, files2, outfiles
    else:
        return False, None, None, None

function_dict = {'xi': run2PCF, 'power':runPower}

def main(z_bin, tag='xi'):
    # status, files1, files2, outfiles = get_file_list(z_bin, tag=tag)
    
    # add all
    files1, files2, outfiles = [], [], []
    status = True
    files1.append(path+'data_%i.fits'%z_bin)
    files2.append(path+'rnd_%i.fits'%z_bin)
    outfiles.append(path+'nbody_output/%s_z%i_all.json'%(tag,z_bin))

    func = function_dict[tag]
    # status, files1, files2, outfiles = get_data(my_rank, tag, z_bin, comm)
    if status:
        with TaskManager(cpus_per_task=256, use_all_cpus=True) as tm:
            #print('Correlation Function')
            n_patches = len(files1)
            for kk in tm.iterate(range(0, n_patches)):
                print('%i - runing file: %s'%(kk, files1[kk]))
                func(files1[kk], files2[kk], outfiles[kk])
                # runPower(files1[kk], files2[kk], outfiles[kk])
    else:
        print('Exiting')
        exit()
    
if __name__ == "__main__":
    args = parser.parse_args()
    main(args.z_bin, tag=args.tag)


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
