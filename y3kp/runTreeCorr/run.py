""" Computes the clustering of clusters

Estimates the two-point correlation function of the redMaPPer Y3 clusters using treecorr. 
The covariance is computed by jackknife of Npatches defined by kmeans.

Args:
    tag (str): tag label to the output files
    n_patches (int): number of patches used in the JK computation
    parallel (bool): to run in parallel using joblib
    nCores (int): number of cores to run

Outputs:
     fname = ./output/kmean{Npatches}/{tag}_l%i-z%i.npz
  
Example:
    vec = np.load(fname)
    list(vec.keys()) # r, mean, sig, cov
"""
import os
import numpy as np
import treecorr

import argparse
from joblib import Parallel, delayed

# local libraries
from set_config import config, config3
from util import load_mock

parser = argparse.ArgumentParser(description='Run TreeCorr in Npatches and estimate a JK covariance')
parser.add_argument('tag', type=str, help='a label to the output files')
# parser.add_argument('--dataset', default='mock', type=str, help='input dataset, available options are mock or desy3')
parser.add_argument('--nPatches', default=20, type=int, help='number of patches used by the JK cov estimator')
parser.add_argument('--is_3d', default=0, type=int, help='computes angular or projected distance correlation function')
parser.add_argument('--nCores', default=2, type=int, help='number of cores')
parser.add_argument('-p', '--parallel', default=1, type=int, help='parallel true (1) or false (0)')


# data_dict = {'mock': mock_fname_base, 'desy3': des_fname_base}

def get_angular_correlation_jk_cov(cat, rcat, n_patches=10, config=config):
    # process data catalog
    dd = treecorr.NNCorrelation(config, var_method='jackknife')
    dd.process(cat)
    
    # process random catalog (it takes longer)
    rr = treecorr.NNCorrelation(config)
    rr.process(rcat)
    
    # process data x random
    dr = treecorr.NNCorrelation(config)
    dr.process(cat, rcat)
    
    # return Landy & Szalay estimator
    xi, varxi = dd.calculateXi(rr=rr, dr=dr)
    
    cov = dd.cov
    r   = np.exp(dd.meanlogr)
    sig = np.sqrt(varxi)
    
    # number of objects
    nobjs = cat.ra.size
    
    # clean memory
    rcat = cat = 0
    rr = dd = 0
    return r, xi, cov, sig, nobjs

def save_jk_covariance(fname, output):
    r, mean, sig, cov, nobj = output
    np.savez(fname, r=r, mean=mean, cov=cov, sig=sig, nobj=nobj)

def pick_config(is_3d):
    if is_3d:
        return config3
    else:
        return config
    
def run_treecorr_jk(outfile, cat, rcat, n_patches=10, is_3d=False):
    # at least 50 clusters in each patch
    #n_patchs_min = int(np.where(cat.ra.size/n_patches<50, cat.ra.size/50, n_patches))
    
    # pick config file
    _config = pick_config(is_3d)
    
    # run treecorr on n_patches
    r, mean, sig, cov, nobj = get_angular_correlation_jk_cov(cat, rcat, config=_config,n_patches=n_patches)
    
    # save results
    save_jk_covariance(outfile, [r, mean, cov, sig, nobj])
    
############## START CODE ################
def main(tag, n_patches, is_3d=False, parallel=False, nCores=10): 
    nbins = 21
    
    print('Load Data')    
    # load data in redshift and lambda bins
    data, outfiles = load_mock(nsize=nbins, n_patches=n_patches, is_3d=is_3d)
    # print('\n'.join(sorted(outfiles)))

    print('Running trecorr on lambda/redshift bins')
    if parallel:
        Parallel(n_jobs=nCores)(delayed(run_treecorr_jk)
                            (outfiles[ii], data[ii][0], data[ii][1], n_patches, is_3d=is_3d)
                            for ii in range(nbins))
    else:
        for ii in range(nbins):
            run_treecorr_jk(outfiles[ii], data[ii][0], data[ii][1], n_patches, is_3d=is_3d)
    
    print('')
    print('Saved outfiles')
    print('\n'.join(sorted(outfiles)))
    
if __name__ == "__main__":
    args = parser.parse_args()
    is_3d = bool(args.is_3d)
        
    main(args.tag, args.nPatches, is_3d=is_3d, 
         parallel=args.parallel, nCores=args.nCores)
