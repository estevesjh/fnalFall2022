""" Computes the clustering of clusters

Estimates the two-point correlation function of the redMaPPer Y3 clusters using treecorr. 
The covariance is computed by jackknife of Npatches defined by kmeans.

Args:
    tag (str): tag label to the output files
    n_patches (int): number of patches used in the JK computation
    parallel (bool): to run in parallel using joblib
    nCores (int): number of cores to run

Outputs:
     fname = ./output/{tag}_kmean{Npatches}.npz
  
Example:
    vec = np.load(fname)
    list(vec.keys()) # r, mean, sig, cov
"""
import os
import numpy as np

from jackEstimatorKmeans import JackKniferKmeans
from util import load_data, SetupFiles, config

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

def run_treecorr_jk(jk, rm, ran, kk, fname_base):
    fname = fname_base.format(kk)
    
    # random group labels
    ran_labels = jk.add_randoms(ran['ra'], ran['dec'])

    # find indices
    idx = jk.get_mask(kk)
    idx1= jk.get_mask(kk, ran_labels)

    # run correlation functions
    _r, _xi = get_angular_correlation(rm[idx], ran[idx1])

    # append results
    np.save(fname, _xi)
    print('Saved results: %s'%fname)
    print('')

def set_jk_estimator(rm, Nk, fname):
    jk =  JackKniferKmeans(rm['ra'], rm['dec'], Nk, fname=fname)
    jk.show_stats()
    jk.write(fname)
    return jk

def get_jk_cov(jk_stats,npatches):
    """
    Get JackKnife covariance and error bar. 
 
    Args:
        jk_stats (array): array with the statistics computed for each JK patch
        npatches (float): number of JK patches 
    """
    mean  = np.mean(jk_stats, axis =0)
    cov_jk  = np.sum((jk_stats - mean[None, :])[:, :, None] *
                 (jk_stats - mean[None, :])[:, None, :], axis=0)
    cov_jk *= ((npatches - 1.)/npatches)  
    sig_jk  = np.sqrt(np.diag(cov_jk))
    return cov_jk, sig_jk, mean

def joinFiles(fname_base, Nk):
    vec = [np.load(fname_base.format(kk)) for kk in range(Nk)]
    return np.array(vec)
    
def save_jk_covariance(fname, output):
    r, cov, sig, mean = output
    np.savez(fname, r=r, mean=mean, cov=cov, sig=sig)


############## START CODE ################
def main(tag, n_patches, parallel=False, nCores=2):
    # setup output filenames
    fl = SetupFiles(tag, n_patches)
    
    print('Load Data')
    rm_all, ran_all = load_data()
    
    print('Apply lambda and redshift bin cut')
    cut = np.random.randint(len(ran_all),size=100000)
    rm, ran = rm_all, ran_all[cut]
    
    print('JK Patches')
    jk =  set_jk_estimator(rm, n_patches, fl.fname_kmeans_centers)
    
    print('Running trecorr on patches')
    if parallel:
        Parallel(n_jobs=nCores)(delayed(run_treecorr_jk)
                            (jk, rm, ran, kk, fl.tmp_fname_base)
                            for kk in range(n_patches))
    else:
        for kk in range(n_patches):
            run_treecorr_jk(jk, rm, ran, kk, fl.tmp_fname_base)
    
    print('Compute JK Covariance')
    vec = joinFiles(fl.tmp_fname_base, n_patches)
    cov, sig, mean = get_jk_cov(vec, vec.shape[0])
    
    # setup radial bin
    rbins = np.logspace(np.log10(config['min_sep']),np.log10(config['max_sep']), config['nbins'])
    
    save_jk_covariance(fl.outfile, [rbins, cov, sig, mean])
    print('Saved outfile: %s'%fl.outfile)

if __name__ == '__main__':
    tag = 'test_all'
    Npatches = 10
    main(tag, Npatches)


# To Do
# - Make ArgParser
# - Define Lambda and Redshift Bin
# - Study the n_patches
# nohup python run.py >log 2>&1 &