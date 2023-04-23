import numpy as np
import glob

from fileLoc import FileLocs
floc = FileLocs(machine='nersc')
path = floc.halo_run_loc
fname_base = path+'nbody_output/tmp/%s_mock_z%i'

tag, zbin = "debug", 0



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

def joinFiles(files):
    r = np.load(files[0])['r']
    vec = [np.load(file)['mean'] for file in files]
    return r, np.array(vec)

def main(tag, zbin):
    fname = path+"'nbody_output/%s_z%i_2PCF"%(tag, zbin)

    outfiles = glob.glob(fname_base%('xi', z_bin)+"_k*")
    n_patches = len(outfiles)
    r, vec = joinFiles(outfiles)
    cov, sig, mean = get_jk_cov(vec, n_patches)

    print('Saved Output: %s'%(fname))
    np.savez(fname, r=r, mean=mean, sig=sig, cov=cov)