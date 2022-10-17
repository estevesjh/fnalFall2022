
# coding: utf-8

# # Create Input File
# ----
# 
# <h4> Description </h4>
# Create new input files to run the random forest algorithm.

import numpy as np
import os
from astropy.table import Table, join
from astropy.io.fits import getdata

def create_newVector(data,col='MI'):
    train = data['Train']
    xall = np.c_[data['xvec'],data['z_true']].T
    yall = (np.c_[data['smass'],data[col]]).T
    return yall[:,train].T,yall[:,~train].T,xall[:,train].T,xall[:,~train].T

path = '/data/des61.a/data/johnny/COSMOS/fnal2022/'
fname = path+'desCosmosML_sample.fits'

if not os.path.isfile(fname):
    names = ['../data/%s_indices_matched'%field for field in ['des','cosmos']]
    rowsDES    = np.load(names[0]+'_01arcsec.npy')
    rowsCOSMOS = np.load(names[1]+'_01arcsec.npy')

    des_deep_field_infile = '/data/des61.a/data/johnny/COSMOS/y3_deep_fields.fits'
    des0 = Table(getdata(des_deep_field_infile))

    cosmo_infile = '/data/des61.a/data/johnny/COSMOS/COSMOS2015_Laigle+_v1.1.fits'
    cosmo0  = Table(getdata(cosmo_infile))

    des   = des0[rowsDES]
    cosmo = cosmo0[rowsCOSMOS]

    des['z_true'] = cosmo['PHOTOZ']
    des['smass']  = cosmo['MASS_BEST']
    des['GALAXY'] = cosmo['TYPE']

    for li in ['I','R','Z']:
        des['M%s'%li] = cosmo['M%s'%li]

    des0    = 0
    cosmos0 = 0

    mask = (des['z_true'] <= 10.)&(des['z_true']>=0.)
    mask &= (des['smass'] <= 13.)&(des['smass']>=0.)

    desM = des[mask].copy()
    
    # load old training
    x_train = np.load('/data/des61.a/data/johnny/DESY3/galpro_files/DES_WF/data/x_train.npy')
    y_train = np.load('/data/des61.a/data/johnny/DESY3/galpro_files/DES_WF/data/y_train.npy')

    x_test = np.load('/data/des61.a/data/johnny/DESY3/galpro_files/DES_WF/data/x_test.npy')
    y_test = np.load('/data/des61.a/data/johnny/DESY3/galpro_files/DES_WF/data/y_test.npy')

    # create an astropy table
    testTrain = Table()
    yall = np.vstack([y_train,y_test])
    xall = np.vstack([x_train,x_test])
    nsize = yall.shape[0]

    testTrain['row'] = np.arange(nsize)
    testTrain['Train'] = False
    testTrain['Train'][np.arange(y_train.shape[0])] = True
    testTrain['z_true'] = yall[:,0]
    testTrain['smass'] = yall[:,1]
    testTrain['xvec'] = xall

    print('Sample size: %i'%nsize)
    print('Training size: %i'%y_train.shape[0])

    joined = join(desM,testTrain,keys=['z_true','smass'])
    joined = joined[np.unique(joined['row'],return_index=True)[1]]
    print('New sample size: %i'%len(joined))
    print('New Training size: %i'%np.count_nonzero(joined['Train']))

    joined.write(fname)
else:
    joined = Table(getdata(fname))

print('Saving Vectors')
mylist = create_newVector(joined,col='MI')

for vec,label in zip(mylist,['y_train','y_test','x_train','x_test']):
    np.save(path+label,vec)
