#! /usr/bin/env python 

import os
import copy


import numpy as np
import numpy.ma as ma
import scipy as sp
import find_modes
import matplotlib.pyplot as plt
import h5py
from numpy import linalg

from kiyopy import parse_ini, utils
import core.fitsGBT
import kiyopy.custom_exceptions as ce
from time_stream import base_single
from time_stream import cal_scale
from time_stream import rotate_pol

class SVD_ForegroundClean(base_single.BaseSingle) :
    '''Pipeline module that flags RFI and other forms of bad data.

    '''

    prefix = 'svd_'
    params_init = {
                   # In multiples of the standard deviation of the whole block
                   # once normalized to the time median.
                   'perform_hanning' : False,
                   'cal_scale' : False,
                   'cal_phase' : False,
                   'cal_rm' : False,
                   # Rotate to XX,XY,YX,YY is True.
                   'rotate' : False,
                   'save_svd': True,
                   'modes' : 0,
                   }
    feedback_title = 'New flags each data Block: '
    
    def action(self, Data):
        '''Prepares Data and flags RFI.
        
        Parameters
        ----------
        Data : DataBlock
            Contains information in a usable format direct from GBT. 

        Returns
        -------
        Data : DataBlock
            The input `Data` with RFI flagged. Will also be cal scaled and
            rotated to XX,YY... if so chosen.

        '''
        params = self.params

        if params["rotate"]:
            if (tuple(Data.field['CRVAL4']) == (1, 2, 3, 4)):
                rotate_pol.rotate(Data, (-5,-7,-8,-6))
                Data.add_history('Rotated to XX,XY,YX,YY')
                print "Rotated to XX,XY,YX,YY"

        if params["cal_scale"] or params["cal_phase"]:
            cal_scale.scale_by_cal(Data, params['cal_scale'], False, False,
                                  False, rotate=params['cal_phase'])
            Data.add_history('Converted to units of noise cal temperture.')
            print 'Converted to units of noise cal temperture.'

        if params['cal_rm']:
            cal_xx = np.mean(Data.data[:, 0, 0, :] - Data.data[:, 0, 1, :], axis=0)
            cal_yy = np.mean(Data.data[:, 3, 0, :] - Data.data[:, 3, 1, :], axis=0)
            if np.median(cal_xx) > 0:
                print np.median(cal_xx) 
                Data.data[:, 0, 0, :] -= cal_xx[None, :]
                Data.data[:, 3, 0, :] -= cal_yy[None, :]
                Data.data[:, 1, 0, :] -= np.sqrt(cal_xx*cal_yy)[None, :]
                Data.data[:, 2, 0, :] -= np.sqrt(cal_xx*cal_yy)[None, :]
            else:
                Data.data[:, 0, 1, :] -= cal_xx[None, :]
                Data.data[:, 3, 1, :] -= cal_yy[None, :]
                Data.data[:, 1, 1, :] -= np.sqrt(cal_xx*cal_yy)[None, :]
                Data.data[:, 2, 1, :] -= np.sqrt(cal_xx*cal_yy)[None, :]
            Data.add_history('Remove noise cal signal.')
            print 'Remove noise cal signal.'

        #return corr_svd(Data, params, self.file_ind)

        # get freq axis and time axis
        Data.calc_time()
        time = Data.time

        freq = (np.arange(Data.data.shape[-1]) - Data.field['CRPIX1'] + 1 )
        freq *= Data.field['CDELT1']
        freq += Data.field['CRVAL1']

        Data.data = svd_spec_time(Data.data, params, self.file_ind, freq, time)

        return Data

def svd_spec_time(data, params, file_ind, freq=None, time=None):

    data[...,:100] = np.inf
    data[...,-100:] = np.inf
    data[...,1640:1740] = np.inf
    data[...,2066:2166] = np.inf
    freq_mask = np.any(np.isfinite(data), axis=(0, 2))
    weights = np.ones(data.shape)
    data_mask = np.logical_not(np.isfinite(data))
    weights[data_mask] = 0.
    data[data_mask] = 0.

    sh = data.shape

    # for XX
    data_svd = data[:,0,:,:].reshape([-1, sh[-1]])[:,freq_mask[0, :]]
    weight_svd = weights[:,0,:,:].reshape([-1, sh[-1]])[:,freq_mask[0, :]]
    vec_t, val, vec_f = linalg.svd(data_svd)
    vec_f = vec_f.T
    sorted_index = np.argsort(val)[::-1]

    vec_t = vec_t[:, sorted_index]
    vec_f = vec_f[:, sorted_index]
    val = val[sorted_index]

    modes = params['modes']
    amps = sp.empty((modes, data_svd.shape[0]))
    for i in np.arange(modes):

        amp =  sp.tensordot(vec_f[:, i], data_svd * weight_svd, axes=(0,1))
        amp /= sp.tensordot(vec_f[:, i], 
                vec_f[:, i][None, :] * weight_svd, axes=(0, 1))

        data_svd -= vec_f[:, i][None, :] * amp[:, None]

        amps[i, :] = amp

    data[:, 0, :, :][...,freq_mask[0, :]] = data_svd.reshape([sh[0], 2, -1])

    if params['save_svd']:

        f_name = params['output_root'] + \
                params['file_middles'][file_ind] + '_svd_XX.hdf5'
        utils.mkparents(f_name)
        f = h5py.File(f_name, 'w')
        f['singular_values'] = val
        f['left_vectors'] = vec_t.T
        f['right_vectors'] = vec_f.T
        #f['outmap_left'] = outmap_left
        f['outmap_right'] = amps
        #f['map_left'] = map1
        f['map_right'] = data[:,0,:,:]
        f['freq_mask'] = freq_mask[0,:]
        f['freq'] = freq
        f['time'] = time

        f.close()

    # for YY
    data_svd = data[:,3,:,:].reshape([-1, sh[-1]])[:,freq_mask[3, :]]
    weight_svd = weights[:,3,:,:].reshape([-1, sh[-1]])[:,freq_mask[3, :]]
    vec_t, val, vec_f = linalg.svd(data_svd)
    vec_f = vec_f.T
    sorted_index = np.argsort(val)[::-1]

    vec_t = vec_t[:, sorted_index]
    vec_f = vec_f[:, sorted_index]
    val = val[sorted_index]

    modes = params['modes']
    amps = sp.empty((modes, data_svd.shape[0]))
    for i in np.arange(modes):

        amp =  sp.tensordot(vec_f[:, i], data_svd * weight_svd, axes=(0,1))
        amp /= sp.tensordot(vec_f[:, i], 
                vec_f[:, i][None, :] * weight_svd, axes=(0, 1))

        data_svd -= vec_f[:, i][None, :] * amp[:, None]

        amps[i, :] = amp

    data[:, 3, :, :][...,freq_mask[3, :]] = data_svd.reshape([sh[0], 2, -1])

    if params['save_svd']:

        f_name = params['output_root'] + \
                params['file_middles'][file_ind] + '_svd_YY.hdf5'
        utils.mkparents(f_name)
        f = h5py.File(f_name, 'w')
        f['singular_values'] = val
        f['left_vectors'] = vec_t.T
        f['right_vectors'] = vec_f.T
        #f['outmap_left'] = outmap_left
        f['outmap_right'] = amps
        #f['map_left'] = map1
        f['map_right'] = data[:,3,:,:]
        f['freq_mask'] = freq_mask[3,:]
        f['freq'] = freq
        f['time'] = time

        f.close()

    data[data_mask] = np.inf

    return data


def corr_svd(Data, params, file_ind):

    # get the svd mode
    Data.data = ma.array(Data.data)
    Data.data[np.logical_not(np.isfinite(Data.data))] = np.ma.masked
    Data.data[...,1640:1740] = np.ma.masked
    Data.data[...,2066:2166] = np.ma.masked
    freq_mask = np.logical_not(np.prod(Data.data.mask, axis=(0, 2)).astype('bool'))
    weights = np.ones(Data.data.shape)
    weights[Data.data.mask]   = 0.
    Data.data[Data.data.mask] = 0.

    Data.calc_time()
    time = Data.time

    freq = (np.arange(Data.data.shape[-1]) - Data.field['CRPIX1'] + 1 )
    freq *= Data.field['CDELT1']
    freq += Data.field['CRVAL1']
    #print Data.field['CRPIX1']
    #print Data.field['CDELT1']
    #print Data.field['CRVAL1']
    #print freq / 1.e6

    # for XX
    if params['save_svd']:
        map1_raw = copy.deepcopy(Data.data[:,0,0,:].T)
        map2_raw = copy.deepcopy(Data.data[:,0,1,:].T)

    corr, weight = find_modes.freq_covariance(
            Data.data[:,0,0,:].T, Data.data[:,0,1,:].T, 
            weights[:,0,0,:].T, weights[:,0,1,:].T, 
            freq_mask[0,:], freq_mask[0, :], no_weight=False)

    svd_result = find_modes.get_freq_svd_modes(corr, corr.shape[0])

    map1, map2, outmap_left, outmap_right = subtract_foregrounds(
            svd_result, 
            Data.data[:,0,0,:].T, Data.data[:,0,1,:].T, 
            weights[:,0,0,:].T, weights[:,0,1,:].T, 
            freq_mask[0,:], 
            0, params['modes'])

    Data.data[:,0,0,:] = map1.T
    Data.data[:,0,1,:] = map2.T

    if params['save_svd']:

        f_name = params['output_root'] + \
                params['file_middles'][file_ind] + '_svd_XX.hdf5'
        utils.mkparents(f_name)
        f = h5py.File(f_name, 'w')
        f['singular_values'] = svd_result[0]
        f['left_vectors'] = svd_result[1]
        f['right_vectors'] = svd_result[2]
        f['outmap_left'] = outmap_left
        f['outmap_right'] = outmap_right
        f['map_left'] = map1
        f['map_right'] = map2
        f['raw_left'] = map1_raw
        f['raw_right'] = map2_raw
        f['freq_mask'] = freq_mask[0,:]
        f['freq'] = freq
        f['time'] = time

        f.close()

    # for YY 
    if params['save_svd']:
        map1_raw = copy.deepcopy(Data.data[:,3,0,:].T)
        map2_raw = copy.deepcopy(Data.data[:,3,1,:].T)

    corr, weight = find_modes.freq_covariance(
            Data.data[:,3,0,:].T, Data.data[:,3,1,:].T, 
            weights[:,3,0,:].T, weights[:,3,1,:].T, 
            freq_mask[3,:], freq_mask[3, :], no_weight=False)

    svd_result = find_modes.get_freq_svd_modes(corr, corr.shape[0])

    map1, map2, outmap_left, outmap_right = subtract_foregrounds(
            svd_result, 
            Data.data[:,3,0,:].T, Data.data[:,3,1,:].T, 
            weights[:,3,0,:].T, weights[:,3,1,:].T, 
            freq_mask[3,:], 
            0, params['modes'])

    Data.data[:,3,0,:] = map1.T
    Data.data[:,3,1,:] = map2.T

    if params['save_svd']:

        f_name = params['output_root'] + \
                params['file_middles'][file_ind] + '_svd_YY.hdf5'
        utils.mkparents(f_name)
        f = h5py.File(f_name, 'w')
        f['singular_values'] = svd_result[0]
        f['left_vectors'] = svd_result[1]
        f['right_vectors'] = svd_result[2]
        f['outmap_left'] = outmap_left
        f['outmap_right'] = outmap_right
        f['map_left'] = map1
        f['map_right'] = map2
        f['raw_left'] = map1_raw
        f['raw_right'] = map2_raw
        f['freq_mask'] = freq_mask[3,:]
        f['freq'] = freq
        f['time'] = time

        f.close()

        #check_svd(f_name)

    return Data

def subtract_foregrounds(svd_result, map1, map2, weight1, weight2, freq_mask,
        n_modes_start, n_modes_stop):
    r"""take the SVD modes from above and clean each LOS with them"""
    svd_modes1 = svd_result[1]
    svd_modes2 = svd_result[2]

    print "subtracting %d to %d modes " % (n_modes_start, n_modes_stop)

    svd_modes1_use = svd_modes1[n_modes_start: n_modes_stop]
    svd_modes2_use = svd_modes2[n_modes_start: n_modes_stop]

    map1, map2, outmap_left, outmap_right = \
            subtract_frequency_modes(map1, map2, weight1, weight2, freq_mask, 
            svd_modes1_use, svd_modes2_use)

    return map1, map2, outmap_left, outmap_right

def subtract_frequency_modes(map1, map2, weight1, weight2, freq_mask, 
        modes1, modes2=None):
    r"""Subtract frequency modes from the map.
    """

    if modes2 == None:
        modes2 = modes1

    # First map.
    outmap_left = sp.empty((len(modes1), ) + map1.shape[1:])

    for mode_index, mode_vector in enumerate(modes1):
        #mode_vector = mode_vector.reshape(freq_mask.shape)

        amp = sp.tensordot(mode_vector, 
                map1[freq_mask, :] * weight1[freq_mask, :], axes=(0,0))
        amp /= sp.tensordot(mode_vector, 
                mode_vector[:, None] * weight1[freq_mask, :], axes=(0,0))

        fitted = mode_vector[:, None] * amp[None, :]
        map1[freq_mask, :] -= fitted

        outmap_left[mode_index, :] = amp

    # Second map.
    outmap_right = sp.empty((len(modes2), ) + map2.shape[1:])

    for mode_index, mode_vector in enumerate(modes2):
        #mode_vector = mode_vector.reshape(freq_mask.shape)

        amp = sp.tensordot(mode_vector, 
                map2[freq_mask, :] * weight2[freq_mask, :], axes=(0,0))
        amp /= sp.tensordot(mode_vector, 
                mode_vector[:, None] * weight2[freq_mask, :], axes=(0,0))
        fitted = mode_vector[:, None] * amp[None, :]
        map2[freq_mask, :] -= fitted

        outmap_right[mode_index, :] = amp

    return map1, map2, outmap_left, outmap_right

def svd_clean(data):

    I = data.data[:, 0, :]

    freq = np.ones(I.shape[0]).astype('bool')

    corr, weight = find_modes.freq_covariance(I, I, None, None, freq, freq, 
            no_weight=True, normalize=False)

    svd_result = find_modes.get_freq_svd_modes(corr, corr.shape[0])

    return svd_result


def check_map(svd_filename):

    svd_file = h5py.File(svd_filename, 'r')
    name = svd_filename.split('/')[-2] + ' ' + \
            svd_filename.split('/')[-1].split('.')[0]

    time = svd_file['time'].value
    freq = svd_file['freq'].value / 1.e6

    time0 = time.min()
    time -= time0

    #map_left  = np.ma.array(svd_file['map_left'].value)
    map = np.ma.array(svd_file['map_right'].value)
    #map.shape = time.shape + (2,) + freq.shape
    map_left  = map[:,0,:]
    map_right = map[:,1,:]

    #map_left  = np.ma.array(svd_file['raw_left'].value)
    #map_right = np.ma.array(svd_file['raw_right'].value)

    #map_left  = map_left[:, :300]
    #map_right = map_right[:, :300]

    fig = plt.figure(figsize=(10, 5))
    ax1 = fig.add_axes([0.1, 0.52, 0.80, 0.38])
    ax2 = fig.add_axes([0.1, 0.10, 0.80, 0.38])
    cax1 = fig.add_axes([0.91, 0.52, 0.01, 0.38])
    cax2 = fig.add_axes([0.91, 0.10, 0.01, 0.38])

    map_left = np.ma.array(map_left)
    map_left[np.logical_not(np.isfinite(map_left))] = np.ma.masked
    #map_left[1640:1740, :] = np.ma.masked
    #map_left[2066:2166, :] = np.ma.masked

    Y, X = np.meshgrid(freq, time)

    map_right = np.ma.array(map_right)
    map_right[np.logical_not(np.isfinite(map_right))] = np.ma.masked

    sigma = np.ma.std(map_left)
    mean = np.ma.mean(map_left)
    vmax = mean + sigma
    vmin = mean - sigma
    im1 = ax1.pcolormesh(X, Y, map_left, vmax=vmax, vmin=vmin)

    sigma = np.ma.std(map_right)
    mean = np.ma.mean(map_right)
    vmax = mean + sigma
    vmin = mean - sigma
    im2 = ax2.pcolormesh(X, Y, map_right, vmax=vmax, vmin=vmin)

    fig.colorbar(im1, ax=ax1, cax=cax1)
    fig.colorbar(im2, ax=ax2, cax=cax2)

    ax1.set_title(name)
    ax1.set_xlim(xmin=time.min(), xmax=time.max())
    ax1.set_ylim(ymin=freq.min(), ymax=freq.max())
    ax1.minorticks_on()
    ax1.tick_params(length=4, width=1, direction='out')
    ax1.tick_params(which='minor', length=2, width=1, direction='out')
    #ax1.set_xlabel('[time]')
    ax1.set_ylabel('Frequency [MHz]\nCal on', 
            horizontalalignment='right', 
            verticalalignment='center',
            multialignment='center')
    ax1.set_xticklabels([])

    ax2.set_xlim(xmin=time.min(), xmax=time.max())
    ax2.set_ylim(ymin=freq.min(), ymax=freq.max())
    ax2.minorticks_on()
    ax2.tick_params(length=4, width=1, direction='out')
    ax2.tick_params(which='minor', length=2, width=1, direction='out')
    ax2.set_xlabel('Time + %f [s]'%time0)
    ax2.set_ylabel('Frequency [MHz]\nCal off',
            horizontalalignment='right', 
            verticalalignment='center',
            multialignment='center')

    cax1.minorticks_on()
    cax1.tick_params(length=4, width=1, direction='out')
    cax1.tick_params(which='minor', length=2, width=1, direction='out')

    cax2.minorticks_on()
    cax2.tick_params(length=4, width=1, direction='out')
    cax2.tick_params(which='minor', length=2, width=1, direction='out')

    #plt.show()
    plt.savefig(svd_filename.replace('hdf5', 'png'), format='png')

    #outmap_left = svd_file['outmap_left']
    #print outmap_left.shape
    #for i in range(outmap_left.shape[0]):
    #    plt.plot(np.arange(outmap_left.shape[1]), outmap_left[i, :], '.-')


    #plt.show()


def check_svd(svd_filename):

    svd_file = h5py.File(svd_filename, 'r')

    svd_result = [svd_file['singular_values'].value, 
                  svd_file['left_vectors'].value, 
                  svd_file['right_vectors'].value]
    freq_mask = svd_file['freq_mask'].value
    freq = svd_file['freq'].value / 1.e6

    #time = svd_file['time'].value
    #print (time - time[0])[-1]
    #print time.shape[0]
    #print (time - time[0])[-1] / np.float(time.shape[0])

    #exit()

    fig1 = plt.figure(figsize=(6,6))
    ax1 = fig1.add_axes([0.15, 0.1, 0.8, 0.85])
    ax1.plot(np.arange(svd_result[0].shape[0]), svd_result[0], '.-')
    ax1.semilogy()
    ax1.set_xlim(xmin=-1, xmax=30)
    ax1.set_ylim(ymin=1.e-1)
    ax1.minorticks_on()
    ax1.tick_params(length=4, width=1.)
    ax1.tick_params(which='minor', length=2, width=1.)
    plt.savefig(svd_filename.replace('.hdf5', '_eigval.png'))

    fig2 = plt.figure(figsize=(6,6))
    #ax2 = fig2.add_axes([0.1, 0.1, 0.85, 0.85])
    n_mode = 4
    H = 0.85
    h = H / float(n_mode)
    W = 0.80
    l = 0.15
    b = 0.1
    #freq = np.arange(freq_mask.shape[0])
    for i in range(n_mode):
        ax2 = fig2.add_axes([l, b + ( n_mode - i - 1 ) * h, W, h])
        #ax2.plot(freq[freq_mask], svd_result[1][i])
        ax2.plot(freq[freq_mask], svd_result[2][i])
        ax2.set_ylim(ymax=0.07, ymin=-0.07)
        ax2.set_xlim(xmin=freq.min(), xmax=freq.max())
        ax2.minorticks_on()
        ax2.tick_params(length=4, width=1.)
        ax2.tick_params(which='minor', length=2, width=1.)
        ax2.set_ylabel('mode %02d'%i)
        if i < n_mode-1:
            ax2.set_xticklabels([])
        else:
            ax2.set_xlabel('Frequency [MHz]')
    plt.savefig(svd_filename.replace('.hdf5', '_eigvec.png'))
    #plt.show()

    ##outmap_left = svd_file['outmap_left'].value
    #outmap_right = svd_file['outmap_right'].value
    #fig3 = plt.figure(figsize=(6,6))
    #n_mode = 4
    #H = 0.85
    #h = H / float(n_mode)
    #W = 0.80
    #l = 0.15
    #b = 0.1
    #for i in range(n_mode):
    #    ax3 = fig3.add_axes([l, b + ( n_mode - i - 1 ) * h, W, h])


    #    #ax3.plot(np.arange(outmap_left.shape[1]), outmap_left[i])
    #    ax3.plot(np.arange(outmap_right.shape[1]), outmap_right[i])

    #    #ax3.set_ylim(ymax=0.07, ymin=-0.07)
    #    #ax3.set_xlim(xmin=freq.min(), xmax=freq.max())
    #    ax3.minorticks_on()
    #    ax3.tick_params(length=4, width=1.)
    #    ax3.tick_params(which='minor', length=2, width=1.)
    #    ax3.set_ylabel('mode %02d'%i)
    #    if i < n_mode-1:
    #        ax3.set_xticklabels([])
    #    else:
    #        ax3.set_xlabel('time')
    #plt.savefig(svd_filename.replace('.hdf5', '_fg.png'))

if __name__=="__main__":

    #output_path = '/project/ycli/data/gbt/cal_removed/'
    #data_name  = 'guppi_56856_wigglez1hr_centre_0011_0001'

    #data = readGBT.load( output_path + data_name)
    #svd_result = svd_clean(data)

    #svd_file = h5py.File(output_path + data_name + '_svd.h5py')
    #svd_file['singular_values'] = svd_result[0]
    #svd_file['left_vectors'] = svd_result[1]
    #svd_file['right_vectors'] = svd_result[2]

    #svd_file.close()

    #output_path = '/home/ycli/data/gbt/map_making/svd/GBT14B_339/'
    output_path = '/project/ycli/data/gbt/map_making/svd_05/GBT14B_339/'
    data_name = '78_wigglez1hr_centre_ralongmap_80'

    check_svd(output_path + data_name + '_svd_XX.hdf5')
    check_svd(output_path + data_name + '_svd_YY.hdf5')
    check_map(output_path + data_name + '_svd_XX.hdf5')
    check_map(output_path + data_name + '_svd_YY.hdf5')

    #output_path = '/project/ycli/data/gbt/map_making/cal_rm/GBT14B_339/'
    #data_name = '78_wigglez1hr_centre_ralongmap_80_data.npy'
    #data = np.load( output_path + data_name )
    #svd_spec_time(data, None, 0)




