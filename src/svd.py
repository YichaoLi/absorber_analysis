#! /usr/bin/env python 

import os
import copy
import gc
import warnings


import numpy as np
import numpy.ma as ma
import scipy as sp
import find_modes
import matplotlib.pyplot as plt
import h5py
from numpy import linalg

from kiyopy import parse_ini, utils
import core.fitsGBT as fitsGBT
import kiyopy.custom_exceptions as ce
from time_stream import base_single
from time_stream import cal_scale
from time_stream import rotate_pol

from mpi4py import MPI

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
                   'save_plot': True,
                   'modes' : 0,
                   }
    feedback_title = 'New flags each data Block: '

    def mpiexecute(self, n_processes=1):
        """
        Process all data with MPI
        To start with MPI, you need to change manager.py 
        calling mpiexecute instead of execute. 
        and do,

        $ mpirun -np 9 --bynode  python  manager.py pipeline.pipe

        """

        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()

        params = self.params
        if rank == 0:
            output_fname = (params['output_root']
                    + params['file_middles'][0] + params['output_end'])
            utils.mkparents(output_fname)

        comm.barrier()

        n_files = len(params['file_middles'])
        for file_ind in range(n_files)[rank::size]:

            self.process_file(file_ind)

        comm.barrier()

    
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

        #Data.data = svd_spec_time(Data.data, params, self.file_ind, freq, time)
        Data.data = corr_svd(Data.data, params, self.file_ind, freq, time)

        return Data

def svd_spec_time(data, params, file_ind, freq=None, time=None):

    #data[...,:100] = np.inf
    #data[...,-100:] = np.inf
    #data[...,1640:1740] = np.inf
    #data[...,2066:2166] = np.inf

    freq_mask = np.any(np.isfinite(data), axis=(0, 2))
    weights = np.ones(data.shape)
    data_mask = np.logical_not(np.isfinite(data))
    weights[data_mask] = 0.
    data[data_mask] = 0.

    #if np.sum(weights) < np.prod(weights.shape) * 0.1:
    #    #print "Warning: too much data masked, no svd performed"
    #    msg = ("WARNING: too much data masked, no svd performed")
    #    warnings.warn(msg)
    #    data[data_mask] = np.inf
    #    return data

    sh = data.shape

    # for XX
    data_svd = data[:,0,:,:].reshape([-1, sh[-1]])[:,freq_mask[0, :]]
    weight_svd = weights[:,0,:,:].reshape([-1, sh[-1]])[:,freq_mask[0, :]]
    # check flag percent
    weight_svd = np.ma.array(weight_svd)
    weight_svd[weight_svd==0] = np.ma.masked
    percent = float(np.ma.count_masked(weight_svd)) / weight_svd.size * 100
    print 'Flag percent XX: %f%%'%percent
    if np.sum(weight_svd) < np.prod(weight_svd.shape) * 0.1 or data_svd.shape[-1]<10:
        #print "Warning: too much data masked, no svd performed"
        msg = ("WARNING: too much data masked for XX, no svd performed")
        warnings.warn(msg)
        data[data_mask] = np.inf
        return data
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
        del amp

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

    if params['save_plot']:
        f_name = params['output_root'] + \
                params['file_middles'][file_ind] + '_svd_XX.hdf5'
        utils.mkparents(f_name)
        check_svd(f_name, [val, vec_t.T, vec_f.T], freq_mask[0,:], freq)
        check_map(f_name, np.ma.array(data[:,0,:,:]), time, freq)

    del data_svd, weight_svd, val, vec_t, vec_f, amps
    gc.collect()

    # for YY
    data_svd = data[:,3,:,:].reshape([-1, sh[-1]])[:,freq_mask[3, :]]
    weight_svd = weights[:,3,:,:].reshape([-1, sh[-1]])[:,freq_mask[3, :]]
    # check flag percent
    weight_svd = np.ma.array(weight_svd)
    weight_svd[weight_svd==0] = np.ma.masked
    percent = float(np.ma.count_masked(weight_svd)) / weight_svd.size * 100
    print 'Flag percent XX: %f%%'%percent
    if np.sum(weight_svd) < np.prod(weight_svd.shape) * 0.1 or data_svd.shape[-1]<10:
        #print "Warning: too much data masked, no svd performed"
        msg = ("WARNING: too much data masked for YY, no svd performed")
        warnings.warn(msg)
        data[data_mask] = np.inf
        return data
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
        del amp

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

    if params['save_plot']:
        f_name = params['output_root'] + \
                params['file_middles'][file_ind] + '_svd_YY.hdf5'
        utils.mkparents(f_name)
        check_svd(f_name, [val, vec_t.T, vec_f.T], freq_mask[0,:], freq)
        check_map(f_name, np.ma.array(data[:,0,:,:]), time, freq)

    del data_svd, weight_svd, val, vec_t, vec_f, amps
    gc.collect()

    data[data_mask] = np.inf

    return data


def corr_svd(data, params, file_ind, freq=None, time=None):

    # get the svd mode
    #Data.data = ma.array(Data.data)
    #Data.data[np.logical_not(np.isfinite(Data.data))] = np.ma.masked

    data[..., :100] = np.inf
    data[...,-100:] = np.inf
    data[...,1640:1740] = np.inf
    data[...,2066:2166] = np.inf

    freq_mask = np.any(np.isfinite(data), axis=(0, 2))
    print data[freq_mask[None, :, None, :]].shape
    time_mask = np.all(np.isfinite(data[...,freq_mask[0,:]]), axis=(2, 3))
    weights = np.ones(data.shape)
    data_mask = np.logical_not(np.isfinite(data))
    weights[data_mask] = 0.
    #data[data_mask] = 0.

    # for XX
    if params['save_svd']:
        map1_raw = copy.deepcopy(data[:,0,0,:].T)
        map2_raw = copy.deepcopy(data[:,0,1,:].T)

    print float(np.sum(weights[:,0,:,:])) / np.prod(weights[:,0,:,:].shape) * 100
    print float(np.sum(weights[:,0,:,:][...,freq_mask[0,:]])) /\
            np.prod(weights[:,0,:,:][...,freq_mask[0,:]].shape) * 100
    print float(np.sum(weights[:,0,:,:][...,freq_mask[0,:]][time_mask[0,:],...])) /\
            np.prod(weights[:,0,:,:][...,freq_mask[0,:]][time_mask[0,:],...].shape) * 100

    if params['save_plot']:
        f_name = params['output_root'] + \
                params['file_middles'][file_ind] + '_raw_XX.hdf5'
        utils.mkparents(f_name)
        check_map(f_name, np.ma.array(data[:,0,:,:]), time, freq)
        f_name = params['output_root'] + \
                params['file_middles'][file_ind] + '_raw_YY.hdf5'
        utils.mkparents(f_name)
        check_map(f_name, np.ma.array(data[:,3,:,:]), time, freq)


    corr, weight = find_modes.freq_covariance(
            data[:,0,0,:].T,    data[:,0,1,:].T, 
            weights[:,0,0,:].T, weights[:,0,1,:].T, 
            freq_mask[0,:], freq_mask[0, :], no_weight=False)

    svd_result = find_modes.get_freq_svd_modes(corr, corr.shape[0])

    map1, map2, outmap_left, outmap_right = subtract_foregrounds(
            svd_result, 
            data[:,0,0,:].T, data[:,0,1,:].T, 
            weights[:,0,0,:].T, weights[:,0,1,:].T, 
            freq_mask[0,:], 
            0, params['modes'])

    data[:,0,0,:] = map1.T
    data[:,0,1,:] = map2.T

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

    if params['save_plot']:
        f_name = params['output_root'] + \
                params['file_middles'][file_ind] + '_svd_XX.hdf5'
        utils.mkparents(f_name)
        check_svd(f_name, svd_result, freq_mask[0,:], freq)
        check_map(f_name, np.ma.array(data[:,0,:,:]), time, freq)

    del corr, weight, svd_result, map1, map2, outmap_left, outmap_right 
    gc.collect()

    # for YY 
    if params['save_svd']:
        map1_raw = copy.deepcopy(data[:,3,0,:].T)
        map2_raw = copy.deepcopy(data[:,3,1,:].T)

    print float(np.sum(weights[:,3,:,:])) / np.prod(weights[:,3,:,:].shape) * 100

    corr, weight = find_modes.freq_covariance(
            data[:,3,0,:].T, data[:,3,1,:].T, 
            weights[:,3,0,:].T, weights[:,3,1,:].T, 
            freq_mask[3,:], freq_mask[3, :], no_weight=False)

    svd_result = find_modes.get_freq_svd_modes(corr, corr.shape[0])

    map1, map2, outmap_left, outmap_right = subtract_foregrounds(
            svd_result, 
            data[:,3,0,:].T, data[:,3,1,:].T, 
            weights[:,3,0,:].T, weights[:,3,1,:].T, 
            freq_mask[3,:], 
            0, params['modes'])

    data[:,3,0,:] = map1.T
    data[:,3,1,:] = map2.T

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

    if params['save_plot']:
        f_name = params['output_root'] + \
                params['file_middles'][file_ind] + '_svd_YY.hdf5'
        utils.mkparents(f_name)
        check_svd(f_name, svd_result, freq_mask[3,:], freq)
        check_map(f_name, np.ma.array(data[:,3,:,:]), time, freq)

    del corr, weight, svd_result, map1, map2, outmap_left, outmap_right
    gc.collect()


    #data[data_mask] = np.inf
    data = np.ma.array(data)
    data[data_mask] = np.ma.masked

    return data

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

def check_spec_fits(fits_filename):

    data_block = [fitsGBT.Reader(fits_filename).read(),]

    name = fits_filename.split('/')[-2] + ' ' + \
            fits_filename.split('/')[-1].split('.')[0]

    for Data in data_block:

        Data.calc_time()
        time = Data.time

        time0 = time.min()
        time -= time0

        freq = (np.arange(Data.data.shape[-1]) - Data.field['CRPIX1'] + 1 )
        freq *= Data.field['CDELT1']
        freq += Data.field['CRVAL1']

        map_left = Data.data[:,0,0,:]
        map_right = Data.data[:,0,1,:]

        fig = plt.figure(figsize=(10, 5))
        ax1 = fig.add_axes([0.1, 0.52, 0.80, 0.38])
        ax2 = fig.add_axes([0.1, 0.10, 0.80, 0.38])

        map_left = np.ma.array(map_left)
        map_left[np.logical_not(np.isfinite(map_left))] = np.ma.masked

        map_right = np.ma.array(map_right)
        map_right[np.logical_not(np.isfinite(map_right))] = np.ma.masked

        ax1.plot(freq, np.ma.mean(map_left, axis=0), 'k.-')
        ax2.plot(freq, np.ma.mean(map_right, axis=0), 'k.-')

        ax1.set_title(name)
        ax1.set_xlim(xmin=freq.min(), xmax=freq.max())
        #ax1.set_ylim(ymin=freq.min(), ymax=freq.max())
        ax1.minorticks_on()
        ax1.tick_params(length=4, width=1, direction='out')
        ax1.tick_params(which='minor', length=2, width=1, direction='out')
        #ax1.set_xlabel('Frequency [MHz]')
        ax1.set_ylabel('Cal on', 
                horizontalalignment='right', 
                verticalalignment='center',
                multialignment='center')
        ax1.set_xticklabels([])

        ax2.set_xlim(xmin=freq.min(), xmax=freq.max())
        #ax2.set_ylim(ymin=freq.min(), ymax=freq.max())
        ax2.minorticks_on()
        ax2.tick_params(length=4, width=1, direction='out')
        ax2.tick_params(which='minor', length=2, width=1, direction='out')
        ax2.set_xlabel('Frequency [MHz]')
        ax2.set_ylabel('Cal off',
                horizontalalignment='right', 
                verticalalignment='center',
                multialignment='center')

        plt.show()

def check_map_fits(fits_filename):

    data_block = [fitsGBT.Reader(fits_filename).read(),]

    name = fits_filename.split('/')[-2] + ' ' + \
            fits_filename.split('/')[-1].split('.')[0]

    for Data in data_block:

        Data.calc_time()
        time = Data.time

        time0 = time.min()
        time -= time0

        freq = (np.arange(Data.data.shape[-1]) - Data.field['CRPIX1'] + 1 )
        freq *= Data.field['CDELT1']
        freq += Data.field['CRVAL1']

        map_left = Data.data[:,0,0,:]
        map_right = Data.data[:,0,1,:]

        fig = plt.figure(figsize=(10, 5))
        ax1 = fig.add_axes([0.1, 0.52, 0.80, 0.38])
        ax2 = fig.add_axes([0.1, 0.10, 0.80, 0.38])
        cax1 = fig.add_axes([0.91, 0.52, 0.01, 0.38])
        cax2 = fig.add_axes([0.91, 0.10, 0.01, 0.38])

        map_left = np.ma.array(map_left)
        map_left[np.logical_not(np.isfinite(map_left))] = np.ma.masked

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

        plt.show()

def check_map(svd_filename, map=None, time=None, freq=None):

    if map == None:
        svd_file = h5py.File(svd_filename, 'r')
        time = svd_file['time'].value
        freq = svd_file['freq'].value / 1.e6
        map = np.ma.array(svd_file['map_right'].value)

    name = svd_filename.split('/')[-2] + ' ' + \
            svd_filename.split('/')[-1].split('.')[0]

    time0 = time.min()
    time -= time0

    #map_left  = np.ma.array(svd_file['map_left'].value)
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
    #vmax = mean + sigma
    #vmin = mean - sigma
    vmax = None
    vmin = None
    im1 = ax1.pcolormesh(X, Y, map_left, vmax=vmax, vmin=vmin)

    sigma = np.ma.std(map_right)
    mean = np.ma.mean(map_right)
    #vmax = mean + sigma
    #vmin = mean - sigma
    vmax = None
    vmin = None
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
    plt.close()

    #outmap_left = svd_file['outmap_left']
    #print outmap_left.shape
    #for i in range(outmap_left.shape[0]):
    #    plt.plot(np.arange(outmap_left.shape[1]), outmap_left[i, :], '.-')


    #plt.show()


def check_svd(svd_filename, svd_result=None, freq_mask=None, freq=None):

    if svd_result == None:
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
    plt.close()
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
    #output_path = '/project/ycli/data/gbt/map_making/svd_05/GBT14B_339/'
    #data_name = '78_wigglez1hr_centre_ralongmap_80'
    #output_path = '/home/p/pen/ycli/workspace/output/map_result_gbt_absorb/svd_03/GBT13B_352/'
    #data_name = '01_wigglez1hr_centre_ralongmap_36_03'

    #output_path = '/project/ycli/data/gbt/map_making/svd_05/GBT13B_352/'
    #output_path = '/project/ycli/data/gbt/map_making/flagged/GBT13B_352/'
    #data_name = '01_wigglez1hr_centre_ralongmap_37'

    output_path = '/project/ycli/data/gbt/map_making/svd_01/GBT10B_036/'
    data_name = '87_wigglez1hr_centre_ralongmap_90-99_00'

    #check_svd(output_path + data_name + '_svd_XX.hdf5')
    #check_svd(output_path + data_name + '_svd_YY.hdf5')
    #check_map(output_path + data_name + '_svd_XX.hdf5')
    #check_map(output_path + data_name + '_svd_YY.hdf5')
    #check_map_fits(output_path + data_name + '.fits')
    check_spec_fits(output_path + data_name + '.fits')

    #output_path = '/project/ycli/data/gbt/map_making/cal_rm/GBT14B_339/'
    #data_name = '78_wigglez1hr_centre_ralongmap_80_data.npy'
    #data = np.load( output_path + data_name )
    #svd_spec_time(data, None, 0)




