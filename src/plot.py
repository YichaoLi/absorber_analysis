import os
import copy
import gc
import warnings

import numpy as np
import numpy.ma as ma
import scipy as sp

import core.fitsGBT as fitsGBT

import matplotlib.pyplot as plt
import h5py

def check_spec_fits(fits_filename):

    data_block = fitsGBT.Reader(fits_filename).read()
    if not hasattr(data_block, 'index'):
        data_block = (data_block,)

    name = fits_filename.split('/')[-2] + ' ' + \
            fits_filename.split('/')[-1].split('.')[0]

    for Data in data_block[-1:]:

        Data.calc_time()
        time = Data.time

        time0 = time.min()
        time -= time0

        freq = (np.arange(Data.data.shape[-1]) - Data.field['CRPIX1'] + 1 )
        freq *= Data.field['CDELT1']
        freq += Data.field['CRVAL1']

        #Data.data[Data.data.mask] = np.inf
        map_left = Data.data[:,0,0,:]
        map_right = Data.data[:,0,1,:]

        fig = plt.figure(figsize=(10, 5))
        ax1 = fig.add_axes([0.1, 0.52, 0.80, 0.38])
        ax2 = fig.add_axes([0.1, 0.10, 0.80, 0.38])

        #map_left = np.ma.array(map_left)
        #map_left[np.logical_not(np.isfinite(map_left))] = np.ma.masked

        #map_right = np.ma.array(map_right)
        #map_right[np.logical_not(np.isfinite(map_right))] = np.ma.masked

        ax1.plot(freq, np.ma.mean(map_left[:100], axis=0), 'r-')
        ax2.plot(freq, np.ma.mean(map_right[:100], axis=0), 'r-')
        #ax1.plot(freq, np.ma.var(map_left[:100], axis=0), 'g.-')
        #ax2.plot(freq, np.ma.var(map_right[:100], axis=0), 'g.-')
        #ax1.plot(time, np.ma.mean(map_left, axis=-1), 'r-')
        #ax2.plot(time, np.ma.mean(map_right, axis=-1), 'r-')
        #ax1.plot(time, np.ma.var(map_left, axis=-1), 'g.-')
        #ax2.plot(time, np.ma.var(map_right, axis=-1), 'g.-')

        ax1.set_title(name)
        #ax1.set_xlim(xmin=time.min(), xmax=time.max())
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

        #ax2.set_xlim(xmin=time.min(), xmax=time.max())
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

    from time_stream import cal_scale
    from time_stream import rotate_pol

    data_block = fitsGBT.Reader(fits_filename).read()
    if not hasattr(data_block, 'index'):
        data_block = (data_block,)


    name = fits_filename.split('/')[-2] + ' ' + \
            fits_filename.split('/')[-1].split('.')[0]

    #for Data in data_block:
    for Data in data_block[-1:]:

        #rotate_pol.rotate(Data, (1,2,3,4))
        #rotate_pol.rotate(Data, (-5, -7, -8, -6))
        #cal_scale.scale_by_cal(Data, True, False, False, False, True)

        Data.calc_time()
        time = Data.time

        time0 = time.min()
        time -= time0

        freq = (np.arange(Data.data.shape[-1]) - Data.field['CRPIX1'] + 1 )
        freq *= Data.field['CDELT1']
        freq += Data.field['CRVAL1']

        Y, X = np.meshgrid(freq, time)

        fig = plt.figure(figsize=(10, 5))
        ax1 = fig.add_axes([0.1, 0.52, 0.80, 0.38])
        ax2 = fig.add_axes([0.1, 0.10, 0.80, 0.38])
        cax1 = fig.add_axes([0.91, 0.52, 0.01, 0.38])
        cax2 = fig.add_axes([0.91, 0.10, 0.01, 0.38])

        #Data.data[Data.data.mask] = np.inf
        #freq_mask = np.any(np.isfinite(Data.data), axis=(0, 2)) 
        #time_mask = np.all(np.isfinite(Data.data[:,0,0,:][...,freq_mask[0,:]]), axis=-1)

        #time_mask = np.all(Data.data.mask, axis=(1, 2, 3))
        #freq_mask = np.any(Data.data.mask[np.logical_not(time_mask),...], axis=0)
        #Data.data.mask[:,freq_mask] = True
        #freq_mask = np.logical_not(np.all(np.all(Data.data.mask, axis=0), axis=1))

        freq_mask = np.all(Data.data.mask, axis=0, keepdims=True)

        #ax1.plot(np.logical_not(freq_mask[0, 0, 0, :]).astype('float'), label='cal on')
        #ax1.plot(np.logical_not(freq_mask[0, 0, 1, :]).astype('float')+2, label='cal off')

        #freq_mask = np.any(np.all(Data.data.mask, axis=0), axis=1)
        #ax1.plot(np.logical_not(freq_mask[0, :]).astype('float')+4, label='cal on + off')

        #plt.legend()
        #plt.show()

        #exit()

        #print freq_mask.shape
        #freq_mask = np.logical_not(freq_mask)
        #print freq_mask.shape
        #map_right = np.zeros(Data.data.shape)[:,0,0,:]
        #map_right[:, freq_mask[0, 0, 0, :]] = 1.

        #freq_mask = np.logical_not(np.all(Data.data.mask, axis=0))
        #map_left  = np.zeros(Data.data.shape)[:,0,0,:]
        #map_left[:, freq_mask[0, 1, :]] = 1.

        map_left = Data.data[:,0,0,:]
        map_right = Data.data[:,0,1,:]
        #map_left = Data.data.mask[:,0,0,:]
        #map_right = Data.data.mask[:,0,1,:]

        #map_left[:,np.logical_not(freq_mask[0,:])] = 5.e5
        #map_left[np.logical_not(time_mask),:] = 5.e5
        #map_left = np.ma.array(map_left)
        #map_left[np.logical_not(np.isfinite(map_left))] = np.ma.masked
        #map_left[map_left.mask] = 1.e99

        #map_right = np.ma.array(map_right)
        #map_right[np.logical_not(np.isfinite(map_right))] = np.ma.masked
        #map_right[map_right.mask] = 1.e99

        sigma = np.ma.std(map_left)
        mean = np.ma.mean(map_left)
        vmax = mean + 3 * sigma
        vmin = mean - 3 * sigma
        im1 = ax1.pcolormesh(X, Y, map_left, vmax=vmax, vmin=vmin)
        #im1 = ax1.pcolormesh(X, Y, map_left)
        #im1 = ax1.pcolormesh(X, Y, map_left, vmax=1, vmin=0)

        sigma = np.ma.std(map_right)
        mean = np.ma.mean(map_right)
        vmax = mean + 3 * sigma
        vmin = mean - 3 * sigma
        im2 = ax2.pcolormesh(X, Y, map_right, vmax=vmax, vmin=vmin)
        #im2 = ax2.pcolormesh(X, Y, map_right)
        #im2 = ax2.pcolormesh(X, Y, map_right, vmax=1, vmin=0)

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

    #output_path = '/project/ycli/data/gbt/map_making/flagged_fg/GBT10B_036/'
    #output_path = '/project/ycli/data/gbt/map_making/svd_03/GBT10B_036/'
    #output_path = '/project/ycli/data/gbt/map_making/svd_03_masked/GBT10B_036/'
    output_path = '/project/ycli/data/gbt/map_making/svd_01_masked/GBT10B_036/'
    #output_path = '/project/ycli/data/gbt/map_making/raw/GBT10B_036/'
    #output_path = '/project/ycli/data/gbt/map_making/split_scans/GBT10B_036/'
    #output_path = '/project/ycli/data/gbt/map_making/flagged2/GBT10B_036/'
    data_name = '90_wigglez1hr_centre_ralongmap_27-36'

    #check_svd(output_path + data_name + '_svd_XX.hdf5')
    #check_svd(output_path + data_name + '_svd_YY.hdf5')
    #check_map(output_path + data_name + '_svd_XX.hdf5')
    #check_map(output_path + data_name + '_svd_YY.hdf5')
    check_map_fits(output_path + data_name + '.fits')
    #check_spec_fits(output_path + data_name + '.fits')


