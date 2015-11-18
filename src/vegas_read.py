#! /bin/env python 

import numpy as np
import scipy as sp
from scipy.signal import medfilt
import pyfits
import matplotlib.pyplot as plt

import h5py

def check_fits(data_path, band_list=[]):

    for band in band_list:

        print '-- '*30
        print '  ', data_path%band
        hdulist = pyfits.open(data_path%band, memmap=False)
        tbdata = hdulist[1].data

        fieldlabel = []
        for i in range(hdulist[1].header['TFIELDS']): 
            fieldlabel.append(hdulist[1].header['TTYPE%d'%(i+1)])
        print fieldlabel
        #for i in range(len(tbdata)):
        #   print tbdata[i][fieldlabel[3]], 
        for i in range(hdulist[1].header['TFIELDS']):
            print hdulist[1].header['TTYPE%d'%(i+1)]
            if hdulist[1].header['TTYPE%d'%(i+1)] == 'DATA':
                continue
            print sp.unique(tbdata.field(fieldlabel[i])).shape
            print tbdata.field(fieldlabel[i]).shape
            print tbdata.field(fieldlabel[i])
            print
        exit()

def read_data(data_path, band_list=[], output='./test'):

    data_dict = h5py.File(output + '.hdf5', 'w')
    #data_dict = {}

    for band in band_list:

        print '-- '*30
        print '  ', data_path%band
        hdulist = pyfits.open(data_path%band, memmap=False)
        tbdata = hdulist[1].data

        #data_dict['band_%s'%band] = {}
        timestamp = tbdata.field('TIMESTAMP')

        unique_timestamp, index = sp.unique(timestamp, return_index=True)
        index = np.sort(index)
        unique_timestamp = timestamp[index]

        object = tbdata.field('OBJECT')[index]
        object = list(object)
        procscan = tbdata.field('PROCSCAN')[index]
        procscan= list(procscan)
        ra = tbdata.field('CRVAL2')[index]
        ra = list(ra)
        dec = tbdata.field('CRVAL3')[index]
        dec = list(dec)
        #print object
        #print procscan
        #print index 
        #print unique_timestamp

        data = np.array(tbdata.field('DATA'))
        expo = np.array(tbdata.field('EXPOSURE'))

        nfreq = data.shape[1]
        freq_crpix = tbdata.field('CRPIX1')[0]
        freq_delta = tbdata.field('CDELT1')[0]
        freq_crval = tbdata.field('CRVAL1')[0]
        freq = ((np.arange(nfreq) + 1) - freq_crpix) * freq_delta + freq_crval

        print
        print "\tfreq_delta: %5.4f, \n\tfreq_crval: %5.4f, \n\tfreq_N: %d"\
                %(freq_delta, freq_crval, nfreq)
        print

        data_dict['band_%s/freq'%(band,) ] = freq
        for i in range(len(index)):
            print '\t--' + unique_timestamp[i] + " RA: %4.2f, Dec %4.2f"%(ra[i], dec[i])
            print '\t  band_%s/%02d_%s_%s'%(band, i, object[i], procscan[i])
            data_dict['band_%s/%02d_%s_%s'%(band, i, object[i], procscan[i])]\
                    = data[timestamp == unique_timestamp[i], :].reshape([-1, 4, 2, nfreq])
            data_dict['band_%s/%02d_%s_%s_RA'%(band, i, object[i], procscan[i])]  = ra[i]
            data_dict['band_%s/%02d_%s_%s_DEC'%(band, i, object[i], procscan[i])] = dec[i]
            data_dict['band_%s/%02d_%s_%s_EXPOSURE'%(band, i, object[i], procscan[i])] = \
                    expo[timestamp == unique_timestamp[i]].reshape([-1, 4, 2])
            print '\t  int time: '+\
                    'XX[%4.2f, %4.2f], YY[%4.2f, %4.2f], XY[%4.2f, %4.2f], YX[%4.2f, %4.2f]'\
                    %tuple(np.sum(expo[timestamp == unique_timestamp[i]].reshape([-1, 4, 2]),
                        axis=0).flatten())

        hdulist.close()
        data_dict.flush()
        print

    data_dict.close()

def read_data_old(data_path, band_list=[]):

    data_dict = {}

    print

    for band in band_list:

        data_dict['band_%s'%band] = {}

        print data_path%band
        hdulist = pyfits.open(data_path%band)
        tbdata = hdulist[1].data

        # get the long.,lat. and timestamp
        longitude = tbdata.field('CRVAL2')
        latitude = tbdata.field('CRVAL3')
        timestamp = tbdata.field('TIMESTAMP')

        procscan = tbdata.field('PROCSCAN')

        unique_timestamp, index = sp.unique(timestamp, return_index=True)
        index = np.sort(index)
        unique_timestamp = timestamp[index]

        #object, index = sp.unique(tbdata.field('OBJECT'), return_index=True)
        #index = np.sort(index)
        object = tbdata.field('OBJECT')[index]
        object = list(object)
        print object, index, unique_timestamp

        # there are 13 different timestamps
        # 0 ~ 3 calibrator observation, [0, 2] ON source, [1, 3] OFF source
        # [4, 5] on target observation 
        # 6 nearby source observation 
        # [7, 8] off target observation
        # 9 ~ 12 calibrator observation, [9, 11] ON source, [10, 12] OFF source

        data = np.array(tbdata.field('DATA'))

        nfreq = data.shape[1]
        freq_crpix = tbdata.field('CRPIX1')[0]
        freq_delta = tbdata.field('CDELT1')[0]
        freq_crval = tbdata.field('CRVAL1')[0]
        freq = ((np.arange(nfreq) + 1) - freq_crpix) * freq_delta + freq_crval

        print "freq_delta: %5.4f, freq_crval: %5.4f, freq_N: %d"%(freq_delta, freq_crval, nfreq)
        print

        #data.shape = [-1, 4, 2, nfreq]
        #print data.shape

        data_dict['band_%s'%band]['freq'] = freq
        data_dict['band_%s'%band]['object'] = [object[0]+'_ON', object[0]+'_OFF'] + object[1:]
        data_dict['band_%s'%band][object[0] + '_ON'] = data[np.logical_or(
            np.logical_or(timestamp == unique_timestamp[0], 
                          timestamp == unique_timestamp[2]),
            np.logical_or(timestamp == unique_timestamp[9], 
                          timestamp == unique_timestamp[11])), :]
        data_dict['band_%s'%band][object[0] + '_OFF'] = data[np.logical_or(
            np.logical_or(timestamp == unique_timestamp[1],
                          timestamp == unique_timestamp[3]),
            np.logical_or(timestamp == unique_timestamp[10], 
                          timestamp == unique_timestamp[12])), :]
        #data_dict['band_%s'%band][object[1]] = data[np.logical_or(
        #    timestamp == unique_timestamp[4], timestamp == unique_timestamp[5]) , :]
        data_dict['band_%s'%band][object[1]] = data[timestamp == unique_timestamp[5], :]
        data_dict['band_%s'%band][object[2]] = data[timestamp == unique_timestamp[6], :]
        data_dict['band_%s'%band][object[3]] = data[timestamp == unique_timestamp[7], :]
        data_dict['band_%s'%band][object[4]] = data[timestamp == unique_timestamp[8], :]

        data_dict['band_%s'%band][object[0] + '_ON'].shape  = (-1, 4, 2, nfreq)
        data_dict['band_%s'%band][object[0] + '_OFF'].shape = (-1, 4, 2, nfreq)
        data_dict['band_%s'%band][object[1]].shape          = (-1, 4, 2, nfreq)
        data_dict['band_%s'%band][object[2]].shape          = (-1, 4, 2, nfreq)
        data_dict['band_%s'%band][object[3]].shape          = (-1, 4, 2, nfreq)
        data_dict['band_%s'%band][object[4]].shape          = (-1, 4, 2, nfreq)

        #plt.figure(figsize=(5,5))
        #for time in sp.unique(timestamp)[4:9]:
        #    print longitude[timestamp==time],latitude[timestamp==time]
        #    plt.plot(longitude[timestamp==time], latitude[timestamp==time], 'o', label=time)
        #plt.legend()
        #plt.show()

    return data_dict
    
    exit()

def sub_cal(source1, cal_unit=False):

    t = source1.shape[0]
    bad_time = np.any(0 * source1.reshape(t, -1), axis=1)
    #print bad_time
    source1 = source1[np.logical_not(bad_time), ...]
    #source1 = np.ma.array(source1)
    #source1[np.logical_not(np.isfinite(source1))] = np.ma.masked
    cal = source1[:,:,1,:] - source1[:,:,0,:]
    #cal[np.logical_not(np.isfinite(cal))] = np.ma.masked
    #cal_mean = np.mean(cal, axis=0)
    cal_mean = np.median(cal, axis=0)
    #cal_mean[cal_mean==0] = np.ma.masked
    #cal_mean[np.logical_not(np.isfinite(cal_mean))] = np.ma.masked
    source1[:,:,1,:] -= cal_mean[None,:,:]

    #cal_level = cal_mean[None,:,None,:] / source1
    #cal_level = np.mean(np.mean(cal_level, axis=0), axis=-1)
    #print cal_level

    if cal_unit:
        cal_mean[cal_mean==0] = np.inf
        source1 /= cal_mean[None,:,None,:]

    return source1, cal_mean

def flag_RFI(data):

    data = np.ma.array(data)

    sigma_raw = np.ma.std(data)

    while 1:
        spec = np.ma.mean(data, axis=0)
        spec_mean = np.ma.mean(spec)
        spec_std = np.ma.std(spec)
        print spec_mean, spec_std
        spec_bad = np.logical_or(spec > spec_mean + 6.* spec_std,
                spec < spec_mean - 6.* spec_std)
        data[:, spec_bad] = np.ma.masked

        sigma = np.ma.std(data)
        if sigma_raw - sigma > 0.1 * sigma_raw:
            sigma_raw = sigma
        else:
            return data

def get_spec(data_path):

    band_list = ['A', 'B', 'C']

    data_dict = read_data(data_path)

    for i in range(len(band_list)):
        band = band_list[i]

        object =  data_dict['band_%s'%band]['object']
        freq = data_dict['band_%s'%band]['freq']/1.e9

        source1, tcal1 = sub_cal(data_dict['band_%s'%band][object[2]])
        source2, tcal2 = sub_cal(data_dict['band_%s'%band][object[3]])
        ref1, tcal_ref1  = sub_cal(data_dict['band_%s'%band][object[4]])
        ref2, tcal_ref2  = sub_cal(data_dict['band_%s'%band][object[5]])

        # integrated over time and cal_on off
        source1 = np.ma.mean(np.ma.mean(source1, axis=2), axis=0)
        source2 = np.ma.mean(np.ma.mean(source2, axis=2), axis=0)
        ref1 = np.ma.mean(np.ma.mean(ref1, axis=2), axis=0)
        ref2 = np.ma.mean(np.ma.mean(ref2, axis=2), axis=0)

        source1_I = source1[0,:] + source1[1,:]
        source2_I = source2[0,:] + source2[1,:]
        ref1_I = ref1[0,:] + ref1[1,:]
        ref2_I = ref2[0,:] + ref2[1,:]

        ref1_I = medfilt(ref1_I, kernel_size=(11,))
        ref2_I = medfilt(ref2_I, kernel_size=(11,))

        print source1.shape
        print ref1.shape

        text_data = [freq[:,None], ref1.T, ref2.T,
                source1.T, (source1_I/ref1_I)[:,None], 
                source2.T, (source2_I/ref2_I)[:,None]]

        text_data = np.concatenate(text_data, axis=1)

        np.savetxt('./Absorber_Band%s.txt'%band, text_data, fmt='%.18e', 
                header='freq \t Ref1 XX \t Ref1 YY \t Ref1 XY \t Ref1 YX'\
                +'\t Ref1 XX \t Ref1 YY \t Ref1 XY \t Ref1 YX'\
                +'\t Target1 XX \t Target1 YY \t Target1 XY \t Target1 YX '\
                +'\t Target1 I / Ref I'\
                +'\t Target2 XX \t Target2 YY \t Target2 XY \t Target2 YX '\
                +'\t Target2 I / Ref I')

def plot_all(data_txt):

    data = np.loadtxt(data_txt)
    plt.step(data[:,0], data[:,13], 'r', linewidth=1.5)
    plt.step(data[:,0], data[:,18], 'b', linewidth=1.5)
    plt.show()

def analysis(data_path, plot='fine'):

    band_list = ['A', 'B', 'C']

    data_dict = read_data(data_path)

    if plot == '1d':
        fig1d = plt.figure(figsize=(8,6))
        ax1d = []
        ax1d.append(fig1d.add_axes([0.07, 0.67, 0.91, 0.27]))
        ax1d.append(fig1d.add_axes([0.07, 0.07, 0.43, 0.54]))
        ax1d.append(fig1d.add_axes([0.55, 0.07, 0.43, 0.54]))
    elif plot == '2d':
        fig = plt.figure(figsize=(10,10))
        ax  = []
        cax = []
        ax.append(fig.add_axes([0.07, 0.67, 0.81, 0.27]))
        ax.append(fig.add_axes([0.07, 0.36, 0.81, 0.27]))
        ax.append(fig.add_axes([0.07, 0.05, 0.81, 0.27]))
        cax.append(fig.add_axes([0.9, 0.67, 0.02, 0.27]))
        cax.append(fig.add_axes([0.9, 0.36, 0.02, 0.27]))
        cax.append(fig.add_axes([0.9, 0.05, 0.02, 0.27]))
    elif plot == 'fine':
        #band_list = ['B', 'B', 'B']
        band_list = ['A',]
        fig = plt.figure(figsize=(10,10))
        ax  = []
        ax.append(fig.add_axes([0.07, 0.69, 0.91, 0.27]))
        ax.append(fig.add_axes([0.07, 0.37, 0.91, 0.27]))
        ax.append(fig.add_axes([0.07, 0.05, 0.91, 0.27]))

        #freq_ref = 0.770510
        #freq_ref = 0.771240
        freq_ref = 0.762460

        ax[0].set_xlim(xmin=-500, xmax=500)
        ax[1].set_xlim(xmin=-100, xmax=100)
        ax[2].set_xlim(xmin=-10,  xmax=10)

        ax[0].set_ylim(ymin=-15, ymax=15)
        ax[1].set_ylim(ymin=-15, ymax=15)
        ax[2].set_ylim(ymin=-15, ymax=15)
        #ax[0].set_ylim(ymin=10, ymax=40)
        #ax[1].set_ylim(ymin=10, ymax=40)
        #ax[2].set_ylim(ymin=10, ymax=40)
    elif plot == 'fine2':
        band_list = ['A',]
        fig = plt.figure(figsize=(8,9))
        ax  = []
        ax.append(fig.add_axes([0.1, 0.760, 0.88, 0.225]))
        ax.append(fig.add_axes([0.1, 0.530, 0.88, 0.225]))
        ax.append(fig.add_axes([0.1, 0.300, 0.88, 0.225]))
        ax.append(fig.add_axes([0.1, 0.070, 0.88, 0.225]))

        freq_ref = 0.762460

        freq_band = 199
        ax[0].set_xlim(xmin=-freq_band, xmax=freq_band)
        ax[1].set_xlim(xmin=-freq_band, xmax=freq_band)
        ax[2].set_xlim(xmin=-freq_band, xmax=freq_band)
        ax[3].set_xlim(xmin=-freq_band, xmax=freq_band)

        #ax[0].set_ylim(ymin=-15, ymax=15)
        #ax[1].set_ylim(ymin=-15, ymax=15)
        #ax[2].set_ylim(ymin=-15, ymax=15)

        ax[0].set_ylim(ymin=21, ymax=26)
        ax[1].set_ylim(ymin=21, ymax=26)

        ax[2].set_ylim(ymin=-2, ymax=2)
        ax[3].set_ylim(ymin=-2, ymax=2)

    elif plot == 'fine3':
        band_list = ['A',]
        fig = plt.figure(figsize=(8,9))
        ax  = []
        ax.append(fig.add_axes([0.1, 0.760, 0.88, 0.225]))
        ax.append(fig.add_axes([0.1, 0.530, 0.88, 0.225]))
        ax.append(fig.add_axes([0.1, 0.300, 0.88, 0.225]))
        ax.append(fig.add_axes([0.1, 0.070, 0.88, 0.225]))

        freq_ref = 0.762460
        #freq_ref = 0.771240

        freq_band = 299
        ax[0].set_xlim(xmin=-freq_band, xmax=freq_band)
        ax[1].set_xlim(xmin=-freq_band, xmax=freq_band)
        ax[2].set_xlim(xmin=-freq_band, xmax=freq_band)
        ax[3].set_xlim(xmin=-freq_band, xmax=freq_band)

    elif plot == 'all':
        band_list = ['A',]
        fig = plt.figure(figsize=(8,16))
        ax  = []
        freq_band = 499
        #freq_ref = 0.780000
        freq_ref = 0.800000
        for i in range(10):
            ax.append(fig.add_axes([0.1, 0.05 + 0.094*i, 0.88, 0.0938]))
            ax[-1].minorticks_on()
            ax[-1].tick_params(length=6, width=1.)
            ax[-1].tick_params(which='minor', length=3, width=1.)
            ax[-1].set_xlim(xmin=-freq_band, xmax=freq_band)
            ax[-1].set_ylim(ymin=2.03, ymax=2.15)
            #ax[-1].set_ylabel('T/T_cal')
            #ax[0].set_xlabel('freq - %8.7f'%(freq_ref) + r'$\times10^{6}$ [kHz]')
            if i != 0:
                ax[-1].set_xticklabels([])
        ax = ax[::-1]

    else:
        pass

    for i in range(len(band_list)):
        band = band_list[i]

        object =  data_dict['band_%s'%band]['object']
        freq = data_dict['band_%s'%band]['freq']/1.e9

        source1, tcal1 = sub_cal(data_dict['band_%s'%band][object[2]])
        source2, tcal2 = sub_cal(data_dict['band_%s'%band][object[3]])
        ref1, tcal_ref1  = sub_cal(data_dict['band_%s'%band][object[4]])
        ref2, tcal_ref2  = sub_cal(data_dict['band_%s'%band][object[5]])

        ref = 0.5*(np.mean(np.mean(ref1, axis=2), axis=0) + 
                np.mean(np.mean(ref2, axis=2), axis=0))
        ref = medfilt(ref, kernel_size=(1, 11))

        ref1 = medfilt(np.mean(np.mean(ref1, axis=2), axis=0), kernel_size=(1, 11))
        ref2 = medfilt(np.mean(np.mean(ref2, axis=2), axis=0), kernel_size=(1, 11))

        I1_ref = np.mean(source1, axis=2) / ref1[None, :, :]
        I2_ref = np.mean(source2, axis=2) / ref2[None, :, :]
        I1 = np.mean(source1, axis=2)
        I2 = np.mean(source2, axis=2)
        #I1 = np.mean(source2, axis=2)
        #I1 = ref[None, :, :]
        I1_ref = I1_ref[:,0,:] + I1_ref[:,1,:]
        I2_ref = I2_ref[:,0,:] + I2_ref[:,1,:]

        print freq[:,None].shape
        print ref.T.shape
        print np.ma.mean(I1, axis=0).T.shape
        text_data = [freq[:,None], ref.T, 
                np.ma.mean(I1, axis=0).T, np.ma.mean(I1_ref, axis=0)[:,None], 
                np.ma.mean(I2, axis=0).T, np.ma.mean(I2_ref, axis=0)[:,None]]

        text_data = np.concatenate(text_data, axis=1)

        np.savetxt('./Absorber_Band%s.txt'%band, text_data, fmt='%.18e', 
                header='freq \t Ref XX \t Ref YY \t Ref XY \t Ref YX'\
                +'\t Target1 XX \t Target1 YY \t Target1 XY \t Target1 YX \t Target1 I / Ref I'\
                +'\t Target2 XX \t Target2 YY \t Target2 XY \t Target2 YX \t Target2 I / Ref I')

        continue

        I1 = I1[:,0,:] + I1[:,1,:]
        I2 = I2[:,0,:] + I2[:,1,:]

        #np.save('./target1_dRef_band%s'%band, np.ma.mean(I1_ref, axis=0))
        #np.save('./target2_dRef_band%s'%band, np.ma.mean(I2_ref, axis=0))

        t = np.arange(I1.shape[0])
        F, T = np.meshgrid(freq, t)

        #I1 = flag_RFI(I1)

        spec = np.ma.mean(I1_ref, axis=0)

        sigma = np.ma.std(spec)
        mean  = np.ma.mean(spec)
        vmax = mean + 5*sigma
        vmin = mean - 5*sigma
        print mean, sigma

        if plot == '2d':
            im = ax[i].pcolormesh(F, T, I1, vmax=vmax, vmin=vmin)
            ax[i].set_ylim(ymin=0, ymax=t.max())
            ax[i].set_xlim(xmin=freq.min(), xmax=freq.max())
            ax[i].set_xlabel('freq [GHz]')
            ax[i].set_ylabel('time')

            fig.colorbar(im, cax=cax[i], ax=ax[i])
        elif plot == '1d':
            nfreq = freq.shape[0]
            print nfreq
            nfreq = (nfreq / 1000) * 1000
            spec_average = np.ma.mean(spec[:nfreq].reshape(1000, -1), axis=1)
            freq_average = np.ma.mean(freq[:nfreq].reshape(1000, -1), axis=1)

            ax1d[i].plot(freq, spec, 'g-', linewidth=0.5)
            ax1d[i].plot(freq_average, spec_average, 'r-', linewidth=1.5)
            ax1d[i].set_ylim(ymin=vmin, ymax=vmax)
            ax1d[i].set_xlim(xmin=freq.min(), xmax=freq.max())
            ax1d[i].set_xlabel('freq [GHz]')
            #ax1d[i].set_ylabel('T/T_cal')
            ax1d[i].set_ylabel('T')

            ax1d[i].minorticks_on()
            ax1d[i].tick_params(length=6, width=1.)
            ax1d[i].tick_params(which='minor', length=3, width=1.)

            if i != 0:
                ax1d[0].axvspan(xmin=freq.min(), xmax=freq.max(), 
                        ymin=vmin, ymax=vmax, color='b', alpha=0.5)
                ax1d[0].axvline(x = freq.min(), ymin=vmin, ymax=vmax, color='b')
                ax1d[0].axvline(x = freq.max(), ymin=vmin, ymax=vmax, color='b')
        elif plot == 'fine':

            freq -= freq_ref
            freq *= 1.e6

            ax[i].plot(freq, spec, 'g-', linewidth=0.5)
            ax[i].minorticks_on()
            ax[i].tick_params(length=6, width=1.)
            ax[i].tick_params(which='minor', length=3, width=1.)
            ax[i].set_ylabel('T/T_cal')
            ax[i].set_xlabel('freq - %8.7f'%(freq_ref) + r'$\times10^{6}$ [kHz]')


        elif plot == 'fine2':

            freq -= freq_ref
            freq *= 1.e6

            plot_range = np.logical_and(freq > -freq_band, freq < freq_band)
            freq = freq[plot_range]

            spec1 = np.ma.mean(I1, axis=0)
            spec1 = spec1[plot_range]
            #ax[0].plot(freq, spec1, 'g-', linewidth=1.5, label='Target1')
            ax[0].step(freq, ref[0, :][plot_range] + ref[1, :][plot_range], 
                    '0.5', where='mid', linewidth=1.5, label='Ref')
            ax[0].step(freq, spec1, 'g-', where='mid', linewidth=1.5, label='Target1')
            ax[0].minorticks_on()
            ax[0].tick_params(length=6, width=1.)
            ax[0].tick_params(which='minor', length=3, width=1.)
            ax[0].set_ylabel('T/T_cal')
            #ax[0].set_xlabel('freq - %8.7f'%(freq_ref) + r'$\times10^{6}$ [kHz]')
            ax[0].set_xticklabels([])
            ax[0].legend(frameon=False, loc=4)

            #spec2 = ref[0, :] + ref[1, :]
            spec2 = np.ma.mean(I2, axis=0)
            spec2 = spec2[plot_range]
            ax[1].step(freq, ref[0, :][plot_range] + ref[1, :][plot_range], 
                    '0.5', where='mid', linewidth=1.5, label='Ref')
            ax[1].step(freq, spec2, 'g-', where='mid', linewidth=1.5, label='Target2')
            ax[1].minorticks_on()
            ax[1].tick_params(length=6, width=1.)
            ax[1].tick_params(which='minor', length=3, width=1.)
            ax[1].set_ylabel('T/T_cal')
            #ax[1].set_xlabel('freq - %8.7f'%(freq_ref) + r'$\times10^{6}$ [kHz]')
            ax[1].set_xticklabels([])
            ax[1].legend(frameon=False, loc=4)

            spec3 = np.ma.mean(I1_ref, axis=0)
            spec3 = spec3[plot_range]
            ax[2].step(freq, spec3, 'g-', where='mid', linewidth=1.5, label='Target1 - Ref')
            ax[2].minorticks_on()
            ax[2].tick_params(length=6, width=1.)
            ax[2].tick_params(which='minor', length=3, width=1.)
            ax[2].set_ylabel('T/T_cal')
            #ax[2].set_xlabel('freq - %8.7f'%(freq_ref) + r'$\times10^{6}$ [kHz]')
            ax[2].set_xticklabels([])
            ax[2].legend(frameon=False, loc=4)

            spec4 = np.ma.mean(I2_ref, axis=0)
            spec4 = spec4[plot_range]
            ax[3].step(freq, spec4, 'g-', where='mid', linewidth=1.5, label='Target2 - Ref')
            ax[3].minorticks_on()
            ax[3].tick_params(length=6, width=1.)
            ax[3].tick_params(which='minor', length=3, width=1.)
            ax[3].set_ylabel('T/T_cal')
            ax[3].set_xlabel('freq - %8.7f'%(freq_ref) + r'$\times10^{6}$ [kHz]')
            ax[3].legend(frameon=False, loc=4)

        elif plot == 'fine3':

            freq -= freq_ref
            freq *= 1.e6

            plot_range = np.logical_and(freq > -freq_band, freq < freq_band)
            freq = freq[plot_range]

            #spec_ref = ref[0, :][plot_range] + ref[1, :][plot_range]

            spec1 = np.ma.mean(I1, axis=0)
            spec1 = spec1[plot_range]
            ax[0].step(freq, spec1, 'b-', where='mid', linewidth=1.5, label='Target1')
            #ax[0].step(freq, spec1 / spec_ref, 'b-', where='mid', 
            #        linewidth=1.5, label='Target1 / Ref')
            ax[0].minorticks_on()
            ax[0].tick_params(length=6, width=1.)
            ax[0].tick_params(which='minor', length=3, width=1.)
            ax[0].set_ylabel('T')
            ax[0].set_xticklabels([])
            ax[0].legend(frameon=False, loc=1)

            spec2 = np.ma.mean(I2, axis=0)
            spec2 = spec2[plot_range]
            ax[1].step(freq, spec2, 'r-', where='mid', linewidth=1.5, label='Target2')
            #ax[1].step(freq, spec2 / spec_ref, 'r-', where='mid', 
            #        linewidth=1.5, label='Target2 / Ref')
            ax[1].minorticks_on()
            ax[1].tick_params(length=6, width=1.)
            ax[1].tick_params(which='minor', length=3, width=1.)
            ax[1].set_ylabel('T')
            ax[1].set_xticklabels([])
            ax[1].legend(frameon=False, loc=1)

            #ax[2].step(freq, ref[0, :][plot_range] + ref[1, :][plot_range], 
            #        '0.5', where='mid', linewidth=1.5, label='(Ref1 + Ref2)/2')
            spec1_ref = np.ma.mean(I1_ref, axis=0)
            spec1_ref = spec1_ref[plot_range]
            ax[2].step(freq, spec1_ref, 'b-', where='mid', 
                    linewidth=1.5, label='Target1 / Ref')
            spec2_ref = np.ma.mean(I2_ref, axis=0)
            spec2_ref = spec2_ref[plot_range]
            ax[2].step(freq, spec2_ref, 'r-', where='mid', 
                    linewidth=1.5, label='Target2 / Ref')
            ax[2].minorticks_on()
            ax[2].tick_params(length=6, width=1.)
            ax[2].tick_params(which='minor', length=3, width=1.)
            ax[2].set_ylabel('T')
            ax[2].set_xticklabels([])
            ax[2].legend(frameon=False, loc=1)

            ax[3].step(freq, tcal1[0, :][plot_range] + tcal1[1, :][plot_range], 
                    'b', where='mid', linewidth=1.5, label='Tcal Target1')
            ax[3].step(freq, tcal2[0, :][plot_range] + tcal2[1, :][plot_range], 
                    'r', where='mid', linewidth=1.5, label='Tcal Target2')
            ax[3].step(freq, tcal_ref1[0, :][plot_range] + tcal_ref1[1, :][plot_range], 
                    'k', where='mid', linewidth=1.5, label='Tcal Ref1')
            ax[3].step(freq, tcal_ref2[0, :][plot_range] + tcal_ref2[1, :][plot_range], 
                    'c', where='mid', linewidth=1.5, label='Tcal Ref2')
            ax[3].minorticks_on()
            ax[3].tick_params(length=6, width=1.)
            ax[3].tick_params(which='minor', length=3, width=1.)
            ax[3].set_ylabel('T')
            #ax[3].set_xticklabels([])
            ax[3].set_xlabel('freq - %8.7f'%(freq_ref) + r'$\times10^{6}$ [kHz]')
            ax[3].legend(frameon=False, loc=1)
        if plot == 'all':

            spec1_ref = np.ma.mean(I1_ref, axis=0)
            spec2_ref = np.ma.mean(I2_ref, axis=0)

            for i in range(10):

                freq_sub = freq - freq_ref + i * 2 * freq_band * 1.e-6
                freq_sub *= 1.e6
                plot_range = np.logical_and(freq_sub > -freq_band, freq_sub < freq_band)

                freq_sub = freq_sub[plot_range]

                spec1 = spec1_ref[plot_range]
                spec2 = spec2_ref[plot_range]

                ax[i].step(freq_sub, spec1, 'b-', where='mid', 
                        linewidth=1.5)
                ax[i].step(freq_sub, spec2, 'r-', where='mid', 
                        linewidth=1.5)
                ax[i].set_ylabel('%4.4fMHz'%(freq_ref * 1.e3 + i * 2 * freq_band * 1.e-3))

    #plt.savefig('./spec.pdf', formate='pdf')
    plt.savefig('./spec.png', formate='png')
    plt.show()

def plot_spec(data_path, file_name, freq_ref=0.762460, freq_band=299):

    data = np.loadtxt(data_path + file_name)

    freq = data[:,0]

    fig = plt.figure(figsize=(8,5))
    ax  = []
    ax.append(fig.add_axes([0.1, 0.55, 0.88, 0.4]))
    ax.append(fig.add_axes([0.1, 0.10, 0.88, 0.4]))

    #freq_ref = 0.762460
    #freq_ref = 0.771240
    #freq_band = 299

    ax[0].set_xlim(xmin=-freq_band, xmax=freq_band)
    ax[1].set_xlim(xmin=-freq_band, xmax=freq_band)

    #ax[0].set_ylim(ymin=1.016, ymax=1.03)
    #ax[1].set_ylim(ymin=1.06, ymax=1.07)

    freq -= freq_ref
    freq *= 1.e6

    plot_range = np.logical_and(freq > -freq_band, freq < freq_band)
    freq = freq[plot_range]

    data = data[plot_range,:]

    mean = np.median(data[:,13])
    std  = np.std(data[:,13])
    ax[0].set_ylim(ymin=mean - 1.5 * std, ymax=mean + 3 * std)

    ax[0].step(freq, data[:,13], 'b-', where='mid', linewidth=1.5, label='Target1/Ref1')
    ax[0].minorticks_on()
    ax[0].tick_params(length=6, width=1.)
    ax[0].tick_params(which='minor', length=3, width=1.)
    ax[0].set_ylabel('T')
    ax[0].set_xticklabels([])
    ax[0].legend(frameon=False, loc=1)

    mean = np.median(data[:,18])
    std  = np.std(data[:,18])
    ax[1].set_ylim(ymin = mean - 1.5 * std, ymax = mean + 3 * std)

    ax[1].step(freq, data[:,18], 'b-', where='mid', linewidth=1.5, label='Target2/Ref2')
    ax[1].minorticks_on()
    ax[1].tick_params(length=6, width=1.)
    ax[1].tick_params(which='minor', length=3, width=1.)
    ax[1].set_ylabel('T')
    #ax[1].set_xticklabels([])
    ax[1].set_xlabel('freq - %8.7f'%(freq_ref) + r'$\times10^{6}$ [kHz]')
    ax[1].legend(frameon=False, loc=1)

    plt.show()

if __name__=="__main__":

    output_path = './new_data/'

    #data_path = '/home/ycli/data/absorber/AGBT_15A_196/AGBT15A_196_01.raw.vegas.%s.fits'
    #output_name = 'AGBT15A_196_01'
    #band_list = ['A', 'B', 'C']

    #data_path = '/project/zuoshifan/AGBT15A_196_02/AGBT15A_196_02.raw.vegas/AGBT15A_196_02.raw.vegas.%s.fits'
    #output_name = 'AGBT15A_196_02'
    #band_list = ['A', ]  + ['C', 'E']

    data_path = '/project/zuoshifan/AGBT15A_196_03/AGBT15A_196_03.raw.vegas/AGBT15A_196_03.raw.vegas.%s.fits'
    output_name = 'AGBT15A_196_03'
    band_list = ['A', ]  + ['C', 'E']

    read_data(data_path, band_list=band_list, output = output_path + output_name)
    #check_fits(data_path, band_list=band_list)


#get_spec(data_path)
#plot_all('./AGBT_15A_196_TXT/Absorber_BandA.txt')
#read_data(data_path)
#analysis(data_path)
#analysis(data_path, plot='1d')
#analysis(data_path, plot='fine2')
#analysis(data_path, plot='fine3')
#analysis(data_path, plot='all')

#freq_ref = 0.771240
#freq_ref = 0.867773
#freq_ref = 0.762460

#plot_spec('./AGBT_15A_196_TXT/', 'Absorber_BandA.txt', freq_ref=freq_ref)
#exit()

