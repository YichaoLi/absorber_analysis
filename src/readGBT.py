#! /bin/env python 

import numpy as np
import scipy as sp
from scipy.signal import medfilt
import pyfits
import h5py

def check_fits(data_path):

    hdulist = pyfits.open(data_path, memmap=False)
    print hdulist

    mheader = hdulist[1].header

    for key in mheader.keys():
        print key, '\t', mheader[key]
    print mheader['STT_IMJD']
    print mheader['STT_SMJD']
    print mheader['STT_OFFS']

    exit()

    for k in range(1, 2):
        tbdata = hdulist[k].data

        fieldlabel = []
        for i in range(hdulist[k].header['TFIELDS']): 
            fieldlabel.append(hdulist[k].header['TTYPE%d'%(i+1)])
        print fieldlabel
        #for i in range(len(tbdata)):
        #    print tbdata[i][fieldlabel[k]]
        for i in range(hdulist[k].header['TFIELDS']):
            print hdulist[k].header['TTYPE%d'%(i+1)]
            print sp.unique(tbdata.field(fieldlabel[i])).shape
            print tbdata.field(fieldlabel[i]).shape
            print tbdata.field(fieldlabel[i])
            print

def read_fits(hdulist):

    mheader = hdulist[0].header
    dheader = hdulist[1].header
    freq = hdulist[1].data[0]['DAT_FREQ'].astype(float)
    delta_t = dheader['TBIN']
    nfreq = dheader['NCHAN']
    delta_f = dheader['CHAN_BW']
    freq0 = mheader['OBSFREQ'] - mheader['OBSBW'] / 2. + delta_f / 2
    time0 = mheader["STT_SMJD"] + mheader["STT_OFFS"]

    print delta_t

    record_ra  = hdulist[1].data.field('RA_SUB')
    record_dec = hdulist[1].data.field('DEC_SUB')

    nrecords = len(hdulist[1].data)
    print nrecords
    record_ind = np.arange(nrecords)
    nrecords_read = len(record_ind)

    ntime_record, npol, nfreq, one = hdulist[1].data[0]["DATA"].shape
    print ntime_record, npol, nfreq, one
    out_data = np.empty((nfreq, npol, nrecords_read, ntime_record), dtype=np.float32)

    for ii in range(nrecords_read):
        index = record_ind[ii] 
        #print record_ra[index], record_dec[index]
        record = hdulist[1].data[index]["DATA"]
        scl = hdulist[1].data[index]["DAT_SCL"]
        scl.shape = (1, npol, nfreq, 1)
        offs = hdulist[1].data[index]["DAT_OFFS"]
        offs.shape = (1, npol, nfreq, 1)
        record *= scl
        record += offs
        # Interpret as unsigned int (for Stokes I only).
        #record = record.view(dtype=np.uint8)
        # Select stokes I and copy.
        out_data[:,0,ii,:] = np.transpose(record[:,0,:,0])
        # Interpret as signed int (except Stokes I).
        record = record.view(dtype=np.int8)
        out_data[:,1,ii,:] = np.transpose(record[:,1,:,0])
        out_data[:,2,ii,:] = np.transpose(record[:,2,:,0])
        out_data[:,3,ii,:] = np.transpose(record[:,3,:,0])
    out_data.shape = (nfreq, npol, nrecords_read * ntime_record)

    #time = np.arange(data.shape[2]) * delta_t + time0
    time = np.arange(record_ind[0]*ntime_record, (record_ind[-1] + 1)*ntime_record)\
            * delta_t + time0
    subint_time = time0 + hdulist[1].data['OFFS_SUB']
    ra = sample_subint(subint_time, hdulist[1].data['RA_SUB'], time)
    dec = sample_subint(subint_time, hdulist[1].data['DEC_SUB'], time)
    az = sample_subint(subint_time, hdulist[1].data['TEL_AZ'], time)
    el = 90. - sample_subint(subint_time, hdulist[1].data['TEL_ZEN'], time)

    return out_data, time, ra, dec, az, el, freq

def sample_subint(sub_time, sub_var, time):

    diff = np.diff(sub_var) / np.diff(sub_time)
    rate = np.mean(diff)
    start_ind = np.argmin(np.abs(time - sub_time[0]))
    return (time - time[start_ind]) * rate + sub_var[0]

def read_raw(data_path, data_name):

    hdulist = pyfits.open(data_path + data_name + '.fits', 'readonly', memmap=False)
    data, time, ra, dec, az, el, freq = read_fits(hdulist)

    timestream_data = Data(data, time, freq, ra, dec)

    return timestream_data

def write(timestream_data, output_path):

    data_file = h5py.File(output_path + '.hdf5', 'w')
    data_file['data'] = timestream_data.data
    data_file['time'] = timestream_data.time
    data_file['freq'] = timestream_data.freq
    data_file['ra']   = timestream_data.ra
    data_file['dec']  = timestream_data.dec
    #data_file['az']   = az 
    #data_file['el']   = el

    data_file.close()

def load(data_path):

    data_file = h5py.File(data_path + '.hdf5', 'r')

    timestream_data = Data(
            data_file['data'], 
            data_file['time'], 
            data_file['ferq'], 
            data_file['ra'], 
            data_file['dec'])

    data_file.close()

    return timestream_data

class TimestreamData(object):


    def __init__(self, data, time, freq, ra, dec):

        self.data = data
        self.time = time
        self.freq = freq
        self.ra   = ra
        self.dec  = dec


if __name__=="__main__":

    data_path  = '/project/ycli/data/gbt/AGBT14B_339/01_20140718/'
    data_name  = 'guppi_56856_wigglez1hr_centre_0011_0001'

    output_path = '/project/ycli/data/gbt/converted/'

    #check_fits(data_path + data_name + '.fits')

    data = read_raw(data_path, data_name)

    write(data, output_path + data_name)
    data = load(output_path + data_name)


