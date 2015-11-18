#! /usr/bin/env python 

import pyfits
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

data_path = "/home/ycli/data/gbt/raw/3c286/"
data_list = [
        '41_3C286_onoff_124-125.fits',
        #'42_3C286_onoff_328-329.fits',
        #'85_3C286_onoff_22-23.fits',
        #'84_3C286_onoff_8-9.fits',
        ]

plt.figure(figsize=(6,6))

spec = None

for data_fits in data_list:

    hdulist = pyfits.open(data_path + data_fits)

    tbdata = hdulist[1].data

    hdulist.close()

    print hdulist
    print tbdata

    for key in hdulist[0].header.keys():
        print key, hdulist[0].header[key]

    for key in hdulist[1].header.keys():
        print key, hdulist[1].header[key]

    fieldlabel = []

    for i in range(hdulist[1].header['TFIELDS']):
        fieldlabel.append(hdulist[1].header['TTYPE%d'%(i+1)])

    for i in range(hdulist[1].header['TFIELDS']):
        print hdulist[1].header['TTYPE%d'%(i+1)]
        print sp.unique(tbdata.field(fieldlabel[i])).shape
        print tbdata.field(fieldlabel[i]).shape
        print tbdata.field(fieldlabel[i])[:16]
        print
    exit()

    data = tbdata.field('DATA')
    data = data.reshape(2, -1, 4, 2, 4096)

    if spec == None:
        spec_crval = tbdata.field('CRVAL1')[0]
        spec_crpix = tbdata.field('CRPIX1')[0]
        spec_cdelt = tbdata.field('CDELT1')[0]
        spec_shape = data.shape[-1]

        spec = ( np.arange(spec_shape) - spec_crpix + 1 ) * spec_cdelt + spec_crval
        spec /= 1.e6

    #quasar_spec_on = data[0,:,0,0,:] - data[1,:,0,0,:]
    #quasar_spec_off = data[0,:,0,1,:] - data[1,:,0,1,:]
    #plt.plot(spec, np.mean(quasar_spec_on, axis=0), '.-')
    #plt.plot(spec, np.mean(quasar_spec_off,axis=0), '.-')

    #plt.plot(tbdata.field('RA')[0], tbdata.field('DEC')[0], 'o-')
    #plt.plot(tbdata.field('RA')[432], tbdata.field('DEC')[432], 'o-')
    #plt.show()


    label = '%s [ RA: %12.6f DEC %12.6f ]'%(tbdata.field('DATE-OBS')[0],
                                            tbdata.field('RA')[0], 
                                            tbdata.field('DEC')[0], )

    quasar_spec = data[0,:,0,:,:] - data[1,:,0,:,:]
    quasar_spec = np.mean(quasar_spec, axis=1)

    #print np.mean(quasar_spec, axis=0)
    #print spec
    plt.plot(spec, np.mean(quasar_spec, axis=0), '.-', label=label)

    spec_selected = spec[np.logical_and(spec>839,spec<840)]
    quasar_selected = np.mean(quasar_spec, axis=0)[np.logical_and(spec>839,spec<840)]
    output = np.append(spec_selected[:,None], quasar_selected[:,None], axis=1)
    np.savetxt('./%s_RA%9.6f_DEC%8.6f'%(tbdata.field('DATE-OBS')[0],
        tbdata.field('RA')[0],tbdata.field('DEC')[0]), output, fmt='%10.6f')


#plt.ylim(ymin=18, ymax=20)
plt.ylim(ymin=15, ymax=38)
plt.xlim(xmin=839, xmax=840)
plt.legend(frameon=False)
plt.grid()
plt.show()
