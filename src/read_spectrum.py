#! /usr/bin/env python 

import pyfits
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

data_path = "/mnt/raid-project/gmrt/kiyo/data/guppi_data/"
data_list = [

        '/GBT14B_339/07_3C286_onoff_24-25.fits',
        '/GBT14B_339/07_3C286_onoff_26-27.fits',
        '/GBT14B_339/10_3C286_onoff_10-11.fits',
        '/GBT14B_339/10_3C286_onoff_12-13.fits',
        '/GBT14B_339/11_3C286_onoff_10-11.fits',
        '/GBT14B_339/11_3C286_onoff_12-13.fits',
        '/GBT14B_339/28_3C286_onoff_39-40.fits',
        '/GBT14B_339/28_3C286_onoff_41-42.fits',
        '/GBT14B_339/31_3C286_onoff_6-7.fits',
        '/GBT14B_339/31_3C286_onoff_8-9.fits',
        '/GBT14B_339/38_3C286_onoff_20-21.fits',
        '/GBT14B_339/38_3C286_onoff_22-23.fits',
        '/GBT14B_339/39_3C286_onoff_20-21.fits',
        '/GBT14B_339/39_3C286_onoff_22-23.fits',
        '/GBT14B_339/39_3C286_onoff_38-39.fits',
        '/GBT14B_339/39_3C286_onoff_40-41.fits',

        #'/GBT13B_352/21_3C286_onoff_23-24.fits',
        #'/GBT13B_352/21_3C286_onoff_25-26.fits',
        #'/GBT13B_352/23_3C286_onoff_22-23.fits',
        #'/GBT13B_352/23_3C286_onoff_24-25.fits',
        #'/GBT13B_352/29_3C286_onoff_12-13.fits',
        #'/GBT13B_352/29_3C286_onoff_14-15.fits',
        #'/GBT13B_352/31_3C286_onoff_11-12.fits',
        #'/GBT13B_352/31_3C286_onoff_13-14.fits',
        #'/GBT13B_352/32_3C286_onoff_37-38.fits',
        #'/GBT13B_352/32_3C286_onoff_39-40.fits',

        #'/GBT12A_418/19_3C286_track_36.fits',
        #'/GBT12A_418/19_3C286_track_37.fits',
        #'/GBT12A_418/19_3C286_track_38.fits',
        #'/GBT12A_418/19_3C286_track_39.fits',
        #'/GBT12A_418/19_3C286_track_40.fits',
        #'/GBT12A_418/19_3C286_track_41.fits',
        #'/GBT12A_418/19_3C286_track_42.fits',
        #'/GBT12A_418/19_3C286_track_43.fits',
        #'/GBT12A_418/19_3C286_track_62.fits',
        #'/GBT12A_418/19_3C286_track_63.fits',
        #'/GBT12A_418/19_3C286_track_64.fits',
        #'/GBT12A_418/19_3C286_track_65.fits',
        #'/GBT12A_418/19_3C286_track_66.fits',
        #'/GBT12A_418/19_3C286_track_67.fits',
        #'/GBT12A_418/19_3C286_track_68.fits',
        #'/GBT12A_418/19_3C286_track_69.fits',
        #'/GBT12A_418/19_3C286_track_71.fits',
        #'/GBT12A_418/19_3C286_track_72.fits',
        #'/GBT12A_418/19_3C286_track_73.fits',
        #'/GBT12A_418/19_3C286_track_74.fits',
        #'/GBT12A_418/19_3C286_track_75.fits',
        #'/GBT12A_418/19_3C286_track_76.fits',
        #'/GBT12A_418/19_3C286_track_77.fits',
        #'/GBT12A_418/19_3C286_track_78.fits',
        #'/GBT12A_418/21_3C286_track_18.fits',
        #'/GBT12A_418/21_3C286_track_19.fits',
        #'/GBT12A_418/21_3C286_track_20.fits',
        #'/GBT12A_418/21_3C286_track_21.fits',
        #'/GBT12A_418/21_3C286_track_22.fits',
        #'/GBT12A_418/21_3C286_track_23.fits',
        #'/GBT12A_418/21_3C286_track_24.fits',
        #'/GBT12A_418/21_3C286_track_25.fits',
        #'/GBT12A_418/21_3C286_track_33.fits',
        #'/GBT12A_418/21_3C286_track_34.fits',
        #'/GBT12A_418/21_3C286_track_35.fits',
        #'/GBT12A_418/21_3C286_track_36.fits',
        #'/GBT12A_418/21_3C286_track_37.fits',
        #'/GBT12A_418/21_3C286_track_38.fits',
        #'/GBT12A_418/21_3C286_track_39.fits',
        #'/GBT12A_418/21_3C286_track_40.fits',
        #'/GBT12A_418/21_3C286_track_52.fits',
        #'/GBT12A_418/21_3C286_track_53.fits',
        #'/GBT12A_418/21_3C286_track_54.fits',
        #'/GBT12A_418/21_3C286_track_55.fits',
        #'/GBT12A_418/21_3C286_track_56.fits',
        #'/GBT12A_418/21_3C286_track_57.fits',
        #'/GBT12A_418/21_3C286_track_58.fits',
        #'/GBT12A_418/21_3C286_track_59.fits',
        ]

plt.figure(figsize=(6,6))

spec = None

for data_fits in data_list:

    hdulist = pyfits.open(data_path + data_fits)

    tbdata = hdulist[1].data

    hdulist.close()

    #print hdulist
    #print tbdata

    print data_fits
    print tbdata.field('DATE-OBS')[0]
    print tbdata.field('DATE-OBS')[-1]
    print 

    continue

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

    continue

    data = tbdata.field('DATA')
    data = data.reshape(2, -1, 4, 2, 4096)

    if spec == None:
        spec_crval = tbdata.field('CRVAL1')[0]
        spec_crpix = tbdata.field('CRPIX1')[0]
        spec_cdelt = tbdata.field('CDELT1')[0]
        spec_shape = data.shape[-1]

        spec = ( np.arange(spec_shape) - spec_crpix ) * spec_cdelt + spec_crval
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

    #spec_selected = spec[np.logical_and(spec>839,spec<840)]
    #quasar_selected = np.mean(quasar_spec, axis=0)[np.logical_and(spec>839,spec<840)]
    #output = np.append(spec_selected.T, quasar_selected)
    #np.savetxt(output, './%s_RA%.6f_DEC%.6f'%(tbdata.field('DATE-OBS')[0],
    #    tbdata.field('RA')[0],tbdata.field('DEC')[0],))


#plt.ylim(ymin=18, ymax=20)
#plt.ylim(ymin=15, ymax=38)
#plt.xlim(xmin=839, xmax=840)
#plt.legend(frameon=False)
#plt.grid()
#plt.show()
