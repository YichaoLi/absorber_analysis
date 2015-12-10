#! /usr/bin/env python 

import pyfits
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from core import fitsGBT

data_path = "/project/ycli/data/gbt/3C286/"
data_list = [
        '/GBT14B_339/20140817_07_3C286_onoff_26-27.fits',
        '/GBT14B_339/20140913_10_3C286_onoff_12-13.fits',
        '/GBT14B_339/20140920_11_3C286_onoff_10-11.fits',
        '/GBT14B_339/20141219_28_3C286_onoff_39-40.fits'
        ]

fig = plt.figure(figsize=(6, 4))
ax = fig.add_axes([0.12, 0.12, 0.8, 0.8])

colorlist = ['r', 'b', 'g', 'k', 'm',]
ii = 0

for data_fits in data_list:

    data_block = fitsGBT.Reader(data_path + data_fits).read()
    
    data_all = []
    i = 0
    for block in data_block:
        print block.data.shape

        data = block.data
        data_all.append(data)
        print block.history
        print data.shape

        spec_crval = block.field['CRVAL1']
        spec_crpix = block.field['CRPIX1']
        spec_cdelt = block.field['CDELT1']
        spec_shape = data.shape[-1]

        print spec_crval 
        print spec_crpix
        print spec_cdelt
        print spec_shape

        spec = ( np.arange(spec_shape) - spec_crpix + 1) * spec_cdelt + spec_crval
        spec /= 1.e6


        #fig, ax = plt.subplots(nrows=4, ncols=2, sharex=True, sharey=True, 
        #        figsize=(10,9))
        #plt.subplots_adjust(left=0.1, right=0.95, bottom=0.05, top=0.90,
        #                    wspace=0.01, hspace=0.01)

        #cax = fig.add_axes([0.2, 0.97, 0.6, 0.01])

        ##cmin = np.min(data)
        ##cmax = np.max(data)
        #cmin = -40
        #cmax = 40

        #for i in range(2):
        #    im = ax[0, i].imshow(data[:,0,i,:].T, 
        #                         interpolation='nearest', 
        #                         origin='lower',
        #                         extent=[0,data.shape[0],spec[0],spec[-1]],
        #                         aspect='auto')
        #    im.set_clim(cmin, cmax)

        #    im = ax[1, i].imshow(data[:,1,i,:].T, 
        #                         interpolation='nearest', 
        #                         origin='lower',
        #                         extent=[0,data.shape[0],spec[0],spec[-1]],
        #                         aspect='auto')
        #    im.set_clim(cmin, cmax)

        #    im = ax[2, i].imshow(data[:,2,i,:].T, 
        #                         interpolation='nearest', 
        #                         origin='lower',
        #                         extent=[0,data.shape[0],spec[0],spec[-1]],
        #                         aspect='auto')
        #    im.set_clim(cmin, cmax)

        #    im = ax[3, i].imshow(data[:,3,i,:].T, 
        #                         interpolation='nearest', 
        #                         origin='lower',
        #                         extent=[0,data.shape[0],spec[0],spec[-1]],
        #                         aspect='auto')
        #    im.set_clim(cmin, cmax)

        #ax[0, 0].set_title('noise on')
        #ax[0, 1].set_title('noise off')
        #ax[0, 0].set_ylabel('freq [MHz] I')
        #ax[1, 0].set_ylabel('freq [MHz] Q')
        #ax[2, 0].set_ylabel('freq [MHz] U')
        #ax[3, 0].set_ylabel('freq [MHz] V')

        #fig.colorbar(im, cax=cax, orientation = 'horizontal')

        #plt.tick_params(length=6, width=1.)
        #plt.tick_params(which='minor', length=3, width=1.)
        #plt.savefig('./png/3c286_block%03d.png'%i, format='png')
        #i += 1
        ##plt.show()

        #plt.close()

    data_quasar = data_all[0] - data_all[1]

    cal = np.mean(data_quasar[:, 0, 0, :] - data_quasar[:, 0, 1, :], axis=0)
    data_quasar[:, 0, 0, :] -= cal[None, :]
    #data_quasar /= cal[None, :]

    #ax.step(spec, np.mean(data_quasar[:, 0, 0, :], axis=0), 
    #        label=data_fits + '\n cal_on')
    ax.step(spec, np.mean(data_quasar[:, 0, 1, :], axis=0), 
            color=colorlist[ii], linewidth=2, where='mid',
            label=data_fits.split('_')[0] + data_fits.split('_')[1] + ' cal_off')
    ax.plot(spec, np.mean(data_quasar[:, 0, 1, :], axis=0), '.',
            color=colorlist[ii])

    ii += 1

    #data_quasar = data_all[0]# - data_all[1]
    ##print spec[1] - spec[0]
    #plt.plot(spec, np.mean(data_quasar[:,0,0,:], axis=0), '.-')
    #plt.plot(spec, np.mean(data_quasar[:,0,1,:], axis=0), '.-')
    ##plt.plot(spec, np.mean(data_quasar[:,1,0,:], axis=0))
    #data_quasar = data_all[1]# - data_all[1]
    #plt.plot(spec, np.mean(data_quasar[:,0,0,:], axis=0), '.-')
    #plt.plot(spec, np.mean(data_quasar[:,0,1,:], axis=0), '.-')
    ##plt.plot(spec, np.mean(data_quasar[:,0,0,:], axis=0))
    ##plt.plot(spec, np.mean(data_quasar[:,1,0,:], axis=0))
    #plt.xlim(xmin=835.0, xmax=845.0)
    ##plt.ylim(ymin=28, ymax=32)
    ##plt.ylim(ymin=1, ymax=4)
    #plt.show()

ax.set_xlim(xmin=838.5, xmax=840.0)
ax.set_ylim(ymin=16, ymax=30)

ax.legend(frameon=False)

ax.minorticks_on()
ax.tick_params(length=4, width=1, direction='out')
ax.tick_params(which='minor', length=2, width=1, direction='out')

plt.show()
