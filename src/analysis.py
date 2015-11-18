import h5py
import numpy as np
import matplotlib.pyplot as plt

def flag_freq(spec, freq, range=[[0,1], ]):

    spec = np.ma.array(spec)

    for r in range:

        bad = np.logical_and(freq > r[0], freq < r[1]) 
        spec[:, bad] = np.ma.masked

    return spec


def flag_RFI(spec):

    spec = np.ma.array(spec)
    spec[np.logical_not(np.isfinite(spec))] = np.ma.masked

    sigma_raw = np.ma.std(spec)

    while 1:

        spec_timemean = np.ma.mean(spec, axis=0)
        spec_mean = np.ma.mean(spec_timemean)
        spec_std = np.ma.std(spec_timemean)
        #print spec_mean, spec_std

        spec_bad = np.logical_or(spec_timemean > spec_mean + 6*spec_std,
                spec_timemean < spec_mean - 6*spec_std)
        spec[:, spec_bad] = np.ma.masked

        sigma = np.ma.std(spec)
        if sigma_raw - sigma > 0.01 * sigma_raw:
            sigma_raw = sigma
        else:
            #bad = spec.mask
            return spec


def pol_rotation(spec):

    spec_I = spec[:, 0, :, :] + spec[:, 1, :, :]

    return spec_I

def check_noise_cal(spec, freq):

    spec = np.ma.array(spec)
    spec[np.logical_not(np.isfinite(spec))] = np.ma.masked

    cal = spec[:,1,:] - spec[:,0,:]

    #cal[np.logical_not(np.isfinite(cal))] = np.ma.masked

    plt.plot(freq, np.ma.mean(cal, axis=0))
    plt.plot(freq, np.ma.mean(spec[:,0,:], axis=0))
    plt.show()

def cal_unit(spec):

    # subtract noise cal, and average over noise on and off

    cal = spec[:,1,:] - spec[:,0,:]

    spec[:,1,:] -= cal

    spec /= cal[:, None, :]

    return np.mean(spec, axis=1)

def check_spec_1d(spec_list, freq, label_list=['',], output='', 
        ymax=0.80, ymin=0.55, flag_freq=None, sigma=None):

    color_list = ['g', 'r', 'k', 'b', 'c']

    fig = plt.figure(figsize=(6,3))
    ax = fig.add_axes([0.14, 0.16, 0.83, 0.80])

    for i in range(len(spec_list)):
        spec = spec_list[i]
        label = label_list[i]

        spec = np.ma.array(spec)
        spec[np.logical_not(np.isfinite(spec))] = np.ma.masked

        #ax.step(freq/1.e6, np.ma.mean(spec[:,0,:], axis=0))
        #ax.step(freq/1.e6, np.ma.mean(spec[:,1,:], axis=0))
        ax.step(freq/1.e6, np.ma.mean(spec, axis=0), color_list[i], 
                where='mid', label=label, linewidth=2)

    if flag_freq != None:
        for flag in flag_freq:
            ax.vlines(flag/1.e6, ymax=ymax, ymin=ymin, colors='r')
    if sigma != None:
        ax.hlines(sigma[0], xmin=freq.min()/1.e6, xmax=freq.max()/1.e6, 
                colors='r', linestyles='dotted', linewidths=2)
        ax.hlines(sigma[1], xmin=freq.min()/1.e6, xmax=freq.max()/1.e6, 
                colors='r', linestyles='dotted', linewidths=2)

    ax.set_xlim(xmin=freq.min()/1.e6, xmax=freq.max()/1.e6)
    ax.set_ylim(ymin=ymin, ymax=ymax)

    ax.set_xlabel('Frequency [MHz]')
    ax.set_ylabel('Cal Unit')

    ax.legend(frameon=False, loc=0)
    ax.minorticks_on()
    ax.tick_params(length=6, width=1., direction='out')
    ax.tick_params(which='minor', length=3, width=1., direction='out')

    plt.savefig(output + '.png', format='png')

    plt.show()


def check_data(data_path, data_name, flag_range=[[0, 1],], output_root=''):

    data = h5py.File(data_path + data_name + '.hdf5', 'r')

    band_list = data.keys()
    for band in band_list[:1]:
        freq = data[band + '/freq'].value
        spec_calon  = None
        spec_caloff = None
        for object in data[band].keys()[::-1]:
            print object
            if object[3:] == '3C286_ON':
                if spec_calon == None:
                    spec_calon = pol_rotation(data[band + '/' + object].value)
                else:
                    spec_calon = np.concatenate([
                        pol_rotation(data[band + '/' + object].value), spec_calon],
                        axis=0)
            if object[3:] == '3C286_OFF':
                if spec_caloff == None:
                    spec_caloff = pol_rotation(data[band + '/' + object].value)
                else:
                    spec_caloff = np.concatenate([
                        pol_rotation(data[band + '/' + object].value), spec_calon],
                        axis=0)
            if object[3:] == 'absorber_source_ON':
                spec = pol_rotation(data[band + '/' + object].value)
                #check_noise_cal(spec, freq)
                #check_spec_1d([spec,], freq, 
                #        label_list=['time averaged spec',],
                #        output= output_root + 'abs_timeavg')
                #spec = flag_freq(spec, freq/1.e6, range=flag_range)
                #spec = flag_RFI(spec)
                #check_spec_1d([spec,], freq, 
                #        label_list=['time averaged spec, RFI flagged',],
                #        output= output_root + 'abs_timeavg_RFIflagged')
                #print spec.shape
            if object[3:] == 'absorber_ref_ON':
                spec_ref = pol_rotation(data[band + '/' + object].value)
            
                #check_noise_cal(spec_ref, freq)
                #check_spec_1d([spec_ref,], freq, 
                #        label_list=['time averaged spec',],
                #        output= output_root + 'ref_timeavg')
                #spec_ref = flag_freq(spec_ref, freq/1.e6, range=flag_range)
                ##spec_ref = flag_RFI(spec_ref)
                #spec_ref.mask = spec.mask
                #check_spec_1d([spec_ref,], freq, 
                #        label_list=['time averaged spec, RFI flagged',],
                #        output= output_root + 'ref_timeavg_RFIflagged')
                #print spec.shape

        spec     = cal_unit(spec)
        spec_ref = cal_unit(spec_ref)
        spec_sub = spec - spec_ref
        spec_sub = flag_freq(spec_sub, freq/1.e6, range=flag_range)
        #spec_sub = flag_RFI(spec_sub)
        check_spec_1d([spec_sub,], freq, 
                label_list=['time averaged spec, RFI flagged, sub ref',],
                ymax=0.2, ymin=0.0,
                output = output_root + 'abs_timeavg_RFIflagged_subref')

        spec_sub.dump(data_path + data_name + '_%s_absorber_spec.npy'%band)
        np.save(data_path + data_name + '_%s_absorber_freq.npy'%band, freq)

        #spec_caloff = cal_unit(spec_caloff)
        #spec_calon = cal_unit(spec_calon)
        #spec_cal = spec_calon - spec_caloff
        #spec_cal = flag_freq(spec_cal, freq/1.e6, range=flag_range)
        #check_spec_1d([spec_cal,], freq, 
        #        label_list=['time averaged spec, 3C286',],
        #        ymax=0.8, ymin=0.3,
        #        output = output_root + 'abs_timeavg_3C286')
        #spec_cal.dump(data_path + data_name + '_%s_3C286_spec.npy'%band)
        #np.save(data_path + data_name + '_%s_3C286_freq.npy'%band, freq)
    data.close()

def spec_smoothing(spec):

    import scipy.signal as signal

    spec_mean = np.ma.mean(spec, axis=0)
    spec_smooth = signal.medfilt(spec_mean, 101)
    spec_smooth = np.ma.array(spec_smooth)
    spec_smooth.mask = spec_mean.mask

    return spec_smooth

def find_absorber(data_path, data_name):

    spec = np.load(data_path + data_name + '_spec.npy')
    freq = np.load(data_path + data_name + '_freq.npy')

    spec_mean = np.ma.mean(spec, axis=0)
    spec_smooth = spec_smoothing(spec)

    spec_mean -= spec_smooth

    std  = np.ma.std(spec_mean)
    mean = np.ma.mean(spec_mean)

    abs_freq = freq[np.logical_and(
        spec_mean < mean - 6 * std, np.logical_not(spec_mean.mask))]

    print abs_freq.shape

    check_spec_1d( [spec, spec_smooth[None, ...]], 
        freq, 
        label_list=['target','smoothed'], 
        ymax=0.3, ymin=0.0, 
        output = './png/%s_fullspec'%(data_name,),)

    for abs in abs_freq:

        delta_freq = 0.15 * 1.e6

        freq_range = np.logical_and(
            freq > abs - delta_freq, freq < abs + delta_freq)

        check_spec_1d( [ (spec - spec_smooth[None, ...])[..., freq_range], 
            ((spec_smooth-spec_smooth)[None, ...])[..., freq_range] ], 
            freq[..., freq_range] - abs, 
            label_list=['3c286 - smoothed [%6.2fMHz]'%(abs/1.e6), ''], 
            ymax=0.06, ymin=-0.06, 
            output = './png/%s_abs%07dkHz'%(data_name, int(abs/1.e3)), 
            flag_freq=[0,], 
            sigma=[mean-std, mean+std])


if __name__=="__main__":

    data_path = './data/'
    data_name = 'AGBT15A_196_02'
    flag_range = [[600, 720], [880, 900], [790, 810], [813, 823]] #[790, 825]
    #check_data(data_path, data_name, flag_range=flag_range, output_root='./png/')
    find_absorber(data_path, data_name + '_band_A_absorber')
    #find_absorber(data_path, data_name + '_band_A_3C286')



