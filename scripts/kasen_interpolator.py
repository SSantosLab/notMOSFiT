import numpy as np
import matplotlib.pyplot as plt
import bisect
import h5py

# Kasen Parameters
mass = np.array([.001]) #, .0025, .005, .01, .02, .025, .03, .04, .05, .075, .1])
mass_s = np.array(['0.001']) #, '0.0025', '0.005', '0.010' ,'0.020', '0.025', '0.030', '0.040', '0.050', '0.075', '0.1' ])
vk = np.array([.03]) #, .05, .1, .2, .3])
vk_s = np.array(['0.03']) #, '0.05', '0.10', '0.20', '0.30'])
Xlan = np.array([1e-1]) #, 1e-2, 1e-3, 1e-4, 1e-5, 1e-9])
Xlan_s = np.array(['1e-1']) #, '1e-2', '1e-3', '1e-4', '1e-5', '1e-9'])


XLAN = np.array([])
VKIN = np.array([])
MASS = np.array([])
TIME  = np.array([])
FREQ = np.array([])
VALS = np.array([])

path = '../Kasen_Kilonova_Models_2017/systematic_kilonova_model_grid/'
for mi, m in enumerate(mass_s):
    for vi, v in enumerate(vk_s):
        for xi, x in enumerate(Xlan_s):
            # Read in files for each parameter
            # -----------------------------------------------------------------
            fname = path + 'knova_d1_n10_m' + m + '_vk' + v + '_fd1.0_Xlan' + x + '.0.h5'
            fin    = h5py.File(fname,'r')

            # frequency in Hz
            nu    = np.array(fin['nu'],dtype='d')
            # array of time in seconds
            # throw out negative values and things > 14.0
            times = np.array(fin['time'])
            times = times/3600.0/24.0
            times = times[(times >= 0.0) & (times <= 14.0)]

            # specific luminosity (ergs/s/Hz) 
            # this is a 2D array, Lnu[times][nu]
            Lnu_all   = np.array(fin['Lnu'],dtype='d')

            for t in times:
                # index corresponding to t
                it = bisect.bisect(times,t)
                # spectrum at this epoch
                Lnu = Lnu_all[it,:]

                for ni,n in enumerate(nu):
                    XLAN = np.append(XLAN, Xlan[xi])
                    VKIN = np.append(VKIN, vk[vi])
                    MASS = np.append(MASS, mass[mi])
                    TIME = np.append(TIME, t)
                    FREQ = np.append(FREQ, n)
                    VALS = np.append(VALS, Lnu[ni])