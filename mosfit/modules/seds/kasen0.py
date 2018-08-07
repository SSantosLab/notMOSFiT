"""Definitions for the `Kasen0` class."""
import numpy as np
from mosfit.modules.seds.sed import SED
import pickle 
from astropy import constants as c
from astropy import units as u
import os



# Important: Only define one ``Module`` class per file.


class Kasen0(SED):
    '''
	Defining the Kasen-simulation based SED

    FOR TYPE 0 == SHOCK HEATED EJECTA
    #Frankencode
    Kamile Lukosiute August 2018
    What is my life
    '''

    # Kasen-calculated parameters
    MASS = np.array([.001, .0025, .005, .01, .02, .025, .03, .04, .05, .075, .1])
    VKIN = np.array([.03, .05, .1, .2, .3])
    XLAN = np.array([1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-9])

    MASS_S = np.array(['0.001', '0.0025', '0.005', '0.010' ,'0.020', '0.025', '0.030', '0.040', '0.050', '0.075', '0.1' ])
    VKIN_S = np.array(['0.03', '0.05', '0.10', '0.20', '0.30'])
    XLAN_S = np.array(['1e-1', '1e-2', '1e-3', '1e-4', '1e-5', '1e-9'])

    C_CONST = c.c.cgs.value

    def __init__(self, **kwargs):
        super(Kasen0, self).__init__(**kwargs)

        # Read in times and frequencies arrays (same for all SEDs)
        self._dir_path = os.path.dirname(os.path.realpath(__file__))
        self._kasen_frequencies = pickle.load( open(os.path.join(self._dir_path, 'kasen_seds/frequency_angstroms.p'), "rb"))
        self._kasen_times = pickle.load( open(os.path.join(self._dir_path, 'kasen_seds/times_days.p'), "rb"))



    def process(self, **kwargs):
        lum_key = self.key('luminosities')
        kwargs = self.prepare_input(lum_key, **kwargs)
        self._luminosities = kwargs[lum_key]

        self._times = kwargs['times_to_proc']
#        times_key = self.key('times_to_proc')
#        kwargs = self.prepare_input(times_key, **kwargs)
#        self._times = kwargs[times_key]
#        self._dense_times = kwargs[self.key('dense_times')]
        print(len(self._times))
        print(len(self._luminosities))
       # print(len(self._dense_times))
        '''
            Okay what have we learned so far? There's some sort of weird data
            format within MOSFiT, where  dense_times is evenly spaced in log 
            space. BUT! When processing luminosities, it's not exactly 1 to 1
            for the times. Instead, each t in dense_time has an array that it
            corresponds to called dense_indices, which tells you which index
            in the dense_lums array that time corresponds to... 
            18:02 EST 2018.08.07
        '''
        '''     print('times')
        print(self._times)
        print(len(self._times))
        print('dense_times')
        print(self._dense_times)
        print(self._dense_times.shape)


        print('lums')
        print(self._luminosities.shape)
        print('indices')
        print(self._indices)
        print(self._indices.shape)
        print('rest times')
        print(self._rest_times)
        print(self._rest_times.shape)



        [
            np.inf
            if self._rest_texplosion > x else (x - self._rest_texplosion)
            for x in self._dense_times
        ]
        print('_times_to_proc')
        print(self._times_to_proc)
        print(len(self._times_to_proc))

        '''

        self._band_indices = kwargs['all_band_indices']
        self._frequencies = kwargs['all_frequencies']


        # Physical parameters from Kasen simulations, provided by
        # neutrinosphere module (thank you, Jessica Metzger)
        self._vk = kwargs[self.key('vk')]
        self._xlan = kwargs[self.key('xlan')]
        self._mass = kwargs[self.key('mejecta')]
        
        # mass fractional weight provided by neutrinoshere module
        self._mass_weight = kwargs[self.key('mass_weight')]
        
        # viewing angle and opening angle 
        self._phi = kwargs[self.key('phi')]
        self._theta = kwargs[self.key('theta')]

        # Total weight function
        # TYPE 0 == SHOCK HEATED SO GEOMETRIC FACTOR IS JUST FOR CONE
        weight_goem = 2*np.cos(self._theta)*(1. - np.cos(self._phi))
        weight = self._mass_weight * weight_goem

        # Some temp vars for speed.
        cc = self.C_CONST

        zp1 = 1.0 + kwargs[self.key('redshift')]
        Azp1 = u.Angstrom.cgs.scale / zp1
        czp1 = cc / zp1


        seds = []
        rest_wavs_dict = {}

        # Find nearest neighbors to the Kasen-calculated simulation
        m_closest = self.MASS_S[(np.abs(self.MASS-self._mass)).argmin()]
        v_closest = self.VKIN_S[(np.abs(self.VKIN-self._vk)).argmin()]
        x_closest = self.XLAN_S[(np.abs(self.XLAN-self._xlan)).argmin()]

        # Open nearest neighbor file
        kasen_seds = pickle.load( open(os.path.join(self._dir_path, 'kasen_seds/knova_d1_n10_m' + m_closest + '_vk' + v_closest + '_fd1.0_Xlan' + x_closest + '.0.p', ), "rb" ))

        # For each time (luminosities as proxy)
        for li, lum in enumerate(self._luminosities):
            bi = self._band_indices[li]
            if bi >= 0:
                rest_wavs = rest_wavs_dict.setdefault(
                    bi, self._sample_wavelengths[bi] * Azp1)
                # bi, self._sample_wavelengths[bi] * 10.)
            else:
                rest_wavs = np.array(  # noqa: F841
                    [czp1 / self._frequencies[li]])

            # Find corresponding closest time

            t_closest_i = (np.abs(self._kasen_times-self._times[li])).argmin()

            # Evaluate the SED at the rest frame frequencies
            sed = np.array([])
            for w in rest_wavs:
                # find index of closest wav
                w_closest_i = np.abs(self._kasen_frequencies-w).argmin()

                sed = np.append(sed, weight * kasen_seds['SEDs'][t_closest_i][w_closest_i] )
        #print(t_closest_i, w_closest_i)
            seds.append(sed)
            seds[-1][np.isnan(seds[-1])] = 0.0
        

        seds = self.add_to_existing_seds(seds, **kwargs)

        return {'sample_wavelengths': self._sample_wavelengths, 'seds': seds}
