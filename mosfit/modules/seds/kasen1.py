"""Definitions for the `Blackbody` class."""
import numpy as np
from mosfit.modules.seds.sed import SED
import pickle 
from astropy import constants as c
from astropy import units as u



# Important: Only define one ``Module`` class per file.


class Kasen1(SED):
    '''
	Defining the Kasen-simulation based SED

    FOR TYPE 1 == TIDAL TAIL EJECTA
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
        super(Kasen1, self).__init__(**kwargs)

        # Read in times and frequencies arrays (same for all SEDs)

	open_path = "./kasen_seds/"

        kasen_frequencies = pickle.load( open( open_path + "frequency_angstroms.p" , "rb" ) )
        kasen_times = pickle.load( open( open_path + "times_days.p", "rb") )


    def process(self, **kwargs):
        open_path = "./kasen_seds/"

	# Physical parameters from Kasen simulations, provided by
        # neutrinosphere module (thank you, Jessica Metzger)
        self._vk = kwargs['vk']
        self._xlan = kwargs['xlan']
        self._mass = kwargs['mejecta']
        
        # mass fractional weight provided by neutrinoshere module
        self._mass_weight = kwargs['mass_weight']
        
        # viewing angle and opening angle 
        self._phi = kwargs['phi']
        self._theta = kwargs['theta']

        # Get times + other important things
        self._times = kwargs[self.key('dense_times')]
        self._band_indices = kwargs['all_band_indices']
        self._frequencies = kwargs['all_frequencies']

        # Total weight function
        # TYPE 1 == TIDAL TAIL SO GEOMETRIC FACTOR IS 1 - geometrical_weight
        weight_goem = 2*np.cos(self._theta)*(1. - np.cos(self._phi))
        weight = self._mass_weight * (1. - weight_goem)


        # Some temp vars for speed.
        cc = self.C_CONST

        zp1 = 1.0 + kwargs[self.key('redshift')]
        Azp1 = u.Angstrom.cgs.scale / zp1
        czp1 = cc / zp1


        seds = []
        rest_wavs_dict = {}

        # Find nearest neighbors to the Kasen-calculated simulation
        m_closest = MASS_S[(np.abs(MASS-self._mass)).argmin()]
        v_closest = VKIN_S[(np.abs(VKIN-self._vk)).argmin()]
        x_closest = XLAN_S[(np.abs(XLAN-self._xlan)).argmin()]

        # Open nearest neighbor file
        kasen_seds = pickle.load( open( save_path + 'knova_d1_n10_m' + m_closest + '_vk' + v_closest + '_fd1.0_Xlan' + x_closest + '.0.p', "rb" ) )

        # For each time
        for ti, t in enumerate(self._times):

            # Find index of closest time: this is the SED we will pull 
            t_closest_i = (np.abs(kasen_times-t)).argmin()

            # Create a rest-frame wavelength array (tbh I don't know how this
            # works, but I can venture a guess that this more or less right)
            bi = self._band_indices[ti]
            if bi >= 0:
                rest_wavs = rest_wavs_dict.setdefault(
                    bi, self._sample_wavelengths[bi] * Azp1)
            else:
                rest_wavs = np.array(  # noqa: F841
                    [czp1 / self._frequencies[ti]])

            # Evaluate the SED at the rest frame frequencies
            sed = np.array([])
            for w in rest_wavs:
                # find index of closest wav
                w_closest_i = np.abs(kasen_frequencies-w).argmin()

                sed = np.append(sed, weight * kasen_seds['SEDs'][t_closest_i][w_closest_i] )


        seds = self.add_to_existing_seds(seds, **kwargs)

        return {'sample_wavelengths': self._sample_wavelengths, 'seds': seds}
