"""Definitions for the `Neutrinosphere0` class."""
import numpy as np
from astrocats.catalog.source import SOURCE
from scipy import interpolate
from scipy.constants import c, G 
from mosfit.modules.energetics.energetic import Energetic
import os

# Important: Only define one ``Module`` class per file.


class Neutrinosphere0(Energetic):
    """Generate Kasen's `v_k`, `X_lan`, `m_ej`, and a mass-corresponding
        weight based on some other parameters.

        TYPE0 == KASEN0  == SHOCK HEATED EJECTA

        Actual physics done by Jessica Metzger <jmetzger@uchicago.edu>, 
        implemented by Kamile Lukosiute <klukosiu@wellesley.edu>

        August 2018
    """

    _REFERENCES = [
        {SOURCE.BIBCODE: '2014ApJ...789L..39W'},
        {SOURCE.BIBCODE: '2014MNRAS.443.3134P'},
        {SOURCE.BIBCODE: '2015MNRAS.452.3894G'},
        {SOURCE.BIBCODE: '1978A&A....67..185T'},
        {SOURCE.BIBCODE: '2015ApJ...808..188P'},
        {SOURCE.BIBCODE: '1996ApJ...471..331Q'},
        {SOURCE.BIBCODE: '2017Natur.551...80K'},
        {SOURCE.BIBCODE: '2015ApJ...815...82L'},
        {SOURCE.BIBCODE: '2018CQGra..35c4001M'}, 
        {SOURCE.BIBCODE: '2017PhRvL.119w1102S'}, 
        {SOURCE.BIBCODE: '2017PhRvD..96l4005B'}, 
        {SOURCE.BIBCODE: '2017MNRAS.472..904L'}, 
        {SOURCE.BIBCODE: '2017ApJ...846..114F'}, 
        {SOURCE.BIBCODE: '2016MNRAS.460.3255R'}, 
        {SOURCE.BIBCODE: '2015PhRvD..91f4059S'}, 
        {SOURCE.BIBCODE: '2015MNRAS.448..541J'}, 
        {SOURCE.BIBCODE: '2014MNRAS.443.3134P'}, 
        {SOURCE.BIBCODE: '2014ApJ...789L..39W'}
    ]

    C = c 		# m/s
    GG = G 		# Global variables are hard okay?  
    MSUN = 2e30
    M = 2.75*MSUN       # kg
    SB = 1.028e-5       # 10^51 erg/s, km^-2, MeV^-4
    M_N = 939.6         # Neutron mass, MeV
    M_P=  938.3         # Proton mass, MeV
    D = 1.293           # Neutron proton mass difference, MeV
    KCONST = 9.3e-48    #m^2 / MeV^2
    ENG = 6.24151e56    #MeV per 10^51 erg/s
    YE0  = .06          # initial electron fraction determined 
                        # from cold beta-equilibrium value for 
                        # neutron stars (2014ApJ...789L..39W)

    # NEUTRINOSPHERE CONSTANTS
    # Temperature of the anti/neutrinosphere
    # MeV
    # From scaled up data in 2014MNRAS.443.3134P
    T_NEUT_DISK = 4.3
    T_NEUT_NS = 7.

    T_ANTI_DISK = 5.4
    T_ANTI_NS = 8.7

    # Neutrinosphere geometry
    # km
    R_NEUT = 40.
    R_ANTI = 27.
    H = 15. 

    # Projected areas for each viewing angle
    # km**2
    A_NEUT_DISK_POL = np.pi*(R_NEUT**2)
    A_NEUT_DISK_EQU = 4*R_NEUT*H

    A_ANTI_DISK_POL = np.pi*(R_ANTI**2)
    A_ANTI_DISK_EQU = 4*R_ANTI*H

    A_NS = np.pi*(H**2)

    # Luminosities for each viewing angle
    # MeV
    L_NEUT_DISK_POL = SB * 4 * (A_NEUT_DISK_POL * T_NEUT_DISK**4) * ENG
    L_NEUT_NS_POL = SB * 4 * (A_NS * T_NEUT_NS**4) * ENG
    L_NEUT_EQU = SB * 4 * (A_NEUT_DISK_EQU * T_NEUT_DISK**4) * ENG
    L_NEUT_POL = L_NEUT_DISK_POL + L_NEUT_NS_POL

    L_ANTI_DISK_POL = SB * 4 * (A_ANTI_DISK_POL * T_ANTI_DISK**4) * ENG
    L_ANTI_NS_POL = SB * 4 * (A_NS * T_ANTI_NS**4) * ENG
    L_ANTI_EQU =  SB * 4 * (A_ANTI_DISK_EQU * T_ANTI_DISK**4) * ENG
    L_ANTI_POL = L_ANTI_DISK_POL + L_ANTI_NS_POL

    # Average anti/neutrino energies for each viewing angle
    # MeV
    # 2015MNRAS.452.3894G
    k = 3.1514
    E_NEUT_POL = ( k * T_NEUT_DISK * T_NEUT_NS * L_NEUT_POL / 
        (L_NEUT_DISK_POL * T_NEUT_NS + L_NEUT_NS_POL * T_NEUT_DISK) ) 
    E_NEUT_EQU = k * T_NEUT_DISK

    E_ANTI_POL = ( k * T_ANTI_DISK * T_ANTI_NS * L_ANTI_POL /
        (L_ANTI_DISK_POL * T_ANTI_NS + L_ANTI_NS_POL * T_ANTI_DISK))
    E_ANTI_EQU = k * T_ANTI_DISK


    # Factor to go from "transition velocity" to "kinetic velocity"
    # 1978A&A....67..185T
    V_FACTOR=((14/9.)*(1/4. + 1/5. - 1/(50*(10**.5))))**.5



    def __init__(self, **kwargs):
        """Initialize module."""
        super(Neutrinosphere0, self).__init__(**kwargs)

        self._dir_path = os.path.dirname(os.path.realpath(__file__))

        # Table giving lanthanide fraction for a given Y_e, s, Msph, and velocity
        self.LRDATA = open(os.path.join(self._dir_path , 'LRdata_modified.txt'),'r').readlines()[1:]
        self.LRDATA = np.array([x.rstrip().split(',') for x in self.LRDATA]).astype(np.float)

        # Exclude data that won't be needed
        self.LRDATA = self.LRDATA[np.where((self.LRDATA[:,1] >= 7.5) & (self.LRDATA[:,1] <= 32))]
        self.LR_MASSES = np.unique(self.LRDATA[:,2])



    def getYe(self, vf, t0, dt):
        '''
        Integrate Ye, electron fraction for dynamical ejecta
        
        COMPONENT 0 == SHOCK HEATED EJECTA

        vf = final velocity of ejecta [c]
        t0 = something???? [ms]
        dt = integration time step [s]
        '''
        # THIS SECTION SPECIFIC TO TYPE 0 EJECTA
        # ---------------------------------------------------------------------
        L_n = self.L_NEUT_POL
        E_n = self.E_NEUT_POL
        L_a = self.L_ANTI_POL
        E_a = self.E_ANTI_POL
        # ---------------------------------------------------------------------
        # anti/neutrino capture cross sections
        W_n = 1+1.02*1.2*1.14*E_n/self.M_N
        W_a = 1-7.22*1.2*1.14*E_a/self.M_P

        sigma_n = self.KCONST * (1.2*E_n**2 + 2*self.D*E_n + self.D**2)*W_n
        sigma_a = self.KCONST * (1.2*E_a**2 - 2*self.D*E_a + self.D**2)*W_a
        
        # anti/neutrino capture reaction rates
        k_n = L_n*sigma_n / (4*np.pi*E_n)
        k_a = L_a*sigma_a / (4*np.pi*E_a)
        k1 = k_n
        k2 = k_n + k_a
        
        # ejected from ~surface of neutron star
        r = 15.*1000  # m
        t = t0*0.001 # s
        v = ( (vf*self.C)**2. + 2.*self.GG*self.M / r)**0.5 # m/s - starting velocity

        reactions = 0
        for n in range(5000): # sum up reactions
            # anti/neutrino luminosity (and thus reaction rates) 
            # increases linearly through first 5ms
            if t < .005:
                k2_new = k2 * (max(0,t)/0.005)
            else: k2_new=k2
            
            #integrate forward reactions, a, v, r, t
            reactions += (k2_new / (r**2)) * dt
            a = -self.GG*self.M/(r**2)
            v += a*dt
            r += v*dt + .5*a*(dt**2)
            t += dt
            if r < 0: # if 0, function will call again w/smaller timestep
                return 0
        
        Ye_eq = k1/k2
        Ye_answer = Ye_eq-(Ye_eq-self.YE0)*(np.exp(-reactions))
        
        return Ye_answer

    def process(self, **kwargs):
        '''
        Starting with the parameters of ejecta mass, coasting velocity 
        (Kasen's transition velocity), and opening polar angle of shock
        heated ejecta, return equivalent spherical mass ejecta, kinetic
        velocity, lanthanide mass fraction, and SED mass-corresponding
        mass_weight.

        FREE PARAMETERS
        mejecta = ejecta mass [M_sun] 
        vcoast = coasting velocity [c]
        phi = opening angle [radians]

        RETURN PARAMETERS - DETERMINES SED
        Msph = spherical mass equivalent [M_sun]
        v_k = kinetic velocity [c]
        Xlan = lanthanide fraction [unitless?]

        COMPONENT 0 == SHOCK HEATED POLAR EJECTA

        '''
        # Free parameters
        self._mejecta = kwargs[self.key('mejecta')]
        self._vcoast  = kwargs[self.key('vcoast')]
        self._phi     = kwargs[self.key('phi')]

        #mass fractions (or volume fractions, same thing)
        shock_mass_fraction = ( 0.5 * ( (2+np.cos(self._phi) ) *
            ( 1 - np.cos(self._phi) )**2 + 
            ( np.sin(self._phi)**2 ) * np.cos(self._phi) ) )

        # THIS SECTION SPECIFIC TO TYPE 0 EJECTA
        # ---------------------------------------------------------------------
        entropy = 30.
        mass_weight = shock_mass_fraction
        Ye = self.getYe(self._vcoast,2,2.5e-6)
        if Ye == 0: #in case integration failed, repeat with smaller timestep
            Ye = self.getYe(self._vcoast,2,1.0e-6)
        # ---------------------------------------------------------------------

        Msph = self._mejecta / mass_weight
        v_k  = self._vcoast * self.V_FACTOR
        
        #cut out unneeded masses, Ye's from Xla table
        lower_mass = 0
        upper_mass = 1
        if Msph > 0.001:
            try: lower_mass = self.LR_MASSES[np.where(self.LR_MASSES < Msph)[0][-1]]
            except IndexError: lower_mass = 0
        if Msph < 0.1:
            try: upper_mass = self.LR_MASSES[np.where(self.LR_MASSES > Msph)[0][0]]
            except IndexError: upper_mass = 1

        lrdatatemp = self.LRDATA[np.where((self.LRDATA[:,0] >= Ye-.04) 
            & (self.LRDATA[:,0] <= Ye+.04) & (self.LRDATA[:,2] >= 
            lower_mass) & (self.LRDATA[:,2] <= upper_mass))]
        
        # Interpolate grid for Xla. If outside grid, it'll return Xla=0
	# Multiplying by 1. converts np.array of one val to float. Infuriating
	# but... fine. 
        Xla = interpolate.griddata(( lrdatatemp[:,0], 
            lrdatatemp[:,1],lrdatatemp[:,2],lrdatatemp[:,3]), 
            lrdatatemp[:,4], (Ye, entropy, Msph, self._vcoast), 
            method='linear',fill_value=0) * 1.

	
        #exclude ones out of range of Kasen's SED grid
        if Xla>.1: Xla=.1
        if Xla<1.0e-9: Xla=1.0e-9
	        
        #spherical mass [Msun], kinetic velocity [c], lanthanide mass fraction, SED mass_weight
        return {'Msph': Msph, 'vk': v_k, 'xlan': Xla, 'mass_weight': mass_weight}


