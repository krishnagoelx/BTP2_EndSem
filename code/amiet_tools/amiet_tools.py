"""
amiet_tools - a Python package for turbulence-aerofoil noise prediction.
https://github.com/fchirono/amiet_tools
Copyright (c) 2020, Fabio Casagrande Hirono

The 'amiet_tools' (AmT) Python package provides a reference implementation of
Amiet's [JSV 41, 1975] model for turbulence-aerofoil interaction noise with
extensions. These functions allow the calculation of the surface pressure jump
developed over the aerofoil surface (i.e. the acoustic source distribution) in
response to incoming turbulence, and of the acoustic field radiated by the
interaction.

Incoming turbulence can be a single sinusoidal gust, or a sum of incoherent
gusts with amplitudes given by a prescribed energy spectrum.


Dependencies:
    - numpy: array processing for numbers, strings, records, and objects;
    - scipy: scientific library.


All dependencies are already included in the Anaconda Python Distribution, a
free and open source distribution of Python. Anaconda 4.8.2 (with Python 3.7)
was used to develop and test AmT, and is recommended for using AmT.


Author:
    Fabio Casagrande Hirono - fchirono@gmail.com


Main Technical References:

    Amiet, R. K., "Acoustic radiation from an airfoil in a turbulent stream",
    Journal of Sound and Vibration, Vol. 41, No. 4:407–420, 1975.

    Blandeau, V., "Aerodynamic Broadband Noise from Contra-Rotating Open
    Rotors", PhD Thesis, Institute of Sound and Vibration Research, University
    of Southampton, Southampton - UK, 2011.

    Casagrande Hirono, F., "Far-Field Microphone Array Techniques for Acoustic
    Characterisation of Aerofoils", PhD Thesis, Institute of Sound and
    Vibration Research, University of Southampton, Southampton - UK, 2018.

    Reboul, G., "Modélisation du bruit à large bande de soufflante de
    turboréacteur", PhD Thesis, Laboratoire de Mécanique des Fluides et
    d’Acoustique - École Centrale de Lyon, Lyon - France, 2010.

    Roger, M., "Broadband noise from lifting surfaces: Analytical modeling and
    experimental validation". In Roberto Camussi, editor, "Noise Sources in
    Turbulent Shear Flows: Fundamentals and Applications". Springer-Verlag,
    2013.

    de Santana, L., "Semi-analytical methodologies for airfoil noise
    prediction", PhD Thesis, Faculty of Engineering Sciences - Katholieke
    Universiteit Leuven, Leuven, Belgium, 2015.
"""

import numpy as np
import scipy.special as ss
import scipy.optimize as so     # for shear layer correction functions
import scipy.integrate as integrate
#import mpmath as mp

# import gc


class TestSetup:
    """
    Class to store test setup variables. Initializes to DARP2016 configuration
    by default.
    """

    def __init__(self):
        # Acoustic characteristics
        self.c0 = 340.	# Speed of sound [m/s]
        self.rho0 =1.2	# Air density [kg/m**3]
        self.p_ref = 20e-6	# Ref acoustic pressure [Pa RMS]

        # turbulent flow properties
        self.Ux = 60			# flow velocity [m/s]
        self.turb_intensity = 0.025	# turbulence intensity = u_rms/U
        self.length_scale = 0.007	# turb length scale [m]

        # shear layer height
        self.z_sl = -0.075	# [m]

        self.flow_dir = 'x'              # mean flow in the +x dir
        self.dipole_axis = 'z'           # airfoil dipoles are pointing 'up' (+z dir)

        self._calc_secondary_vars()


    def _calc_secondary_vars(self):
        # calculate other setup variables from initial ones
        self.Mach = self.Ux/self.c0
        self.beta = np.sqrt(1-self.Mach**2)
        self.flow_param = (self.flow_dir, self.Mach)

    def export_values(self):
        return (self.c0, self.rho0, self.p_ref, self.Ux, self.turb_intensity,
                self.length_scale, self.z_sl, self.Mach, self.beta,
                self.flow_param, self.dipole_axis)


class AirfoilGeom:
    """
    Class to store aerofoil geometry
    """

    def __init__(self):
        # Aerofoil geometry
        self.b = 0.075	# airfoil half chord [m]
        self.d = 0.225	# airfoil half span [m]
        self.Nx = 100	# number of chordwise points (non-uniform sampl)
        self.Ny = 101	# number of spanwise points (uniform sampl)

        self._calc_grid()


    def _calc_grid(self):
        self.XYZ, self.dx, self.dy = create_airf_mesh(self.b, self.d, self.Nx,
                                                      self.Ny)

    def export_values(self):
        return (self.b, self.d, self.Nx, self.Ny, self.XYZ, self.dx,
                self.dy)


class FrequencyVars:
    """
    Class to store frequency-related variables
    """
    def __init__(self, f0, testSetup):

        self.f0 = f0                # frequency [Hz]

        c0 = testSetup.c0           # speed of sound [m/s]
        Mach = testSetup.Mach       # Mach number
        beta = testSetup.beta

        self.k0 = 2*np.pi*f0/c0     # acoustic wavenumber
        self.Kx =self.k0/Mach       # gust/hydrodynamic chordwise wavenumber

        # gust/hydrodynamic spanwise critical wavenumber
        self.Ky_crit = self.Kx*Mach/beta

    def export_values(self):
        return (self.k0, self.Kx, self.Ky_crit)


# %% *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# load test setup from file

def loadTestSetup(*args):
    """
    Load test variable values for calculations, either from default values or
    from a given .txt test configuration setup file. File must contain the
    following variable values, in order:

        c0                  speed of sound [m/s]
        rho0                density of air [kg/m**3]
        p_ref               reference pressure [Pa RMS]
        Ux                  mean flow velocity [m/s]
        turb_intensity      turbulent flow intensity [ = u_rms/Ux]
        length_scale        turbulence integral length scale [m]
        z_sl                shear layer height (from aerofoil)

    Empty lines and 'comments' (starting with '#') are ignored.

    path_to_file: str [Optional]
        Relative path to setup file
    """

    # if called without path to file, load default testSetup (DARP2016)
    if len(args)==0:
        return TestSetup()

    else:
        # initialize new instance of testSetup
        testSetupFromFile = TestSetup()

        varList = ['c0', 'rho0', 'p_ref', 'Ux', 'turb_intensity',
                   'length_scale', 'z_sl']
        i=0

        # open file name given by 'args[0]'
        with open(args[0]) as f:
            # get list with file lines as strings
            all_lines = f.readlines()

            # for each line...
            for line in all_lines:

                # skip comments and empty lines
                if line[0] in ['#', '\n']:
                    pass

                else:
                    words = line.split('\t')
                    # take 1st element as value (ignore comments)
                    exec('testSetupFromFile.' + varList[i] + '=' + words[0])
                    i+=1

        # calculate other variables from previous ones
        testSetupFromFile._calc_secondary_vars()

        return testSetupFromFile



def loadAirfoilGeom(*args):
    """
    Load airfoil geometry values for calculations, either from default values or
    from a given .txt airfoil geometry file. File must contain the
    following variable values, in order:

        b       Airfoil semichord [m]
        d       Airfoil semispan [m]
        Nx      Number of chordwise points (non-uniform sampling)
        Ny      Number of spanwise points (uniform sampling)

    Empty lines and 'comments' (starting with '#') are ignored.

    path_to_file: str [Optional]
        Relative path to airfoil geometry file
    """

    # if called without path to file, load default geometry (DARP2016)
    if len(args)==0:
        return AirfoilGeom()

    else:
        # initialize new instance of testSetup
        airfoilGeomFromFile = AirfoilGeom()

        varList = ['b', 'd', 'Nx', 'Ny']
        i=0

        #path_to_file = '../DARP2016_AirfoilGeom.txt'
        with open(args[0]) as f:
            # get list with file lines as strings
            all_lines = f.readlines()

            # for each line...
            for line in all_lines:

                # skip comments and empty lines
                if line[0] in ['#', '\n']:
                    pass

                else:
                    words = line.split('\t')
                    # take 1st element as value (ignore comments)
                    exec('airfoilGeomFromFile.' + varList[i] + '=' + words[0])
                    i+=1

        # calculate other variables from previous ones
        airfoilGeomFromFile._calc_grid()

        return airfoilGeomFromFile



def DARP2016_MicArray():
    """
    Returns the microphone coordinates for the 'near-field' planar microphone
    array used in the 2016 experiments with the DARP open-jet wind tunnel at
    the ISVR, Univ. of Southampton, UK [Casagrande Hirono, 2018].

    The array is a Underbrink multi-arm spiral array with 36 electret
    microphones, with 7 spiral arms containing 5 mics each and one central mic.
    The array has a circular aperture (diameter) of 0.5 m.

    The calibration factors were obtained using a 1 kHz, 1 Pa RMS calibrator.
    Multiply the raw mic data by its corresponding factor to obtain a
    calibrated signal in Pascals.
    """

    M = 36                      # number of mics
    array_height = -0.49        # [m] (ref. to airfoil height at z=0)

    # mic coordinates (corrected for DARP2016 configuration)
    XYZ_array = np.array([[  0.     ,  0.025  ,  0.08477,  0.12044,  0.18311,  0.19394,
                             0.01559,  0.08549,  0.16173,  0.19659,  0.24426, -0.00556,
                             0.02184,  0.08124,  0.06203,  0.11065, -0.02252, -0.05825,
                            -0.06043, -0.11924, -0.10628, -0.02252, -0.09449, -0.15659,
                            -0.21072, -0.24318, -0.00556, -0.05957, -0.13484, -0.14352,
                            -0.19696,  0.01559,  0.02021, -0.01155,  0.03174, -0.00242],
                           [-0.     , -0.     ,  0.04175,  0.11082,  0.10542,  0.15776,
                            -0.01955, -0.04024, -0.02507, -0.07743, -0.05327, -0.02437,
                            -0.09193, -0.14208, -0.20198, -0.22418, -0.01085, -0.0744 ,
                            -0.1521 , -0.17443, -0.22628,  0.01085, -0.00084, -0.04759,
                            -0.01553, -0.05799,  0.02437,  0.07335,  0.09276,  0.15506,
                             0.15397,  0.01955,  0.09231,  0.16326,  0.20889,  0.24999],
                          array_height*np.ones(M)])

    # calibration factors
    array_cal = np.array([73.92182641429085,    96.84446743391487,  85.48777846463159,
                          85.24410968090712,    83.63917149322562,  68.94090765134432,
                          79.2385037527723,     112.77357210746612, 84.8483307868491,
                          87.18956628936178,    97.75046920293282,  89.2829545690508,
                          79.51644155562396,    90.39403884030057,  80.71754629014218,
                          89.4418210091059,     98.33634233056068,  79.2212022850229,
                          91.25543447201031,    89.55040012572815,  85.77495667666254,
                          82.74418222820202,    84.63061055646973,  77.01568014644964,
                          95.52764533324982,    92.16734812591154,  95.27123074600838,
                          87.93335310521428,    96.65066131188675,  93.58564782091074,
                          78.1446818728945,     101.3047738767648,  83.68569643491034,
                          84.7981031520437,     94.40796508430756,  83.52266614867919])

    return XYZ_array, array_cal


# %% *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# --->>> wrappers to create surface pressure cross-spectral matrix

def calc_airfoil_Sqq(testSetup, airfoilGeom, frequencyVars, Ky_vec, Phi):
    """
    Calculates the aerofoil surface pressure jump cross-spectral density matrix
    (CSM).

    Parameters
    ----------
    testSetup : instance of TestSetup class
        Instance containing test setup data.

    airfoilGeom : instance of AirfoilGeom class
        Instance containing airfoil geometry data.

    frequencyVars : instance of FrequencyVars class
        Instance containing frequency data.

    Ky_vec : (N_ky,) array_like
        1D array containing range of negative and positive spanwise gust
        wavenumber values.

    Phi : (N_ky,) array_like
        1D array containing the turbulent wavenumber energy spectrum for values
        of 'Ky_vec'.

    Returns
    -------
    Sqq : (Nx*Ny, Nx*Ny) array_like
        2D array containing the cross-spectral density of the airfoil surface
        pressure.

    Sqq_dxy : (Nx*Ny, Nx*Ny) array_like
        2D array containing the surface area weights "(dx*dy)*(dx' * dy')" for
        all pairs (x, y), (x', y') of airfoil points.

    Notes
    -----
    'Ky_vec' can be calculated using 'amiet_tools.ky_vec' function.

    'Phi' can be calculated with 'amiet_tools.Phi_2D' function.

    For surface pressure analyses, such as cross-spectrum phase or coherence
    lengths, use 'Sqq' only. For acoustic radiation analysis, the function
    'amiet_tools.calc_radiated_Spp' applies the equivalent areas 'Sqq_dxy' to
    numerically calculate the integration over the airfoil surface.
    """

    # export variables to current namespace
    (c0, rho0, p_ref, Ux, turb_intensity, length_scale, z_sl, Mach, beta,
     flow_param, dipole_axis) = testSetup.export_values()
    (b, d, Nx, Ny, XYZ_airfoil, dx, dy) = airfoilGeom.export_values()
    (k0, Kx, Ky_crit) = frequencyVars.export_values()

    # Surface area weighting matrix for applying to Sqq
    dxy = np.ones((Ny, Nx))*dx[np.newaxis, :]*dy
    Sqq_dxy = np.outer(dxy, dxy)

    # gust chordwise wavenumber
    Kx = k0/Mach

    # gust spanwise wavenumber interval
    dky = Ky_vec[1]-Ky_vec[0]

    # source CSM
    Sqq = np.zeros((Nx*Ny, Nx*Ny), 'complex')

    for kyi in range(Ky_vec.shape[0]):

        # sinusoidal gust peak value
        w0 = np.sqrt(Phi[kyi])

        # Pressure 'jump' over the airfoil (for single gust)
        delta_p1 = delta_p(rho0, b, w0, Ux, Kx, Ky_vec[kyi], XYZ_airfoil[0:2],
                           Mach)

        # reshape for vector calculation
        delta_p1_calc = delta_p1.reshape(Nx*Ny)

        Sqq[:, :] += np.outer(delta_p1_calc, delta_p1_calc.conj())

    Sqq *= (Ux*dky)

    return Sqq, Sqq_dxy


def calc_radiated_Spp(testSetup, airfoilGeom, frequencyVars, Ky_vec, Phi, G):
    """
    Calculates the cross-spectral density matrix (CSM) of the acoustic field
    radiated by the airfoil.

    Parameters
    ----------
    testSetup : instance of TestSetup class
        Instance containing test setup data.

    airfoilGeom : instance of AirfoilGeom class
        Instance containing airfoil geometry data.

    frequencyVars : instance of FrequencyVars class
        Instance containing frequency data.

    Ky_vec : (N_ky,) array_like
        1D array containing range of negative and positive spanwise gust
        wavenumber values. Calculate with 'amiet_tools.ky_vec' function.

    Phi : (N_ky,) array_like
        1D array containing the turbulent wavenumber energy spectrum for values
        of 'Ky_vec'. Calculate with 'amiet_tools.Phi_2D' function.

    G : (M, Nx*Ny) array_like
        2D matrix of complex transfer function between 'M' observer locations
        and 'Nx*Ny' points over the airfoil surface.

    Returns
    -------
    Spp : (M, M) array_like
        2D matrix of cross-spectral density of the acoustic field seen at 'M'
        observers/microphones.

    Notes
    -----
    'G' can be calculated using 'amiet_tools.dipole3D' for dipole sources in
    steady or moving medium, or 'amiet_tools.dipole_shear' for dipole sources
    in a shear layer medium.

    """

    Sqq, Sqq_dxy = calc_airfoil_Sqq(testSetup, airfoilGeom, frequencyVars, Ky_vec, Phi)

    # apply weights for surface area
    Sqq *= Sqq_dxy

    return (G @ Sqq @ G.conj().T)*4*np.pi


# %% *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# --->>> Discretization of aerofoil surface

def chord_sampling(b, Nx=200, exp_length=2):
    """
    Calculates 'Nx' points non-uniformly sampled over the half-open
    interval (-b, b], with the leading edge at '-b'.

    Parameters
    ----------
    b : float
        Airfoil semichord, in meters.

    Nx : int, optional
        Number of points used in sampling.

    exp_length : float, optional
        Length of exponential interval to be sampled.

    Returns
    -------
    x : (Nx,) array_like
        1D array of non-uniformly sampled points in half-open interval (-b, +b]

    dx : (Nx,) array_like
        1D vector of non-uniform sample intervals.

    Notes
    -----
    The function samples an exponential function over the interval
    [-exp_length/2, exp_length/2] uniformly at 'Nx' points, which are remapped
    to the actual chord interval [-b, b] and the leading edge then removed.
    This type of sampling assigns more points around the leading edge.

    Higher values of 'exp_length' provide more non-uniform sampling, while
    lower values provide more uniform sampling.
    """

    # calculates exponential curve
    x = np.exp(np.linspace(-exp_length/2, exp_length/2, Nx+1))

    # normalise to [0, 1] interval
    x = x-x.min()
    x = x/x.max()

    # normalise to [-b, b] interval
    x = (x*2*b)-b

    # calculate dx (for numerically integrating chord functions)
    dx = np.diff(x)

    # remove leading edge singularity
    x = x[1:]

    return x, dx


def create_airf_mesh(b, d, Nx=100, Ny=101):
    """
    Creates a (3, Ny, Nx) mesh containing the airfoil surface points coordinates.

    Parameters
    ----------
    b : float
        Airfoil semichord, in meters.

    d : float
        Airfoil semispan, in meters.

    Nx : int, optional
        Number of points used in sampling the chord (non-uniformly).

    Ny : int, optional
        Number of points used in sampling the span (uniformly).

    Returns
    -------
    XYZ_airfoil : (3, Ny, Nx) array_like
        3D array containing the coordinates of each point on the sampled
        airfoil surface.

    dx : (Nx,) array_like
        1D vector of non-uniform chord sample intervals.

    dy : float
       Span sample interval.

    Notes
    -----
    The airfoil 'z' coordinate is always set to zero.
    """

    x_airfoil, dx = chord_sampling(b, Nx)

    y_airfoil, dy = np.linspace(-d, d, Ny, retstep=True)
    #dy = y_airfoil[1] - y_airfoil[0]

    XY_airfoil = np.meshgrid(x_airfoil, y_airfoil)
    Z_airfoil = np.zeros(XY_airfoil[0].shape)
    XYZ_airfoil = np.array([XY_airfoil[0], XY_airfoil[1], Z_airfoil])

    return XYZ_airfoil, dx, dy


# %% *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# --->>> Sound propagation functions

def dipole3D(xyz_source, xyz_obs, k0, dipole_axis='z', flow_param=None,
             far_field=False):
    """
    Calculates a (M, N)-shaped matrix of dipole transfer functions
    between 'N' sources and 'M' observers/microphones at a single frequency.

    Parameters
    ----------
    xyz_source : (3, N) array_like
        Array containing the (x, y, z) coordinates of the 'N' sources.

    xyz_obs : (3, M) array_like
        Array containing the (x, y, z) coordinates of the 'M' observers /
        microphones.

    k0 : float
        Acoustic wavenumber / spatial frequency. Can be obtained from the
        temporal frequency 'f' [in Hz] and the speed of sound 'c0' [in m/s]
        as 'k0 = 2*pi*f/c0'.

    dipole_axis : {'x', 'y', 'z'}
        String indicating the direction aligned with the dipole axis; default
        is 'z'.

    flow_param : (flow_dir = {'x', 'y', 'z'}, Mach = float), optional
        Tuple containing the flow parameters: 'flow_dir' is the flow direction
        and can assume the strings 'x', 'y' or 'z'; 'Mach' is the Mach number
        of the flow, calculated by 'Mach = U/c0', where 'U' is the flow
        velocity [in m/s].

    far_field : boolean, optional
        Boolean variable determining if returns far-field approximation;
        default is 'False'

    Returns
    -------
    G : (M, N) array_like
        Matrix of complex transfer functions between sources and observers.

    Notes
    -----
    To calculate a vector of acoustic pressures 'p' at the observer locations
    as induced by a vector of source strengths 'q' at the source locations,
    simply do

    >>> p = G @ q

    """

    # check if source/observer coordinates have appropriate shape
    # ( i.e. [3, N_points])
    if (xyz_source.ndim == 1):
        xyz_source = np.array([xyz_source]).transpose()

    if (xyz_obs.ndim == 1):
        xyz_obs = np.array([xyz_obs]).transpose()

    # M = xyz_obs.shape[1]
    # N = xyz_source.shape[1]

    # if calculating Greens function for a steady medium (no flow):
    if not flow_param:

        # matrix of mic-source distances (Euclidean distance)
        r = np.sqrt(((xyz_obs[:, :, np.newaxis]
                      - xyz_source[:, np.newaxis, :])**2).sum(0))

        # dictionary with direction-to-axis mapping
        dir_keys = {'x': 0, 'y': 1, 'z': 2}

        # read dictionary entry for 'flow_dir'
        i_dip = dir_keys[dipole_axis]

        # dipole directivity term (cosine)
        dip_cos = (xyz_obs[i_dip, :, np.newaxis]-xyz_source[i_dip])/r

        # matrix of Green's functions
        G_dipole = (1j*k0 + 1/r)*dip_cos*np.exp(-1j*k0*r)/(4*np.pi*r)

        return G_dipole

    # if calculating Greens function for a convected medium :
    else:

        # parse 'flow_param' tuple
        flow_dir, Mach = flow_param

        # flow correction factor
        beta = np.sqrt(1.-Mach**2)

        # dictionary with direction-to-axis mapping
        dir_keys = {'x': 0, 'y': 1, 'z': 2}

        # read dictionary entry for 'flow_dir'
        i_flow = dir_keys[flow_dir]

        # apply '1/beta' factor to all coordinates
        xyz_obsB = xyz_obs/beta
        xyz_sourceB = xyz_source/beta

        # apply extra '1/beta' factor to flow direction
        xyz_obsB[i_flow] = xyz_obsB[i_flow]/beta
        xyz_sourceB[i_flow] = xyz_sourceB[i_flow]/beta

        # matrix of delta-x in flow direction (uses broadcasting)
        xB = xyz_obsB[i_flow, :, np.newaxis] - xyz_sourceB[i_flow, :]

        # matrix of transformed mic-source distances (Euclidean distance)
        rB = np.sqrt(((xyz_obsB[:, :, np.newaxis]
                       - xyz_sourceB[:, np.newaxis, :])**2).sum(0))

        # read dictionary entry for 'flow_dir'
        i_dip = dir_keys[dipole_axis]

        # dipole directivity term (cosine)
        dip_cos = (xyz_obsB[i_dip, :, np.newaxis]-xyz_sourceB[i_dip])/(rB*beta)

        # if using far-field approximation...
        if far_field:
            # matrix of convected far-field greens functions
            sigma = np.sqrt(xyz_obs[0, :, np.newaxis]**2
                            + (beta**2)*(xyz_obs[1, :, np.newaxis]**2
                                         + xyz_obs[2, :, np.newaxis]**2))

            G_dipole = ((1j*k0*xyz_obs[2, :, np.newaxis]/(4*np.pi*(sigma**2)))
                        * np.exp(1j*k0*(Mach*xyz_obs[0, :, np.newaxis]-sigma)
                                 / (beta**2))
                        * np.exp(1j*k0*(xyz_obs[0, :, np.newaxis]-Mach*sigma)
                                 * xyz_source[0]/(sigma*(beta**2)))
                        * np.exp(1j*k0*(xyz_obs[1, :, np.newaxis]/sigma)
                                 *xyz_source[1]))

        else:
            # Matrix of convected Green's functions
            G_dipole = ((1j*k0 + 1./rB)*dip_cos*np.exp(1j*k0*Mach*xB)
                        * np.exp(-1j*k0*rB)/(4*np.pi*(beta**2)*rB))

        return G_dipole


def dipole_shear(XYZ_source, XYZ_obs, XYZ_sl, T_sl, k0, c0, Mach):
    """
    Calculates a (M, N)-shaped matrix of dipole transfer functions
    between 'N' sources and 'M' observers/microphones at a single frequency
    including shear layer refraction effects. The mean flow is assumed to be
    in the '+x' direction with velocity 'Ux = Mach*c0', and the dipoles are
    assumed to be located with their axes in the '+z' direction.

    Parameters
    ----------
    XYZ_source : (3, N) array_like
        Array containing the (x, y, z) coordinates of the 'N' sources.

    XYZ_obs : (3, M) array_like
        Array containing the (x, y, z) coordinates of the 'M' observers /
        microphones.

    XYZ_sl : (3, M, N) array_like
        Array containing the (x, y, z) coordinates of the shear-layer crossing
        point for an acoustic ray leaving the 'n'-th source and reaching
        the 'm'-th observer.

    T_sl : (3, M, N) array_like
        Array containing the total propagation time for an acoustic ray leaving
        the 'n'-th source and reaching the 'm'-th observer.

    k0 : float
        Acoustic wavenumber / spatial frequency. Can be obtained from the speed
        of sound 'c0' [m/s] and the angular frequency 'omega0' [rad/s] (or
        temporal frequency 'f0' [in Hz]) as 'k0 = 2*np.pi*f0/c0 = omega0/c0'.

    c0 : float
        Speed of sound in free air [m/s].

    Mach : float
        Mean flow Mach number; the flow velocity is 'Ux = Mach*c0'.

    Returns
    -------
    G : (M, N) array_like
        Matrix of complex transfer functions between sources and observers,
        including shear layer refraction effects.

    Notes
    -----
    This code uses a ray-acoustics approximation to predict the refraction
    effects at the shear layer crossing: for every source-observer pair, it
    shoots an acoustic ray that leaves the source, gets refracted at the shear
    layer and reaches the observer. The ray obeys the convected wave equation
    within the convecting region (i.e. before the shear layer), and the
    standard wave equation outside the convecting region.

    The 'XYZ_sl' and the 'T_sl' variables can be obtained from the
    'amiet_tools.ShearLayer_matrix' function.

    To calculate a vector of acoustic pressures 'p' at the observer locations
    as induced by a vector of source strengths 'q' at the source locations,
    simply do

    >>> p = G @ q

    """

    # check if source/observer coordinates have appropriate shape
    # ( i.e. [3, N_points])
    if (XYZ_source.ndim == 1):
        XYZ_source = np.array([XYZ_source]).transpose()

    if (XYZ_obs.ndim == 1):
        XYZ_obs = np.array([XYZ_obs]).transpose()

    # flow-corrected source-to-shear-layer propag distance
    sigma_sl = _sigma(XYZ_sl - XYZ_source[:, np.newaxis, :], Mach)

    # dipole in-flow (cosine) directivity2
    dip_cos = (XYZ_sl[2]-XYZ_source[2, np.newaxis, :])/sigma_sl

    beta2 = 1-Mach**2
    rbar_sl = sigma_sl/beta2

    # geometric shear-layer-to-mic distance
    r_lm = r(XYZ_sl - XYZ_obs[:, :, np.newaxis])

    omega0 = k0*c0

    G_dip_sl = ((1j*k0 + 1/(rbar_sl + r_lm))*dip_cos
                * np.exp(-1j*omega0*T_sl)/(4*np.pi*(sigma_sl + r_lm)))

    return G_dip_sl


# %% *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# --->>> Shear Layer Correction Functions

def r(x):
    """
    Cartesian radius of a point 'x' in 3D space

    Parameters
    ----------
    x : (3,) array_like
        1D vector containing the (x, y, z) coordinates of a point.

    Returns
    -------
    r : float
        Radius of point 'x' relative to origin of coordinate system
    """
    return np.sqrt(x[0]**2 + x[1]**2 + x[2]**2)


def r_bar(x, Mach):
    """
    Flow-transformed cartesian radius of a point 'x' in 3D space, with mean
    flow in the '+x' direction.

    Parameters
    ----------
    x : (3,) array_like
        1D vector containing the (x, y, z) coordinates of a point.

    Mach : float
        Mach number of the mean flow (assumed in '+x' direction).

    Returns
    -------
    r_bar : float
        Flow-transformed radius of point 'x' relative to origin of coordinate
        system.

    Notes
    -----
    This function considers full flow-transformed coordinates [Chapman, JSV
    233, 2000] to calculate the output.
    """
    beta = np.sqrt(1-Mach**2)
    return np.sqrt((x[0]/beta**2)**2 + (x[1]/beta)**2 + (x[2]/beta)**2)


def _sigma(x, Mach):
    """
    Flow-corrected cartesian radius of a point 'x' in 3D space, with mean flow
    in the '+x' direction.

    Parameters
    ----------
    x : (3,) array_like
        1D vector containing the (x, y, z) coordinates of a point.

    Mach : float
        Mach number of the mean flow (assumed in '+x' direction).

    Returns
    -------
    sigma : float
        Flow-corrected radius of point 'x' relative to origin of coordinate
        system

    Notes
    -----
    This function does NOT considers full flow-transformed coordinates
    [Chapman, JSV 233, 2000].
    """
    beta2 = 1-Mach**2
    return np.sqrt(x[0]**2 + beta2*(x[1]**2 + x[2]**2))


def t_sound(x1, x2, c0):
    """
    Calculates the time taken for sound to move from points 'x1' to 'x2' at
    speed 'c0'.

    Parameters
    ----------
    x1 : (3,) array_like
        1D vector containing the (x, y, z) coordinates of initial point.

    x2 : (3,) array_like
        1D vector containing the (x, y, z) coordinates of final point.

    c0 : float
        Speed of sound.

    Returns
    -------
    t : float
        Time taken for sound to travel from point 'x1' to point 'x2' at speed
        'c0'.
    """
    return r(x1-x2)/c0


def t_convect(x1, x2, Ux, c0):
    """
    Propagation time for sound in convected medium, with mean flow in the '+x'
    direction.

    Parameters
    -----
    x1 : (3,) array_like
        1D vector containing the (x, y, z) coordinates of initial point.

    x2 : (3,) array_like
        1D vector containing the (x, y, z) coordinates of final point.

    Ux : float
        Mean flow velocity, assumed in '+x' direction.

    c0 : float
        Speed of sound.

    Returns
    -------
    t_convect : float
        Time taken for sound to travel from point 'x1' to point 'x2' at speed
        'c0' while subject to convection effects in the '+x' direction.
    """

    Mach = Ux/c0
    beta2 = 1-Mach**2  # beta squared

    return (-(x2[0]-x1[0])*Mach + _sigma(x2-x1, Mach))/(c0*beta2)


def t_total(x_layer, x_source, x_mic, Ux, c0):
    """
    Total propagation time for sound to move from source to mic and through a shear layer.

    Parameters
    ----------
    x_layer : (3,) array_like
        1D vector containing the (x, y, z) coordinates of shear layer crossing
        point.

    x_source : (3,) array_like
        1D vector containing the (x, y, z) coordinates of final point.

    x_mic : (3,) array_like
        1D vector containing the (x, y, z) coordinates of final point.

    Ux : float
        Mean flow velocity, assumed in '+x' direction.

    c0 : float
        Speed of sound.

    Returns
    -------
    t_total : float
        Time taken for sound to travel from source, through shear layer, to mic.

    Notes
    -----
    The shear layer crossing point 'x_layer' can be obtained from
    'amiet_tools.ShearLayer_X' function.
    """
    return t_convect(x_source, x_layer, Ux, c0) + t_sound(x_layer, x_mic, c0)


def constr_xl(XYZ_sl, XYZ_s, XYZ_m, Ux, c0):
    """
    Shear layer solver constraint in 'xl'
    """
    Mach = Ux/c0
    beta2 = 1-Mach**2

    return np.abs((XYZ_sl[0]-XYZ_s[0])/_sigma(XYZ_sl-XYZ_s, Mach) - Mach
                  - beta2*(XYZ_m[0]-XYZ_sl[0])/r(XYZ_m-XYZ_sl))


def constr_yl(XYZ_sl, XYZ_s, XYZ_m, Ux, c0):
    """
    Shear layer solver constraint in 'yl'
    """
    Mach = Ux/c0

    return np.abs((XYZ_sl[1]-XYZ_s[1])/_sigma(XYZ_sl-XYZ_s, Mach)
                  - (XYZ_m[1]-XYZ_sl[1])/r(XYZ_m-XYZ_sl))


def ShearLayer_X(XYZ_s, XYZ_m, Ux, c0, z_sl):
    """
    Calculates the shear layer crossing point of an acoustic ray emitted from a
    source at 'xs' to a microphone at 'xm'.

    Parameters
    ----------
    XYZ_s : (3,) array_like
        1D vector containing the (x, y, z) coordinates of source (within mean flow)

    XYZ_m : (3,) array_like
        1D vector containing the (x, y, z) coordinates of microphone (outside mean flow)

    Ux : float
        Mean flow velocity, assumed in '+x' direction.

    c0 : float
        Speed of sound.

    z_sl : float
        Shear layer height.

    Returns
    -------
    XYZ_sl : (3,) array_like
        1D vector containing the (x, y, z) coordinates of the shear layer
        crossing point.

    Notes
    -----
    This code uses a numerical minimization routine to calculate the shear
    layer crossing point 'XYZ_sl' that minimizes the total travel time that an
    acoustic ray takes to propagate from a source point 'XYZ_s' within the
    mean flow to a microphone 'XYZ_m' outside the mean flow.
    """

    # optimization constraints
    cons = ({'type': 'eq', 'fun': lambda XYZ_sl: XYZ_sl[2]-z_sl},
            {'type': 'eq', 'fun':
             lambda XYZ_sl: constr_xl(XYZ_sl, XYZ_s, XYZ_m, Ux, c0)},
            {'type': 'eq', 'fun':
             lambda XYZ_sl: constr_yl(XYZ_sl, XYZ_s, XYZ_m, Ux, c0)})

    # initial guess (straight line between source and mic)
    XYZ_0 = ((XYZ_m + XYZ_s)*((z_sl-XYZ_s[2])/(XYZ_m[2] - XYZ_s[2])))

    # optimize and get result
    XYZ_sl_opt = so.minimize(t_total, XYZ_0, args=(XYZ_s, XYZ_m, Ux, c0),
                             method='SLSQP', constraints=cons)
    XYZ_sl = XYZ_sl_opt.x

    return XYZ_sl


def ShearLayer_matrix(XYZ_s, XYZ_o, z_sl, Ux, c0):
    """ Returns two matrices containing the propagation times and the shear
    layer crossing points for each source-observer pair.

    Parameters
    ----------
    XYZ_s : (3, N) array_like
        2D vector containing the (x, y, z) coordinates of 'N' source points
        (within mean flow)

    XYZ_o : (3, M) array_like
        2D vector containing the (x, y, z) coordinates of 'M' observer points
        (outside the mean flow)

    z_sl : float
        Shear layer height.

    Ux : float
        Mean flow velocity, assumed in '+x' direction.

    c0 : float
        Speed of sound.

    Returns
    -------
    T : (M, N) array_like
        2D array containing the total propagation times of the acoustic signal
        from each source 'n' to each observer 'm'.

    XYZ_sl : (3, M, N) array_like
        3D array containing the (x, y, z) coordinates of the shear layer
        crossing point for each source 'n' to each observer 'm'.

    Notes
    -----
    This is a convenience function, essentially two nested 'for'-loops around
    the 'ShearLayer_X' and 'ShearLayer_t' functions.
    """

    # ensure source/observer coordinates have appropriate shape
    # ( i.e. [3, N_points])
    if (XYZ_s.ndim == 1):
        XYZ_s = np.array([XYZ_s]).transpose()

    if (XYZ_o.ndim == 1):
        XYZ_o = np.array([XYZ_o]).transpose()

    # check if shear layer is located between sources and obs
    sl_height_error = "Shear layer is not located between all sources and observers"
    assert (np.prod(np.sign(XYZ_o[2] - z_sl))
            == np.prod(np.sign(z_sl - XYZ_s[2]))), sl_height_error

    # number of obs and sources
    M = XYZ_o.shape[1]
    N = XYZ_s.shape[1]

    XYZ_sl = np.zeros((3, M, N))     # shear layer crossing point
    T = np.zeros((M, N))            # propag time

    for n in range(N):
        for m in range(M):
            # shear layer crossing point
            XYZ_sl[:, m, n] = ShearLayer_X(XYZ_s[:, n], XYZ_o[:, m], Ux, c0,
                                           z_sl)
            # total propag time
            T[m, n] = t_total(XYZ_sl[:, m, n], XYZ_s[:, n],
                              XYZ_o[:, m], Ux, c0)

    return T, XYZ_sl


def ShearLayer_Corr(XYZ_s, XYZ_sl, XYZ_m, Ux, c0):
    """
    Calculates a corrected position for a microphone measuring sound through
    a shear layer using Amiet's method [JSV 58, 1978].

    Parameters
    ----------
    x_layer : (3,) array_like
        1D vector containing the (x, y, z) coordinates of shear layer crossing
        point.

    x_source : (3,) array_like
        1D vector containing the (x, y, z) coordinates of final point.

    x_mic : (3,) array_like
        1D vector containing the (x, y, z) coordinates of final point.

    Ux : float
        Mean flow velocity, assumed in '+x' direction.

    c0 : float
        Speed of sound.

    Returns
    -------
    t_total : float
        Time taken for sound to travel from source, through shear layer, to mic.

    Notes
    -----
    The corrected position is calculated with Amiet's shear layer correction
    method [JSV 58, 1978], and is defined as the point the acoustic ray would
    have reached if there was no shear layer (i.e. flow everywhere). The
    distance is corrected so that 'r(xc-xr) = r(xm-xs)'.
    """

    # calculate travel time inside flow (source to shear layer)
    tf = t_convect(XYZ_s, XYZ_sl, Ux, c0)

    # determine ray phase velocity in flow (direction and magnitude)
    cp_ray = (XYZ_sl-XYZ_s)/tf

    # travel time for corrected position
    tc = r(XYZ_m-XYZ_s)/c0

    # corrected position
    XYZ_c = XYZ_s + cp_ray*tc

    # retarded source position
    XYZ_r = XYZ_s + np.array([Ux, 0, 0])*tc

    return XYZ_c, XYZ_r


# %% *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# --->>> Flat plate aeroacoustic response functions


# *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# New versions, using numpy only - faster, no memory leak

def fr_integrand_re(x):
    """ Creates the argument to the Fresnel integral."""
    return (np.exp(1j*x)/np.sqrt(x)).real


def fr_integrand_im(x):
    """ Creates the argument to the complex conjugate Fresnel integral."""
    return (np.exp(1j*x)/np.sqrt(x)).imag


def fr_int(zeta):
    """
    Calculates the Fresnel integral of 'zeta'

    Parameters
    ----------
    zeta : (Nz,) array_like
        1D array of parameter 'zeta' for integration.

    Returns
    -------
    E : (Nz,) array_like
        1D array with results of Fresnel integral of each value of 'zeta'

    Notes
    -----
    Its complex-conjugate version can be obtained from the
    'amiet_tools.fr_int_cc' function.
    """

    # Check if zeta is array or float
    if type(zeta) is np.ndarray:
        E = np.zeros(zeta.shape, 'complex')

        # Calculate Fresnel integral for all non-zero values of zeta
        for i in range(zeta.size):
            if zeta[i] != 0:
                E[i] = (integrate.quad(fr_integrand_re, 0, zeta[i])[0]
                        + 1j*integrate.quad(fr_integrand_im, 0, zeta[i])[0])

    else:
        if zeta!=0:
            E = (integrate.quad(fr_integrand_re, 0, zeta)[0]
                 + 1j*integrate.quad(fr_integrand_im, 0, zeta)[0])

    return (1/np.sqrt(2*np.pi))*E


def fr_int_cc(zeta):
    """
    Calculates the complex-conjugate Fresnel integral of 'zeta'

    Parameters
    ----------
    zeta : (Nz,) array_like
        1D array of parameter 'zeta' for integration.

    Returns
    -------
    E_conj : (Nz,) array_like
        1D array with results of complex-conjugate Fresnel integral of each
        value of 'zeta'

    Notes
    -----
    Its non-complex-conjugate version can be obtained from the
    'amiet_tools.fr_int' function.
    """

   # Check if zeta is array or float
    if type(zeta) is np.ndarray:
        E_conj = np.zeros(zeta.shape, 'complex')

        # Calculate complex-conjugate Fresnel integral for all non-zero values
        # of zeta
        for i in range(zeta.size):
            if zeta[i] != 0:
                E_conj[i] = (integrate.quad(fr_integrand_re, 0, zeta[i])[0]
                             - 1j*integrate.quad(fr_integrand_im, 0, zeta[i])[0])

    else:
        if zeta!=0:
            E_conj = (integrate.quad(fr_integrand_re, 0, zeta)[0]
                      - 1j*integrate.quad(fr_integrand_im, 0, zeta)[0])

    return (1/np.sqrt(2*np.pi))*E_conj


def delta_p(rho0, b, w0, Ux, Kx, ky, xy, Mach):
    """
    Calculates the pressure jump response 'delta_p' for a single turbulent gust.

    Parameters
    ----------
    rho0 : float
        Density of air.

    b : float
        Airfoil semichord.

    w0 : float
        Gust amplitude.

    Ux : float
        Mean flow velocity, assumed in '+x' direction.

    Kx : float
        Chordwise turbulent gust wavenumber.

    ky : float
        Spanwise turbulent gust wavenumber.

    xy : ({2, 3}, Ny, Nx) array_like
        2D array containing (x, y) coordinates of airfoil surface mesh.

    Mach : float
        Mean flow Mach number.

    Returns
    -------
    delta_p : (Ny, Nx) array_like
        Surface pressure jump over airfoil surface mesh in response to a single
        turbulent gust with wavenumbers (Kx, ky) and amplitude 'w0'.
    """

    # pressure difference over the whole airfoil surface
    delta_p = np.zeros(xy[0].shape, 'complex')

    if xy.ndim == 3:
        # unsteady lift over the chord line (mid-span)
        g_x = np.zeros(xy[0][0].shape, 'complex')

        # calculates the unsteady lift over the chord
        g_x = g_LE(xy[0][0], Kx, ky, Mach, b)

        # broadcasts a copy of 'g_x' to 'delta_p'
        delta_p = g_x[np.newaxis, :]

    elif xy.ndim == 2:
        # unsteady lift over the chord line (mid-span)
        g_x = np.zeros(xy[0].shape, 'complex')

        # calculates the unsteady lift over the chord
        delta_p = g_LE(xy[0], Kx, ky, Mach, b)

    # adds the constants and the 'k_y' oscillating component
    delta_p = 2*np.pi*rho0*w0*delta_p*np.exp(-1j*ky*xy[1])

    return delta_p


def g_LE(xs, Kx, ky, Mach, b):
    """
    Airfoil non-dimensional chordwise pressure jump in response to a single gust.

    Parameters
    ----------
    xs : (Ny, Nx) or (Nx,) array_like
        Airfoil surface mesh chordwise coordinates.

    Kx : float
        Chordwise turbulent gust wavenumber.

    ky : float
        Spanwise turbulent gust wavenumber.

    Mach : float
        Mean flow Mach number.

    b : float
        Airfoil semichord.

    Returns
    -------
    g_LE : (Ny, Nx) array_like
        Non-dimensional chordwise surface pressure jump over airfoil surface
        mesh in response to a single turbulent gust with wavenumbers (Kx, ky)
        and amplitude 'w0'.

    Notes
    -----
    This function provides the airfoil responses for either subcritical or
    supercritical gusts. For critical gusts, the airfoil response is
    interpolated from slightly sub- and slightly supercritical responses.
    """

    beta = np.sqrt(1-Mach**2)
    ky_critical = Kx*Mach/beta

    # p_diff < 0: supercritical
    # p_diff > 0: subcritical
    p_diff = (np.abs(ky) - ky_critical)/ky_critical

    # supercritical gusts
    if p_diff < -1e-3:
        g = g_LE_super(xs, Kx, ky, Mach, b)

    # subcritical gusts
    elif p_diff > 1e-3:
        g = g_LE_sub(xs, Kx, ky, Mach, b)

    # critical gusts (interpolate between super- and subcritical)
    else:

        # get gusts 1% above and below critical ky
        ky_sp = ky*0.99
        ky_sb = ky*1.01

        g_sp = g_LE_super(xs, Kx, ky_sp, Mach, b)
        g_sb = g_LE_sub(xs, Kx, ky_sb, Mach, b)

        g = (g_sp + g_sb)/2.

    return g


def g_LE_super(xs, Kx, ky, Mach, b):
    """
    Returns airfoil non-dimensional pressure jump for supercritical gusts.

    Parameters
    ----------
    xs : (Ny, Nx) or (Nx,) array_like
        Airfoil surface mesh chordwise coordinates.

    Kx : float
        Chordwise turbulent gust wavenumber.

    ky : float
        Spanwise turbulent gust wavenumber.

    Mach : float
        Mean flow Mach number.

    b : float
        Airfoil semichord.

    Returns
    -------
    g_LE_super : (Ny, Nx) array_like
        Non-dimensional chordwise surface pressure jump over airfoil surface
        mesh in response to a single supercritical turbulent gust with
        wavenumbers (Kx, ky)

    Notes
    -----
    This function includes two terms of the Schwarzchild technique; the first
    term contains the solution for a infinite-chord airfoil with a leading edge
    but no trailing edge, while the second term contains a correction factor
    for a infinite-chord airfoil with a trailing edge but no leading edge.
    """

    beta = np.sqrt(1-Mach**2)
    mu_h = Kx*b/(beta**2)
    mu_a = mu_h*Mach

    kappa = np.sqrt(mu_a**2 - (ky*b/beta)**2)

    g1_sp = (np.exp(-1j*((kappa - mu_a*Mach)*((xs/b) + 1) + np.pi/4))
             / (np.pi*np.sqrt(np.pi*((xs/b) + 1)*(Kx*b + (beta**2)*kappa))))

    g2_sp = -(np.exp(-1j*((kappa - mu_a*Mach)*((xs/b) + 1) + np.pi/4))
              * (1-(1+1j)*fr_int_cc(2*kappa*(1-xs/b)))
              / (np.pi*np.sqrt(2*np.pi*(Kx*b + (beta**2)*kappa))))

    return g1_sp + g2_sp


def g_LE_sub(xs, Kx, ky, Mach, b):
    """
    Returns airfoil non-dimensional pressure jump for subcritical gusts.

    Parameters
    ----------
    xs : (Ny, Nx) or (Nx,) array_like
        Airfoil surface mesh chordwise coordinates.

    Kx : float
        Chordwise turbulent gust wavenumber.

    ky : float
        Spanwise turbulent gust wavenumber.

    Mach : float
        Mean flow Mach number.

    b : float
        Airfoil semichord.

    Returns
    -------
    g_LE_sub : (Ny, Nx) array_like
        Non-dimensional chordwise surface pressure jump over airfoil surface
        mesh in response to a single subcritical turbulent gust with
        wavenumbers (Kx, ky)

    Notes
    -----
    This function includes two terms of the Schwarzchild technique; the first
    term contains the solution for a infinite-chord airfoil with a leading edge
    but no trailing edge, while the second term contains a correction factor
    for a infinite-chord airfoil with a trailing edge but no leading edge.
    """

    beta = np.sqrt(1-Mach**2)
    mu_h = Kx*b/(beta**2)
    mu_a = mu_h*Mach

    kappa1 = np.sqrt(((ky*b/beta)**2) - mu_a**2)

    g1_sb = (np.exp((-kappa1 + 1j*mu_a*Mach)*((xs/b) + 1))*np.exp(-1j*np.pi/4)
             / (np.pi*np.sqrt(np.pi*((xs/b) + 1)
                              * (Kx*b - 1j*(beta**2)*kappa1))))

    g2_sb = -(np.exp((-kappa1 + 1j*mu_a*Mach)*((xs/b) + 1))
              * np.exp(-1j*np.pi/4)*(1 - ss.erf(2*kappa1*(1-xs/b)))
              / (np.pi*np.sqrt(2*np.pi*(Kx*b - 1j*(beta**2)*kappa1))))

    return g1_sb + g2_sb


def L_LE(x, sigma, Kx, ky, Mach, b):
    """
    Returns the effective lift functions - i.e. chordwise integrated surface pressures

    Parameters
    ----------
    x : (M,) array_like
        1D array of observer locations 'x'-coordinates

    sigma : (M,) array_like
        1D array of observer locations flow-corrected distances

    Kx : float
        Chordwise turbulent gust wavenumber.

    ky : float
        Spanwise turbulent gust wavenumber.

    Mach : float
        Mean flow Mach number.

    b : float
        Airfoil semichord.


    Returns
    -------
    L_LE : (M,) array_like
        Effective lift function for all observer locations.

    Notes
    -----
    These functions are the chordwise integrated surface pressures, and are
    parts of the far-field-approximated model for airfoil-turbulente noise.
    """

    beta = np.sqrt(1-Mach**2)
    ky_critical = Kx*Mach/beta

    # percentage difference in ky
    # p_diff < 0: supercritical / p_diff > 0: subcritical
    p_diff = (np.abs(ky) - ky_critical)/ky_critical

    # supercritical gusts
    if p_diff < -1e-3:
        L = L_LE_super(x, sigma, Kx, ky, Mach, b)

    # subcritical gusts
    elif p_diff > 1e-3:
        L = L_LE_sub(x, sigma, Kx, ky, Mach, b)

    # critical gusts (interpolate between super- and subcritical)
    else:
        # get gusts 1% above and below critical ky
        ky_sp = ky*0.99
        ky_sb = ky*1.01

        L_sp = L_LE_super(x, sigma, Kx, ky_sp, Mach, b)
        L_sb = L_LE_sub(x, sigma, Kx, ky_sb, Mach, b)

        L = (L_sp + L_sb)/2.

    return L



def L_LE_super(x, sigma, Kx, Ky, Mach, b):
    """
    Returns the effective lift functions for supercritical gusts

    Parameters
    ----------

    x : (M,) array_like
        1D array of observer locations 'x'-coordinates

    sigma : (M,) array_like
        1D array of observer locations flow-corrected distances

    Kx : float
        Chordwise turbulent gust wavenumber.

    ky : float
        Spanwise turbulent gust wavenumber.

    Mach : float
        Mean flow Mach number.

    b : float
        Airfoil semichord.


    Returns
    -------

    Notes
    -----
    These functions are the chordwise integrated surface pressures, and are
    parts of the far-field-approximated model for airfoil-turbulente noise.
    """

    beta = np.sqrt(1-Mach**2)
    mu_h = Kx*b/(beta**2)
    mu_a = mu_h*Mach

    kappa = np.sqrt(mu_a**2 - (Ky*b/beta)**2)
    H1 = kappa - mu_a*x/sigma
    H2 = mu_a*(Mach - x*sigma) - np.pi/4

    L1 = ((1/np.pi)*np.sqrt(2/((Kx*b + (beta**2)*kappa)*H1))
          * fr_int_cc(2*H1)*np.exp(1j*H2))

    L2 = ((np.exp(1j*H2)
          / (np.pi*H1*np.sqrt(2*np.pi*(Kx*b + (beta**2)*kappa))))
          * (1j*(1 - np.exp(-2j*H1))
             + (1 - 1j)*(fr_int_cc(4*kappa)
                         - np.sqrt(2*kappa/(kappa + mu_a*x/sigma))
                         * np.exp(-2j*H1)
                         * fr_int_cc(2*(kappa + mu_a*x/sigma)))))

    return L1+L2


def L_LE_sub(x, sigma, Kx, Ky, Mach, b):
    """
    Returns the effective lift functions for subcritical gusts

    Parameters
    ----------

    x : (M,) array_like
        1D array of observer locations 'x'-coordinates

    sigma : (M,) array_like
        1D array of observer locations flow-corrected distances

    Kx : float
        Chordwise turbulent gust wavenumber.

    ky : float
        Spanwise turbulent gust wavenumber.

    Mach : float
        Mean flow Mach number.

    b : float
        Airfoil semichord.

    Returns
    -------

    Notes
    -----
    These functions are the chordwise integrated surface pressures, and are
    parts of the far-field-approximated model for airfoil-turbulente noise.
    """

    beta = np.sqrt(1-Mach**2)
    mu_h = Kx*b/(beta**2)
    mu_a = mu_h*Mach

    kappa1 = np.sqrt((Ky*b/beta)**2 - (mu_a**2))
    H2 = mu_a*(Mach - x*sigma) - np.pi/4
    H3 = kappa1 - 1j*mu_a*x/sigma

    L1 = ((1/np.pi)*np.sqrt(2/((Kx*b - 1j*(beta**2)*kappa1)
                               * (1j*kappa1 - mu_a*x/sigma)))
          * fr_int(2*(1j*kappa1 - mu_a*x/sigma))*np.exp(1j*H2))

    L2 = ((1j*np.exp(1j*H2)
           / (np.pi*H3*np.sqrt(2*np.pi*(Kx*b - 1j*(beta**2)*kappa1))))
          * (1 - np.exp(-2*H3) - ss.erf(np.sqrt(4*kappa1))
              + 2*np.exp(-2*H3)*np.sqrt(kappa1/(1j*kappa1 + mu_a*x/sigma))
              * fr_int(2*(1j*kappa1 - mu_a*x/sigma))))

    return L1+L2


# %% *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# Hydrodynamic spanwise wavenumber sampling

def ky_vector(b, d, k0, Mach, beta, method='AcRad', xs_ref=None):
    """
    Returns a vector of spanwise gust wavenumbers for acoustic calculations

    Parameters
    ----------
    b : float
        Airfoil semi chord.

    d : float
        Airfoil semi span.

    k0 : float
        Acoustic wavenumber 'k0'. Can be obtained from the
        temporal frequency 'f' [in Hz] and the speed of sound 'c0' [in m/s]
        as 'k0 = 2*pi*f/c0'.'

    Mach : float
        Mean flow Mach number.

    beta : float
        Prandtl–Glauert parameter (=sqrt(1-M**2))

    method : {'AcRad', 'SurfPressure'}, optional
        Calculation method to use. Defaults to 'AcRad'.

    xs_ref : float, optional
        Chordwise coordinate of reference point, defined in interval (-b, +b].
        Used in 'SurfPressure' mode, not required for 'AcRad' mode. Defaults to
        None.

    Returns
    -------
    Ky : (Nk,) array_like
        1D array containing spanwise gust wavenumbers in range [-ky_max, +ky_max],
        with center sample at ky=0

    Notes
    -----
    Returns a vector of equally-spaced spanwise hydrodynamic (gust) wavenumber
    values for calculations of airfoil response, either for calculating
    the airfoil acoustic radiation (method = 'AcRad') or for calculating the
    airfoil surface pressure cross-spectrum (method = 'SurfPressure').

    'AcRad' mode returns a shorter range of gusts for acoustic radiation
    calculations, and 'SurfPressure' returns a larger range for unsteady
    surface pressure calculations (excessive for acoustic radiation
    calculations).
    """

    # Assert 'method' string for valid inputs
    method_error = "'method' not recognized; please use either 'AcRad' or 'SurfPressure'"
    assert method in ['AcRad', 'SurfPressure'], method_error

    # critical hydrodynamic spanwise wavenumber
    ky_crit = k0/beta

    # width of ky sinc function
    sinc_width = 2*np.pi/(2*d)

    # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # for acoustic radiation calculations:
    if method == 'AcRad':

        if ky_crit < 2*np.pi/d:
            # 'low freq' - include some subcritical gusts (up to 1st sidelobe
            # of sinc function)
            N_ky = 41           # value obtained empirically
            ky_T = 2*np.pi/d    # main lobe + 1st sidelobes in sinc function
            Ky = np.linspace(-ky_T, ky_T, N_ky)

        else:
            # 'high freq' - restrict to supercritical gusts only

            # get ky with spacing equal to approx. 1/8 width of sinc function
            N_ky = int(np.ceil(2*ky_crit/sinc_width)*8)+1
            Ky = np.linspace(-ky_crit, ky_crit, N_ky)

    # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # for surface pressure cross-spectra calculations
    elif method == 'SurfPressure':

        # find ky that is at -20 dB at reference chord point
        ky_20dBAtt = ky_att(xs_ref, b, Mach, k0, Att=-20)

        # largest ky under consideration (25% above ky_20dBAtt, for safety)
        ky_max = 1.25*ky_20dBAtt

        # get ky with spacing equal to approx. 1/8 width of sinc function
        N_ky = int(np.ceil(2*ky_max/sinc_width)*8)+1

        Ky = np.linspace(-ky_max, ky_max, N_ky)

    # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

    return Ky


def ky_att(xs, b, Mach, k0, Att=-20):
    """
    Returns the spanwise gust wavenumber 'ky_att' with response at 'xs' attenuated by 'Att' decibels

    Parameters
    ----------
    xs : float
        Chordwise coordinate of reference point, defined in interval (-b, +b].

    b : float
        Airfoil semi chord.

    Mach : float
        Mean flow Mach number.

    k0 : float
        Acoustic wavenumber 'k0'. Can be obtained from the
        temporal frequency 'f' [in Hz] and the speed of sound 'c0' [in m/s]
        as 'k0 = 2*pi*f/c0'.

    Att : float, optional
        Level of attenuation of the surface pressure at point 'xs', in decibels.
        Defaults to -20 dB.

    Returns
    -------
    ky_att : float
        Subcritical gust spanwise wavenumber 'ky_att' such that the aerofoil
        response at point 'xs' is 'Att' dB reduced.
    """

    beta = np.sqrt(1-Mach**2)

    # critical gust spanwise wavenumber
    ky_crit = k0/beta

    term1 = -(beta**2)*np.log(10**(Att/20))/(k0*(xs + b))

    return ky_crit*np.sqrt(term1**2 + 1)


# %% *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# turbulent velocity spectra

def Phi_2D(Kx, ky_vec, Ux, turb_intensity, length_scale, model='K'):
    """
    Returns 2D isotropic turbulence energy spectrum in wavenumber domain.

    Parameters
    ----------
    Kx : (N_kx,) array_like or float
        1D array of chordwise gust wavenumbers.

    ky_vec : (N_ky,) array_like
        1D array of spanwise gust wavenumbers.

    Ux : float
        Mean flow velocity, assumed in '+x' direction.

    turb_intensity : float
        Turbulence intensity: sqrt(w_meanSquared/(Ux**2))

    length_scale : float
        Turbulence integral length scale, in meters

    model : {'K', 'L'}
        Type of spectrum: 'K' for von Karman spectrum, or 'L' for Liepmann spectrum.

    Returns
    -------
    Phi : (N_kx, N_ky) or (N_ky,) array_like
        2D or 1D array containing the values of two-dimensional turbulence
        energy for each wavenumber, according to von Karman or Liepmann
        spectrum.
    """

    u_mean2 = (Ux*turb_intensity)**2

    if type(Kx) is not np.ndarray:
        Kx = np.asarray([Kx])

    # von Karman model (Amiet 1975)
    if model == 'K':
        ke = (np.sqrt(np.pi)/length_scale)*(ss.gamma(5./6)/ss.gamma(1./3))

        kxe2_ye2 = (Kx[:, np.newaxis]/ke)**2 + (ky_vec/ke)**2

        return (4./(9*np.pi))*(u_mean2/(ke**2))*kxe2_ye2/((1+kxe2_ye2)**(7./3))

    # 2D Liepmann turbulence spectrum
    elif model == 'L':

        ls2 = length_scale**2

        return ((u_mean2*ls2/(4*np.pi))
                * ((1+ls2*(4*Kx[:, np.newaxis]**2 + ky_vec**2)))
                / (1+ls2*(Kx[:, np.newaxis]**2 + ky_vec**2))**(5./2))


# %% *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# --->>> Other assorted functions


def rect_grid(grid_sides, point_spacings):
    """
    Returns the 2D coordinates for a uniformly spaced rectangular grid.

    Parameters
    ----------
    grid_sides : (2,) array_like
        1D array containing the lengths 'Lx' and 'Ly' of the grid in the 'x'
        and 'y' directions, respectively.

    point_spacings : (2,) array_like
        1D array containing the spacings 'dx' and 'dy' between the points in
        the 'x' and 'y' directions, respectively.

    Returns
    -------
    XY_grid : (2, Nx*Ny) array_like
        2D array containing the (x, y) coordinates of all points in the grid.
        The grid contains 'Nx = Lx/dx+1' points in the 'x' direction, and
        'Ny = Ly/dy+1' points in the 'y' direction.
    """

    # number of points on each side = Dx/dx + 1
    N_points = np.array([round(grid_sides[0]/point_spacings[0] + 1),
                         round(grid_sides[1]/point_spacings[1] + 1)],
                        dtype='int')

    x_points = np.linspace(-grid_sides[0]/2., grid_sides[0]/2., N_points[0])
    y_points = np.linspace(-grid_sides[1]/2., grid_sides[1]/2., N_points[1])

    X_points, Y_points = np.meshgrid(x_points, y_points)

    return np.array([X_points.flatten(), Y_points.flatten()])


def read_ffarray_lvm(filename, n_columns = 13):
    """
    Reads a .lvm file containing time-domain data acquired from LabView.

    Parameters
    ----------
    filename : string
        Name of the '.lvm' file to be read.

    N_columns : int, optional
        Number of columns to read (time samples + (N-1) signals). Defaults to 13.

    Returns
    -------
    t : (N,) array_like
        1D array containing the time samples.

    mics : (M, N) array_like
        2D array containing the 'M' microphone signals.

    Notes
    -----
    The .lvm file is assuemd to contain the time samples in the first column
    and each microphone signal in the remaining columns. Default is time signal
    plus 12 microphones.
    """

    # open file as text
    lvm_file = open(filename, 'r')

    # count the number of lines (has to read the whole file,
    # hence closing and opening again)
    n_lines = len(list(lvm_file))
    lvm_file.close()

    t = np.zeros(n_lines)
    mics = np.zeros((n_lines, n_columns-1))

    # read line ignoring the '\r\n' at the end,
    # and using '\t' as column separators
    lvm_file = open(filename, 'r')
    for line in range(n_lines):
        current_line = lvm_file.readline().split('\r\n')[0].split('\t')
        t[line] = float(current_line[0])
        mics[line, :] = [float(x) for x in current_line][1:]

    # close file
    lvm_file.close()

    return t, mics.T


# %% *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# Old aeroacoustics versions, using mpmath - slower, leaks memory when running
# for a long time

# import mpmath as mp

# def fr_integrand(x):
#     """ Creates the argument to the Fresnel integral."""
#     return mp.exp(1j*x)/mp.sqrt(x)


# def fr_integrand_cc(x):
#     """ Creates the argument to the complex conjugate Fresnel integral."""
#     return mp.exp(-1j*x)/mp.sqrt(x)


# def fr_int(zeta):
#     """
#     Calculates the Fresnel integral of 'zeta'

#     Parameters
#     ----------
#     zeta : (Nz,) array_like
#         1D array of parameter 'zeta' for integration.

#     Returns
#     -------
#     E_conj : (Nz,) array_like
#         1D array with results of Fresnel integral of each value of 'zeta'

#     Notes
#     -----
#     This function uses the module 'mpmath' for high precision quadrature of the
#     integrand.

#     Its complex-conjugate version can be obtained from the
#     'amiet_tools.fr_int_cc' function.
#     """

#     # Check if zeta is array or float
#     if type(zeta) is np.ndarray:
#         E_conj = np.zeros(zeta.shape, 'complex')

#         for i in range(zeta.size):
#             E_conj[i] = ((1/np.sqrt(2*np.pi))
#                          * (mp.quad(fr_integrand, [0, zeta[i]])))
#     else:
#         E_conj = ((1/np.sqrt(2*np.pi))
#                   * (mp.quad(fr_integrand, [0, zeta])))

#     return E_conj


# def fr_int_cc(zeta):
#     """
#     Calculates the complex-conjugate Fresnel integral of 'zeta'

#     Parameters
#     ----------
#     zeta : (Nz,) array_like
#         1D array of parameter 'zeta' for integration.

#     Returns
#     -------
#     E_conj : (Nz,) array_like
#         1D array with results of complex-conjugate Fresnel integral of each
#         value of 'zeta'

#     Notes
#     -----
#     This function uses the module 'mpmath' for high precision quadrature of the
#     integrand.

#     Its non-complex-conjugate version can be obtained from the
#     'amiet_tools.fr_int' function.
#     """

#     # Check if zeta is array or float
#     if type(zeta) is np.ndarray:
#         E_conj = np.zeros(zeta.shape, 'complex')

#         for i in range(zeta.size):
#             E_conj[i] = ((1/np.sqrt(2*np.pi))
#                          * (mp.quad(fr_integrand_cc, [0, zeta[i]])))
#     else:
#         E_conj = ((1/np.sqrt(2*np.pi))
#                   * (mp.quad(fr_integrand_cc, [0, zeta])))

#     return E_conj


# def g_LE(xs, Kx, ky, Mach, b):
#     """
#     Airfoil non-dimensional chordwise pressure jump in response to a single gust.

#     Parameters
#     ----------
#     xs : (Ny, Nx) or (Nx,) array_like
#         Airfoil surface mesh chordwise coordinates.

#     Kx : float
#         Chordwise turbulent gust wavenumber.

#     ky : float
#         Spanwise turbulent gust wavenumber.

#     Mach : float
#         Mean flow Mach number.

#     b : float
#         Airfoil semichord.

#     Returns
#     -------
#     g_LE : (Ny, Nx) array_like
#         Non-dimensional chordwise surface pressure jump over airfoil surface
#         mesh in response to a single turbulent gust with wavenumbers (Kx, ky)
#         and amplitude 'w0'.

#     Notes
#     -----
#     This function provides the airfoil responses for either subcritical or
#     supercritical gusts. For critical gusts, the airfoil response is
#     interpolated from slightly sub- and slightly supercritical responses.
#     """

#     beta = np.sqrt(1-Mach**2)
#     ky_critical = Kx*Mach/beta

#     # p_diff < 0: supercritical
#     # p_diff > 0: subcritical
#     p_diff = (np.abs(ky) - ky_critical)/ky_critical

#     # supercritical gusts
#     if p_diff < -1e-3:
#         g = g_LE_super(xs, Kx, ky, Mach, b)

#     # subcritical gusts
#     elif p_diff > 1e-3:
#         g = g_LE_sub(xs, Kx, ky, Mach, b)

#     # critical gusts (interpolate between super- and subcritical)
#     else:

#         # get gusts 1% above and below critical ky
#         ky_sp = ky*0.99
#         ky_sb = ky*1.01

#         g_sp = g_LE_super(xs, Kx, ky_sp, Mach, b)
#         g_sb = g_LE_sub(xs, Kx, ky_sb, Mach, b)

#         g = (g_sp + g_sb)/2.

#     # if single mp_complex, convert to float
#     if type(g) is mp.ctx_mp_python.mpc:
#         # convert to single float
#         return float(g.real) + 1j*float(g.imag)

#     # if single float, return as is
#     elif type(g) is np.complex128:
#         return g

#     # if array of mp_complex, convert to np.ndarray
#     elif type(g) is np.ndarray and g.dtype is np.dtype(mp.ctx_mp_python.mpc):
#         return np.array([float(x.real) + 1j*float(x.imag) for x in g])

#     # if array of complex floats, just return as is
#     elif type(g) is np.ndarray and g.dtype is np.dtype(np.complex128):
#         return g


# def L_LE(x, sigma, Kx, ky, Mach, b):
#     """
#     Returns the effective lift functions - i.e. chordwise integrated surface pressures

#     Parameters
#     ----------
#     x : (M,) array_like
#         1D array of observer locations 'x'-coordinates

#     sigma : (M,) array_like
#         1D array of observer locations flow-corrected distances

#     Kx : float
#         Chordwise turbulent gust wavenumber.

#     ky : float
#         Spanwise turbulent gust wavenumber.

#     Mach : float
#         Mean flow Mach number.

#     b : float
#         Airfoil semichord.


#     Returns
#     -------
#     L_LE : (M,) array_like
#         Effective lift function for all observer locations.

#     Notes
#     -----
#     These functions are the chordwise integrated surface pressures, and are
#     parts of the far-field-approximated model for airfoil-turbulente noise.
#     """

#     beta = np.sqrt(1-Mach**2)
#     ky_critical = Kx*Mach/beta

#     # percentage difference in ky
#     # p_diff < 0: supercritical / p_diff > 0: subcritical
#     p_diff = (np.abs(ky) - ky_critical)/ky_critical

#     # supercritical gusts
#     if p_diff < -1e-3:
#         L = L_LE_super(x, sigma, Kx, ky, Mach, b)

#     # subcritical gusts
#     elif p_diff > 1e-3:
#         L = L_LE_sub(x, sigma, Kx, ky, Mach, b)

#     # critical gusts (interpolate between super- and subcritical)
#     else:
#         # get gusts 1% above and below critical ky
#         ky_sp = ky*0.99
#         ky_sb = ky*1.01

#         L_sp = L_LE_super(x, sigma, Kx, ky_sp, Mach, b)
#         L_sb = L_LE_sub(x, sigma, Kx, ky_sb, Mach, b)

#         L = (L_sp + L_sb)/2.

#     # if single mp_complex, convert to float
#     if type(L) is mp.ctx_mp_python.mpc:
#         # convert to single float
#         return float(L.real) + 1j*float(L.imag)

#     # if single float, return as is
#     elif type(L) is np.complex128:
#         return L

#     # if array of mp_complex, convert to np.ndarray
#     elif type(L) is np.ndarray and L.dtype is np.dtype(mp.ctx_mp_python.mpc):
#         return np.array([float(x.real) + 1j*float(x.imag) for x in L])

#     # if array of complex floats, just return as is
#     elif type(L) is np.ndarray and L.dtype is np.dtype(np.complex128):
#         return L
