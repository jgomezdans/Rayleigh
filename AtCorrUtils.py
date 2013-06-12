#!/usr/bin/env python
"""
This file contains a number of support functions to perform
Atmospheric correction over inland waters. They will be used
by a number of different atmospheric correction algorithms.
"""
import sys
import os
import numpy as np
from osgeo import gdal
import matplotlib.pyplot as plt

def get_o3conc ( doy, year, fname="o3_conc.txt" ):
    """Get O3 concentration for a given year and DoY for the region of interest. 
    
    If you need other areas, you will need to provide a new file for this!!!! We 
    use a very simple linear interpolator of the time series given in the file. File
    format is two columns: date and concentration in dobson units (???TODO check!).
    Time format is YYYYDDD.
    
    Parameters
    ----------
    doy: int
        The Day of Year
    year: int
        The year (tested between 1980 and 2012, but might change that if needed)
    fname: str, optional
        The filename with the timeseries data.
    """
    if  ( year < 1980 ) or ( year > 2012 ):
        print "Only years between 1980 are considered with current data file"
        raise ValueError
    d = np.loadtxt( fname )
    missing_dates = d[ d[:,1] == 0, 0 ] # Missing dates in format YYYYDDD
    present_dates = d[ d[:,1] != 0, 0 ]
    present_meas = d[ d[:, 1] != 0, 1 ] # the actual measurements
    target_period = float ( "%4d%03d" % ( year, doy ) )
    o3_conc_for_today = np.interp ( target_period, present_dates, present_meas )
    return o3_conc_for_today

def get_angles ( fname ):
    """Extract acquisition geometry from MLT file
    
    Parameters
    ----------
    fname: str
        A filename
    """
    if fname.find ( "MTL.txt" ) < 0:
        print "Filename doesn't have 'MTL.txt'. Are you sure this is the right file?"
        raise IOerror
    fp = open( fname, 'r') 
    
    for line in fp: 
    
        if ( line.find ("SUN_ELEVATION") >= 0 ):
            s = line.split("=") 
            theta_0 = float(s[1]) 
        
        elif ( line.find ("SUN_AZIMUTH") >= 0 ):
            s = line.split("=") 
            phi_0 = float(s[1])
    solar_angles=(theta_0,phi_0)
    return solar_angles

def get_lambdai(fname):
    """Returns the mean wavelength value in nanometers.
    
    Parameters
    ----------
    fname: str
        The metadata filename. Assumes USGS convention
    """

    if "LE" in fname:
        return np.array ( [ 482.5, 565, 660., 837.5, 1650, 11450, 2220 ] )
    if "LT" in fname:
        return np.array ( [ 485.2, 560, 660, 830., 1650, 11450., 2215 ] )
            

def get_doy(fname):
    """Finds the day of the year (doy) from Landsat filename
    
    Parameters
    ----------
    fname: str
        The metadata filename. Assumes USGS convention

    """
    fname = os.path.basename ( fname )
    try:
        print fname
        doy = int ( fname[13:16] )
        year = int ( fname[9:13] )
    except:
        print "Problem with elucidating DoY from filename!"
    return doy, year

def get_date(doy,year):
    """ Getting date from doy format"""
    import datetime
    return datetime.datetime(year, 1, 1) + datetime.timedelta(doy - 1)


def solar_distance ( doy ):
    """Calculates the solar distance as a function of day of year (DoY).
    
    Code taken from 6S, subroutine `VARSOL`
    
    Parameters
    ----------
    doy: int 
        Day of year
    """
    #assert ( doy > 0 ) & ( doy < 367 )
    om = (0.9856*( doy -4))*np.pi/180.
    dsol = 1. / ( (1.-.01673*np.cos(om))**2)
    return dsol

def extraterrestrial_radiation( wvl, doy ):
    """Need a function to calculate **spectral** extraterrestrial radiation as a function of Earth-Sun
    distance and wavelength. I think a formula is available in Liang's book, but haven't got it to
    hand...."""
    # TODO
    import urllib
    import os
    # Check if the coeffs are local. If they aren't, download them
    if not os.path.exists ( "solirr.dat" ):
        url="https://gist.github.com/jgomezdans/5460835/raw/f55b43be870fb8809e480f9198c943fbb53047d4/solirr.dat"
        data = np.loadtxt( urllib.urlopen ( url ) )
        np.savetxt ( "solirr.dat", data )
        solirr = data
    else:
        solirr = np.loadtxt ( "solirr.dat" )
    sun_earth_distance = solar_distance ( doy )
    i = np.argmin( np.abs(solirr[:,0] - wvl))
    ETo = solirr[i,1] * sun_earth_distance
    return ETo

def ozone ( lambdai, theta_0, theta, o3_column ):
    """Calculation of the ozone transmittance $T_{O3}(\lambda)$
    We need to get the coefficients which are available in 
    https://gist.github.com/jgomezdans/5443793/raw/fc5e78fb0317d8900dd2db1613511d7060c1c2bd/o3_coeff.txt

    Parameters
    ----------
    lambdai: wavelength in nm
    theta_0: Sun zenith angle (SZA) in degrees
    theta: View zenith angle (VZA) in degrees
    o3_column: Total columnar O3 concentration in Dobson units
    """
    import urllib
    import os
    # Check if the coeffs are local. If they aren't, download them
    if not os.path.exists ( "o3_coeff.txt" ):
        url="http://gist.github.com/jgomezdans/5443793/raw/fc5e78fb0317d8900dd2db1613511d7060c1c2bd/o3_coeff.txt"
        data = np.loadtxt( urllib.urlopen ( url ) )
        np.savetxt ( "o3_coeff.txt", data )
        o3_coeff = data
    else:
        o3_coeff = np.loadtxt ( "o3_coeff.txt" )
    o3_column = o3_column/1000. # Convert to atm cm
    # Calculate mu terms...
    mu = 1./np.cos ( np.deg2rad(theta) )
    mu_0 = 1./np.cos ( np.deg2rad( theta_0 ) )
    iloc = int( np.ceil ( lambdai - o3_coeff[:,0].min()) )
    if lambdai >= 974:
        tau_O3 = 0.0
    else:
        tau_O3 = o3_coeff[ iloc, 1] * o3_column
    T_O3 = np.exp ( -tau_O3*(1./mu + 1./mu_0 ))
    return T_O3

def rayleigh_optical_depth ( lambdai, h0 ):
    """
    Calculates the rayleigh optical depth. Uses an approximation from Hansen and Travis, but 
    it is not clear where the sensor height function comes from?

    Note that we might also consider tau_r as being given as in Bodhaine et al (1999) (Eq. 30)
    tau_r(\lambda in microns) = 0.0021520*( (1.0455996-341.29061*(1./lambda**2) - 0.90230850*lambda**2)/ 
          (1 + 0.0027059889*(1./lambda**2) - 85.968563*lambda**2))

    This is for sea level, latitude of 45deg and 1 atm. 
    Parameters
    ----------
    lambdai: wavelength ***in microns***
    h0: sensor height in km
    """
    lambdai = lambdai/1000.
    Hr = np.exp ( -0.1188*h0 - 0.0011*h0*h0 )
    tau_r = Hr*( 0.00859*(1/(lambdai**4))*(1+0.0113*(1/(lambdai**2)) + 0.00013*(1/(lambdai**4))))
    return tau_r

def rayleigh_phase_functions ( theta_0, theta, phi_0, phi ):
    """Water phase functions require the acquisition geometry in degrees"""
    theta_0 = np.deg2rad ( theta_0 )
    theta = np.deg2rad ( theta )
    phi_0 = np.deg2rad ( phi_0 )
    phi = np.deg2rad ( phi )
    mu = np.cos ( theta )
    mu_0 = np.cos ( theta_0 )
    cos_gamma_down = np.cos ( theta_0 )*np.cos(theta) - np.sin(theta_0)*np.sin(theta)*np.cos(phi-phi_0)
    cos_gamma_up = - np.cos ( theta_0 )*np.cos(theta) - np.sin(theta_0)*np.sin(theta)*np.cos(phi-phi_0)
    
    P_gamma_down = (3./4.)*(1 + cos_gamma_down*cos_gamma_down )
    P_gamma_up = (3./4.) * (1 + (2*mu*mu_0 + cos_gamma_up)**2 )

    return P_gamma_down, P_gamma_up

def fresnel_reflectance ( theta, theta_0, m=1.396 ):
    mu = np.cos ( np.deg2rad ( theta ) )
    mu0 = np.cos ( np.deg2rad ( theta_0 ) )
    
    y = (1./m)*np.sqrt ( m*m + mu*mu - 1 )
    y0 = (1./m)*np.sqrt ( m*m + mu0*mu0 - 1 )
    rho_mu = 1. - 2.*mu*y*m*( 1./(mu + m*y)**2 + 1./(m*mu + y)**2)
    rho_mu0 = 1. - 2.*mu0*y0*m*( 1./(mu0 + m*y0)**2 + 1./(m*mu0 + y0)**2)
    return rho_mu, rho_mu0

def rayleigh_scattering ( theta, theta_0, phi, phi_0, lambdai, o3_conc, h0, doy ):
    # Calculate the terms in Eq. 2 of Wang et al
    # First, the relevant ETr... Note this is not implemented
    F0 = extraterrestrial_radiation ( lambdai, doy)
    
    # Ozone transmittance
    T_O3 = ozone ( lambdai, theta_0, theta, o3_conc )
    # Rayleigh optical depth
    tau_rayleigh = rayleigh_optical_depth ( lambdai, h0 )
    
    # The Rayleigh up and down phase functions
    P_gamma_down, P_gamma_up = rayleigh_phase_functions ( theta_0, theta, phi_0, phi )
    # Fresnel terms
    rho_mu, rho_mu0 = fresnel_reflectance ( theta, theta_0 )
    # Now put together Eq. 2...
    mu = np.cos ( np.deg2rad ( theta ) )
    rayleigh_radiance = (1./(4*np.pi*mu))
    rayleigh_radiance = rayleigh_radiance*(F0*T_O3*tau_rayleigh)
    rayleigh_radiance = rayleigh_radiance*(P_gamma_down + (rho_mu + rho_mu0)*P_gamma_up)
    diffuse = np.exp ( - ( (rayleigh_radiance/2.) + T_O3 )*mu )
            
    return rayleigh_radiance, diffuse

def water_leaving_rad_b7 ( diff_coeff, fname, rayleigh_radiance ):
    """Implements equations 9 and 10 of Wang and su puta madre.
    
    """
    fname = fname.replace ( "MTL.txt", "B7_TOARAD.tif" )
    g = gdal.Open ( fname )
    total_rad = g.ReadAsArray()
    water_leaving_rad = ( total_rad - rayleigh_radiance )/diff_coeff
    return water_leaving_rad

def nearest_neighbour_interpolation ( clear_water, green_aot ):
    """Nearest neighbour interpolation of green AOT for turbid water pixels.
    
    `clear_water` is a 2D array with values 2 for turbid water pixels and
    1 for clear water pixels. The aim is to interpolate the values of
    `green_aot` where `clear_water` is 1 for the locations where `clear_water` 
    is 2. This is a quick'n'ready nearest neighbour interpolator, which suffers
    from not being very efficient when many pixels have to be interpolated.
    
    Parameters
    ----------
    clear_water: array
        A 2D array with 1 indicating clear water pixels and 2 indicating turbid water pixels
    green_aot: array
        A 2D array that stores the AOT in the green for clear water pixels
        
    Returns:
    green_aot_intp: array
        A 2D array with the aerosol scattering contribution for both turbid and non-turbid pixels
    """

    yc, xc = np.nonzero ( clear_water == 2 )
    y, x = np.nonzero ( clear_water == 1 )
    
    v = [ np.hypot ( x - xc[i], y - yc[i]).argmin() for i in xrange(len(xc))]
    v = np.array ( v )
    green_aot_intp = green_aot*1.
    green_aot_intp [yc, xc] = green_aot [y[v], x[v]]
    return green_aot_intp

def aerosol_correction ( tau_diff, fname, l_rayleigh, doy, lambdas, theta_i ):
    """Aerosol correction from Wang et al. (2007).
    
    The following function implements the aerosol correctin described in Wang
    et al. (2007).The procedure is as follows:
    
    1. Use band 7 to find "clear water pixels" using a threshold
    2. Using the green band, and assuming that water leaving radiance is 
       fairly constant for clear water pixels, estimate the AOD at the green
       for these pixels.
    3. Interpolate green AOD for other non-clear water pixels
    4. Extrapolate the green AOD to the other visible bands, and solve for 
       water leaving radiance
    5. Additionally, convert to reflectance?
    
    Parameters
    ----------
    tau_diff: array
        The diffuse transmittance of the atmosphere for all TM bands
    fname: str
        A metadata (_MTL.txt) filename. 
    l_rayleigh: array
        Rayleigh radiance for all TM bands
    doy: int
        Day of year
    lambdas: array
        Centre frequencies of all TM bands
    theta_i: float
        The solar zenith angle (in degrees?)
        
    Returns
    --------
    water_leaving_rad_vis: array
        Water leaving radiance for bands 1, 2, 3 (B, G, R)
     
    NOTE: filenames follow a convention!!!!
    """
    water_leaving_rad = water_leaving_rad_b7 ( tau_diff[-1], fname, l_rayleigh[-1] )
    clear_water = np.where ( np.logical_and ( water_leaving_rad >= -0.05, \
        water_leaving_rad <= 0.2), 1, 0 ) 
    clear_water =  np.where ( np.logical_and ( water_leaving_rad > 0.2, \
        water_leaving_rad <= 0.9), 2, clear_water )
    clear_water =  np.where ( np.abs( water_leaving_rad )> 0.9, 0, clear_water )
    
    et_rad = extraterrestrial_radiation( doy, lambdas[1] )
    distance = solar_distance ( doy )
    green_refl = 0.054 

    for i in xrange ( 3 ):
        fname = fname.replace ( "MTL.txt", "B%d_TOARAD.tif" % ( i+1 ) )
        g = gdal.Open ( fname )
        if i == 0:
            total_rad = np.zeros ( ( 3, g.RasterXSize, g.RasterYSize ) )
        if i == 1:
            green_rad_toa = g.ReadAsArray()
        total_rad[i, :, : ] = g.ReadAsArray()
            
    
    green_rad = green_refl*np.cos(theta_i)*et_rad/(np.pi*distance)
    aerosol_corr = green_rad_toa*0.
    aerosol_corr[clear_water==1] = green_rad_toa[clear_water==1] - \
        tau_diff[1]*green_rad - l_rayleigh[1]
    aerosol_corr = nearest_neighbour_interpolation ( clear_water, aerosol_corr )
    s_factor = np.zeros(3)
    s_factor = np.array ( [ extraterrestrial_radiation( doy, lambdas[i] )/et_rad \
        for i in xrange(3) ] )
    water_leaving_rad_vis = np.zeros( ( 3, green_rad_toa.shape[0], green_rad_toa.shape[1] ) )
    for i in xrange(3):
        water_leaving_rad_vis[ i, :, : ] = ( total_rad[i,:,:] - l_rayleigh[i] - \
            aerosol_corr*s_factor[i] )/tau_diff[i]
    
    
    return water_leaving_rad_vis    