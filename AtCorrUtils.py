#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

def solar_distance ( doy ):
  """Calculates the solar distance as a function of
  day of year (DoY). Code taken from 6S, subroutine
  VARSOL"""
  assert ( doy > 0 ) & ( doy < 367 )
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
    tau_O3 = data[data[:,0]==lambdai, 1] * o3_column
    T_O3 = np.exp ( tau_O3*(1./mu + 1./np.cos(mu_0) ) )
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
    
    Hr = np.exp ( -0.1188*h0 - 0.0011*h0*h0 )
    tau_r = Hr*( 0.00859*(1/(lambdai**4))*(1+0.0113*(1/(lambdai**2)) + 0.00013*(1/(lambdai**4))))
    return tau_r

def rayleigh_phase_functions ( theta_0, theta, phi_0, phi ):
    """Wather phase functions require the acquisition geometry in degrees"""
    theta_0 = np.deg2rad ( theta_0 )
    theta = np.deg2rad ( theta )
    phi_0 = np.deg2rad ( phi_0 )
    phi = np.deg2rad ( phi )
    mu = np.cos ( theta )
    mu_0 = np.cos ( theta_0 )
    cos_gamma_down = np.cos ( theta_0 )*np.cos(theta) - np.sin(theta_0)*np.sin(theta)*np.cos(phi-pi_0)
    cos_gamma_up = - np.cos ( theta_0 )*np.cos(theta) - np.sin(theta_0)*np.sin(theta)*np.cos(phi-pi_0)
    
    P_gamma_down = (3./4.)*(1 + cos_gamma_down*cos_gamma_down )
    P_gamma_up = (3./4.) * (1 + (2*mu*mu_0 + cos_gamma_up)**2 )

    return P_gamma_down, P_gamma_up

def fresnel_reflectance ( theta, theta_0, m=1.396 ):
    mu = np.cos ( np.deg2rad ( theta ) )
    mu_0 = np.cos ( np.deg2rad ( theta_0 ) )
    
    y = (1./m)*np.sqrt ( m*m + mu*mu - 1 )
    y0 = (1./m)*np.sqrt ( m*m + mu0*mu0 - 1 )
    rho_mu = 1. - 2.*mu*y*m*( 1./(mu + m*y)**2 + 1./(m*mu*y)**2)
    rho_mu0 = 1. - 2.*mu0*y0*m*( 1./(mu0 + m*y0)**2 + 1./(m*mu0*y0)**2)
    return rho_mu, rho_mu0

def rayleigh_scattering ( theta, theta_0, phi, phi_0, lambdai, o3_conc, h0, doy ):
    # Calculate the terms in Eq. 2 of Wang et al
    # First, the relevant ETr... Note this is not implemented
    F0 = extraterrestrial_radiation ( doy, lambdai )
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
    return rayleigh_radiance
