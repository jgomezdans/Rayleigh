#!/usr/bin/env python

"""
Wang et al (2007) atmospheric correction for inland waters. We assume that this method will fail 
dramatically over turbid waters. Uses a lot of functions from `AtCorrUtils.py`
"""
import numpy as np
import matplotlib.pyplot as plt

from AtCorrUtils import *

def parse_command_line ( ):
    """Parse command line arguments"""
    import argparse
    
    parser = argparse.ArgumentParser( description='Inverting DALEC with' + \
        ' MCMC.\n By J Gomez-Dans, NCEO/UCL')
    parser.add_argument ( '-f', '--fname', action="store", required=True, \
            help="Input MTL filename" )
    parser.add_argument ( '-H', '--height', action="store", required=True, \
            help="Input height above seal level [km]" )

    args = vars(parser.parse_args())
    return ( args['fname'], float(args['height']) )

if __name__ == "__main__":
    # Parse the command line to get the MTL file
    fname, height = parse_command_line ()
    # Get acquisition geometry
    theta_i, phi_i = get_angles ( fname )
    # Get wavebands
    lambdas = get_lambdai ( fname )
    # Get time 
    doy, year = get_doy( fname )
    # Get O3 concentration
    o3_conc = get_o3conc ( doy, year )
    L_rayleigh = []
    tau_diff = []
    # Rayleigh scattering and O3 pissing 
    for lambdai in lambdas:
        lray, tau_diffu = rayleigh_scattering ( theta_i, 0.0, phi_i, 0.0, lambdai, o3_conc, height, doy )
        L_rayleigh.append ( lray )
        tau_diff.append ( tau_diffu )
    L_rayleigh = np.array ( L_rayleigh )
    tau_diff = np.array ( tau_diff )
    # Aerosol stuff starts here...
    water_leaving_radiance = water_leaving_rad_b7 ( tau_diff[-1], fname, L_rayleigh[-1] )