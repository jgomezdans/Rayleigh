#!/usr/bin/env python

"""
Wang et al (2007) atmospheric correction for inland waters. We assume that this method will fail 
dramatically over turbid waters. Uses a lot of functions from `AtCorrUtils.py`
"""
import numpy as np
import matplotlib.pyplot as plt


def parse_command_line ( ):
    """Parse command line arguments"""
    import argparse
    
    parser = argparse.ArgumentParser( description='Inverting DALEC with' + \
        ' MCMC.\n By J Gomez-Dans, NCEO/UCL')
    parser.add_argument ('-a','--autoregressive', action="store_true", \
            default=False, help="Use an OU likelihood formalism" )
    parser.add_argument('-i','--init_file', action="store", help='Init file', \
            required=True)
    parser.add_argument('-n','--num_years', help='How many years to invert', \
            required=True, type=int, action="store" )
    parser.add_argument ('-t','--truncate', action="store_true", \
            default=False, help="Truncate prior")
    args = vars(parser.parse_args())
    return ( args['autoregressive'], args['init_file'], args['num_years'], \
            args['truncate'] )

if __name__ == "__main__":
    