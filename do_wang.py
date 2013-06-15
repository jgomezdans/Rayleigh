#!/usr/bin/env python
"""
SYNOPSIS

    do_wang.py [-h,--help] [-v,--verbose] [-i,--input] [-O,--o3] [-H,--height] [-r,--roi]
    

DESCRIPTION

  This program receives a metadata file in USGS format, and it then applies
  the Wang atmospheric correction method based on a minimimal set of input
  parameters. It will take the original data, subset it if required using a
  geographical box in projection units (UTM). It will produce datasets with
  the TOA radiance (e.g. `LE72040312010347EDC00_ROI_B2_TOARAD.tif`), as well
  as a datafile with the visible bands with atmospheric correction  (filename
  is e.g. `LE72040312010347EDC00_WANG_VIS_WLRAD.tif`). 

EXAMPLES
            $ ./do_atcorr.py -H 0.5 -i data/LE72040312010347EDC00_MTL.txt \
                --roi 578745.000,4650765.000,608535.000,4618935.000 -v
                
            Sat Jun 15 15:19:56 2013
            1 data/LE72040312010347EDC00_ROI_B1.vrt data/LE72040312010347EDC00_B1.TIF
            Input file size is 8081, 7151
            Computed -srcwin 3822 2685 993 1061 from projected window.
            2 data/LE72040312010347EDC00_ROI_B2.vrt data/LE72040312010347EDC00_B2.TIF
            Input file size is 8081, 7151
            Computed -srcwin 3822 2685 993 1061 from projected window.
            3 data/LE72040312010347EDC00_ROI_B3.vrt data/LE72040312010347EDC00_B3.TIF
            Input file size is 8081, 7151
            Computed -srcwin 3822 2685 993 1061 from projected window.
            4 data/LE72040312010347EDC00_ROI_B4.vrt data/LE72040312010347EDC00_B4.TIF
            Input file size is 8081, 7151
            Computed -srcwin 3822 2685 993 1061 from projected window.
            5 data/LE72040312010347EDC00_ROI_B5.vrt data/LE72040312010347EDC00_B5.TIF
            Input file size is 8081, 7151
            Computed -srcwin 3822 2685 993 1061 from projected window.
            7 data/LE72040312010347EDC00_ROI_B7.vrt data/LE72040312010347EDC00_B7.TIF
            Input file size is 8081, 7151
            Computed -srcwin 3822 2685 993 1061 from projected window.
            Start reading data/LE72040312010347EDC00_ROI_B1.vrt
            Using blocksize 128 x 128
            Creating output data/LE72040312010347EDC00_ROI_B1_TOARAD.tif
            Start reading data/LE72040312010347EDC00_ROI_B2.vrt
            Using blocksize 128 x 128
            Creating output data/LE72040312010347EDC00_ROI_B2_TOARAD.tif
            Start reading data/LE72040312010347EDC00_ROI_B3.vrt
            Using blocksize 128 x 128
            Creating output data/LE72040312010347EDC00_ROI_B3_TOARAD.tif
            Start reading data/LE72040312010347EDC00_ROI_B4.vrt
            Using blocksize 128 x 128
            Creating output data/LE72040312010347EDC00_ROI_B4_TOARAD.tif
            Start reading data/LE72040312010347EDC00_ROI_B5.vrt
            Using blocksize 128 x 128
            Creating output data/LE72040312010347EDC00_ROI_B5_TOARAD.tif
            Start reading data/LE72040312010347EDC00_ROI_B7.vrt
            Using blocksize 128 x 128
            Creating output data/LE72040312010347EDC00_ROI_B7_TOARAD.tif
                    Theta_i=22.545988, Phi_i=160.273075
            Lambdas:  [   482.5    565.     660.     837.5   1650.   11450.    2220. ]
            LE72040312010347EDC00_MTL.txt
            Doy: 347, Year: 2010
            Using default O3 conc file, O3 conc: 265.000000
            Starting interpolation...
            Interpolation done...
            Creating output data/LE72040312010347EDC00_WANG_VIS_WLRAD.tif
            Sat Jun 15 15:23:22 2013
            TOTAL TIME IN MINUTES: 3.43263668219


EXIT STATUS

    -1 if numpy and/or GDAL aren't present

AUTHOR

    J Gomez-Dans <j.gomez-dans@ucl.ac.uk>

NOTE
    
    The program has not been verified. Needs severe testing!!!!!!!!
"""

    


import sys
import os
import traceback
import optparse
import time

try:
    import numpy as np
except ImportError:
    print "You need numpy installed"
    sys.exit ( -1 )
try:
  from osgeo import gdal
except ImportError:
    print "You need GDAL installed"
    sys.exit ( -1 )
    
    
from AtCorrUtils import *
        
GDAL_OPTS = ["COMPRESS=LZW", "INTERLEAVE=PIXEL", "TILED=YES",\
        "SPARSE_OK=TRUE", "BIGTIFF=YES" ]

def process_metadata ( fname ):
    """A function to extract the relelvant metadata from the
    USGS control file. Returns dicionaries with LMAX, LMIN,
    QCAL_LMIN and QCAL_LMAX for each of the bands of interest."""

    fp = open( fname, 'r') # Open metadata file
    lmax = {} # Dicts to store constants
    lmin = {}
    qc_lmax = {}
    qc_lmin = {}
    gain = {}
    bias = {}

    for line in fp: # 
      # Check for LMAX and LMIN strings
      # Note that parse logic is identical to the first case
      # This version of the code works, but is rather inelegant!
      if ( line.find ("RADIANCE_MULT_BAND") >= 0 ):
          s = line.split("=") # Split by equal sign
          the_band = int(s[0].split("_")[3]) # Band number as integer
          if the_band in [1,2,3,4,5,7]: # Is this one of the bands we want?
              gain[the_band] = float ( s[-1] ) # Get constant as float
      elif ( line.find ("RADIANCE_ADD_BAND") >= 0 ):
          s = line.split("=") # Split by equal sign
          the_band = int(s[0].split("_")[3]) # Band number as integer
          if the_band in [1,2,3,4,5,7]: # Is this one of the bands we want?
              bias[the_band] = float ( s[-1] ) # Get constant as float
      elif ( line.find ("QUANTIZE_CAL_MAX_BAND") >= 0 ):
          s = line.split("=") # Split by equal sign
          the_band = int(s[0].split("_")[4]) # Band number as integer
          if the_band in [1,2,3,4,5,7]: # Is this one of the bands we want?
              qc_lmax[the_band] = float ( s[-1] ) # Get constant as float
      elif ( line.find ("QUANTIZE_CAL_MIN_BAND") >= 0 ):
          s = line.split("=") # Split by equal sign
          the_band = int(s[0].split("_")[4]) # Band number as integer
          if the_band in [1,2,3,4,5,7]: # Is this one of the bands we want?
              qc_lmin[the_band] = float ( s[-1] ) # Get constant as float
      elif ( line.find ("RADIANCE_MAXIMUM_BAND") >= 0 ):
          s = line.split("=") # Split by equal sign
          the_band = int(s[0].split("_")[3]) # Band number as integer
          if the_band in [1,2,3,4,5,7]: # Is this one of the bands we want?
              lmax[the_band] = float ( s[-1] ) # Get constant as float
      elif ( line.find ("RADIANCE_MINIMUM_BAND") >= 0 ):
          s = line.split("=") # Split by equal sign
          the_band = int(s[0].split("_")[3]) # Band number as integer
          if the_band in [1,2,3,4,5,7]: # Is this one of the bands we want?
              lmin[the_band] = float ( s[-1] ) # Get constant as float

    return ( bias, gain, lmax, lmin, qc_lmax, qc_lmin )



def get_metadata ( fname ):
    """This function takes `fname`, a filename (opionally with a path), and
    and works out the associated metadata file"""
    original_fname = os.path.basename ( fname )
    metadata_fname = original_fname.split("_")[0] + "_MTL.txt"
    metadata_fname = os.path.join ( os.path.dirname ( fname ),metadata_fname )
    return metadata_fname

def extract_chunk ( the_file, n_blocks=1 ):
    
    """A function that extracts a chunk from a datafile"""
    ds_config = {}
    g = gdal.Open ( the_file )
    block_size = g.GetRasterBand(1).GetBlockSize()
    nx = g.RasterXSize
    ny = g.RasterYSize
    the_bands = np.arange ( g.RasterCount ) + 1
    proj = g.GetProjectionRef()
    geoT = g.GetGeoTransform()
    ds_config['nx'] = nx
    ds_config['ny'] = ny
    ds_config['nb'] = g.RasterCount
    ds_config['geoT'] = geoT
    ds_config['proj'] = proj
    block_size = [ block_size[0]*n_blocks, block_size[1]*n_blocks ]
    # store these numbers in variables that may change later
    nx_valid = block_size[0]
    ny_valid = block_size[1]
    # find total x and y blocks to be read
    nx_blocks = (int)((nx + block_size[0] - 1) / block_size[0]);
    ny_blocks = (int)((ny + block_size[1] - 1) / block_size[1]);
    buf_size = block_size[0]*block_size[1]
    print("Using blocksize %s x %s" %(block_size[0], block_size[1]))
    sys.stdout.flush()
    ################################################################
    # start looping through blocks of data
    ################################################################
    # loop through X-lines
    for X in xrange( nx_blocks ):
        # change the block size of the final piece
        if X == nx_blocks-1:
             nx_valid = nx - X * block_size[0]
             buf_size = nx_valid*ny_valid

        # find X offset
        this_X = X*block_size[0]

        # reset buffer size for start of Y loop
        ny_valid = block_size[1]
        buf_size = nx_valid*ny_valid

        # loop through Y lines
        for Y in xrange( ny_blocks ):
            # change the block size of the final piece
            if Y == ny_blocks-1:
                ny_valid = ny - Y * block_size[1]
                buf_size = nx_valid*ny_valid

            # find Y offset
            this_Y = Y*block_size[1]
            
            buf = g.ReadRaster(this_X, this_Y, nx_valid, ny_valid, \
                buf_xsize=nx_valid, buf_ysize=ny_valid,  \
                band_list= the_bands )
            a = np.frombuffer(buf, dtype=np.uint8 )
            if len(the_bands) > 1:
                data_in = a.reshape(( len(the_bands), ny_valid, nx_valid))
            else:
                data_in = a.reshape(( ny_valid, nx_valid))

            yield ( ds_config, this_X, this_Y, nx_valid, ny_valid, data_in )

def convert_to_radiance ( fname, band, gain, bias, lmax, lmin, qc_lmax, qc_lmin, verbose ):
    """
    This function reads the input file chunk by chunk using ``extract_chunk``,
    and converts from DN to radiance using the methodology given by 
    <http://landsat.usgs.gov/how_is_radiance_calculated.php>. 
    """
    # This variable is used to create the output dataset when the first 
    # chunk of data is read.
    first_time = True
    # `extract_chunk` can be used to efficiently read chunks of data
    # from a GDAL dataset. It uses a "generator", which means that we
    # can iterate over it using e.g. a for loop
    if verbose:
        print "Start reading %s" % fname
    # The output filaname.
    if fname.find ( ".TIF" ) >= 0:
        output_fname = fname.replace(".TIF", "_TOARAD.tif" )
    elif fname.find (".vrt" ) >= 0:
        output_fname = fname.replace(".vrt", "_TOARAD.tif" )
    for ( ds_config, this_X, this_Y, nx_valid, ny_valid, data_in ) in \
        extract_chunk ( fname ):
        if first_time:
            if verbose:
                print "Creating output %s" % output_fname
            # Create output dataset if `first_time` is true
            drv = gdal.GetDriverByName ( "GTiff" )
            dst_ds = drv.Create ( output_fname, ds_config['nx'], \
                ds_config['ny'], 1, gdal.GDT_Float32, options=GDAL_OPTS )
            dst_ds.SetGeoTransform ( ds_config['geoT'] )
            dst_ds.SetProjection ( ds_config['proj'] )
            first_time = False
        
        # Create the output array. Not needed, but we can
        # conveniently set the output type here to float32
        
        radiance = np.zeros_like ( data_in, dtype=np.float32 )
        # DN to radiance conversion if we have a sensible DN
        passer = np.logical_and ( qc_lmin < data_in, data_in < qc_lmax )
        radiance = np.where ( passer, \
             (( lmax - lmin )/(qc_lmax-qc_lmin))*\
             ( (data_in) - qc_lmin) + lmin, \
              -999 )
        # Write the output chunk
        dst_ds.GetRasterBand ( 1 ).WriteArray( radiance, \
                    xoff=this_X, yoff=this_Y)
    # Add some useful metadata to the dataset
    dst_ds.GetRasterBand( 1 ).SetNoDataValue ( -999 )
    dst_ds.GetRasterBand( 1 ).SetMetadata ( {"Band": "%d" % band, \
                                    "Units": "W/(m2*ster*um)", \
                                    "Data": "TOA Radiance" } )
    # Need to do this to flush the dataset to disk
    dst_ds = None
    
def main ( subsets = None ):
    """The main function""" 
  

    global options
    global args
    
    metadata_file = get_metadata ( options.input_f )
    ( gain, bias, lmax, lmin, qc_lmax, qc_lmin ) = process_metadata ( metadata_file )
    original_fname = os.path.basename ( options.input_f)
    prefix = original_fname.split("_")[0] 
    fname_prefix = os.path.join ( os.path.dirname ( options.input_f ), prefix )
    for (i, the_band) in enumerate ( [1,2,3,4,5,7] ):
        if subsets is None:
            input_file = fname_prefix + "_B%d.TIF" % the_band
        else:
            input_file = subsets[i]
        convert_to_radiance ( input_file, the_band, \
                    gain[the_band], bias[the_band], \
                    lmax[the_band], lmin[the_band], qc_lmax[the_band], \
                    qc_lmin[the_band], options.verbose )
    
def subset_datasets ( metadata_file, cutline ):
    """Subset origianl datafiles for given ROI
    
    This function subsets the original datafiles according to a
    user-provided ROI in raster coordinates (typical UTM, units m).
    """
    import subprocess
    ulx, uly, lrx, lry = [ float(x) for x in cutline.split(",") ]
    dn_fnames = []
    for band in [ 1, 2, 3, 4, 5, 7 ]:
        print band,
        new_fname = metadata_file.replace("_MTL.txt", "_ROI_B%d.vrt" % band )
        input_fname = metadata_file.replace("_MTL.txt", "_B%d.TIF" % band )
        print new_fname, input_fname
        retval = subprocess.call ("gdal_translate -of VRT -projwin %f %f %f %f %s %s" % \
                ( ulx, uly, lrx, lry, input_fname, new_fname ), shell=True )
        if retval != 0:
            print "Problem selecting ROI %f %f %f %f in %s" % \
                ( ulx, uly, lrx, lry, input_fname)
            raise IOError
        dn_fnames.append ( new_fname )
    return dn_fnames
def do_wang_atcorr ( fname, height, o3_conc, verbose=False ):
    
    
    # Get acquisition geometry
    theta_i, phi_i = get_angles ( fname )
    if verbose:
        print "\tTheta_i=%f, Phi_i=%f" % ( theta_i, phi_i )
    # Get wavebands
    lambdas = get_lambdai ( fname )
    if verbose:
        print "Lambdas: ", lambdas
    # Get time 
    doy, year = get_doy( fname )
    if verbose:
        print "Doy: %d, Year: %d" % ( doy, year )
    # Get O3 concentration
    if o3_conc is None:
        o3_conc = get_o3conc ( doy, year )
        if verbose:
            print "Using default O3 conc file, O3 conc: %f" % o3_conc
                        
        
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
    water_leaving_radiance_vis = aerosol_correction ( tau_diff, fname, L_rayleigh, doy, \
        lambdas, theta_i, verbose=verbose )
    output_fname = fname.replace("MTL.txt", "WANG_VIS_WLRAD.tif" )
    if verbose:
        print "Creating output %s" % output_fname
    afname = fname.replace ( "MTL.txt", "B1_TOARAD.tif" )
    if not os.path.exists ( afname ):
        afname = fname.replace ( "MTL.txt", "ROI_B1_TOARAD.tif"  )
    g = gdal.Open ( afname )
    if g is None:
        raise IOError, "Can't find file %s" % fname

    # Create output dataset if `first_time` is true
    drv = gdal.GetDriverByName ( "GTiff" )
    dst_ds = drv.Create ( output_fname, g.RasterXSize, \
                g.RasterYSize, 3, gdal.GDT_Float32, options=GDAL_OPTS )
    dst_ds.SetGeoTransform ( g.GetGeoTransform() )
    dst_ds.SetProjection ( g.GetProjectionRef() )
    for band in xrange(3):
        dst_ds.GetRasterBand(band+1).WriteArray ( water_leaving_radiance_vis[band, :, :] )
    dst_ds = None

if __name__ == '__main__':

    start_time = time.time()
    parser = optparse.OptionParser(formatter=optparse.TitledHelpFormatter(), \
            usage=globals()['__doc__'])
    parser.add_option ('-v', '--verbose', action='store_true', \
            default=False, help='verbose output')
    parser.add_option ('-i', '--input', action='store', dest="input_f",\
            type=str, help='Input USGS metadata filename')
    parser.add_option ('-O', '--o3', action='store', dest="o3_conc", \
            type=str, help="O3 concentration" )
    parser.add_option ('-H', '--height', action="store", dest="height", \
            type=float, help="Scene height above mean sea level (km)")
    parser.add_option('-r', '--roi', action="store", dest="cut", \
            type=str, help="ROI specification (ulx uly lrx lry)" )
        
    (options, args) = parser.parse_args()
    if options.verbose: print time.asctime()
    if options.cut is not None:
        the_fnames = subset_datasets ( options.input_f, options.cut )
        main ( subsets = the_fnames )
    else:
        main()
    do_wang_atcorr ( options.input_f, options.height, options.o3_conc, \
        verbose = options.verbose )
    if options.verbose: print time.asctime()
    if options.verbose: print 'TOTAL TIME IN MINUTES:',
    if options.verbose: print (time.time() - start_time) / 60.0
