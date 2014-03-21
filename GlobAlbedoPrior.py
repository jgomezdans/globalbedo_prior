#!/usr/bin/env python

"""
Globalbedo prior generation
============================

This script produces the GlobAlbedo prior from MODIS MCD43 data.
The prior is calculated in two stages: a first, single date
stage, which provides a weighted average of the kernel weights
for a particular day and a second stage that performs some temporal
smoothing and interpolates the result to daily. Two versions of the
prior are produced: a snow and a snow-free version.

"""
import argparse
import sys
import os
import calendar
import logging
import time

try:
    import numpy as np
except ImportError:
    print "You need to have numpy installed!"

try:
    from osgeo import gdal
except ImportError:
    print "You need to have the GDAL Python bindings installed!"

from ga_utils import *

# Authors etc
__author__ = "P Lewis & J Gomez-Dans (NCEO&UCL)"
__copyright__ = "(c) 2014"
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "J Gomez-Dans"
__email__ = "j.gomez-dans@ucl.ac.uk"
__status__ = "Development"

# Set up logging
LOG = logging.getLogger( __name__ )
OUT_HDLR = logging.StreamHandler( sys.stdout )
OUT_HDLR.setFormatter( logging.Formatter( '%(asctime)s %(message)s') )
OUT_HDLR.setLevel( logging.INFO )
LOG.addHandler( OUT_HDLR )
LOG.setLevel( logging.INFO )

class GlobAlbedoPrior ( object ):
    """
    """
    def __init__ ( self, tile, data_dir, output_dir, bands=[1,2,3,4,5,6,7], no_snow=True ):
        """This configures where the data are stored, the tile and the output
        directory. It also allows the choice of bands to process, which by
        default is the 7 MODIS bands.

        Parameters
        ----------
        tile: str 
            The tile in MODIS-speak, e.g. "h17v04"
        data_dir: str
            The top directory where input data are stored
        output_dir: str
            The directory where the output will be saved to
        bands: array-like
            The set of bands to be processed. By default, all 7 bands.
        no_snow: bool
            A flag indicating whether to choose the snow free or the snow albedo results
        """
        self.tile = tile
        self.data_dir = data_dir
        self.output_dir = output_dir
        self.bands = bands
        # Now, get the list of filenames that we are going to use...
        self.get_modis_fnames ()
        self.no_snow = no_snow
    
    def get_modis_fnames ( self ):
        """This method gets the MODIS filenames and stores them in dictionaries
        that can be accessed by DoY. This is quite useful for Stage 1 prior
        creation.
        """
        self.fnames_mcd43a1 = {}
        self.fnames_mcd43a2 = {}
        for doy in xrange ( 1, 367, 8 ):
            pattern = "MCD43A1.A????%03d.%s.005.*.hdf" % ( doy, self.tile )
            self.fnames_mcd43a1[doy] = [ f for f in locate( pattern, \
                root=self.data_dir ) ]
            pattern = "MCD43A2.A????%03d.%s.005.*.hdf" % ( doy, self.tile )
            self.fnames_mcd43a2[doy] = [ f for f in locate( pattern, \
                root=self.data_dir ) ]
            
    def _interpret_qa ( self, qa_data ):
        """A method to interpret the QA according to the GlobAlbedo specs"""
        # First stage: consider the QA flag, first 3 bits
        qa = np.bitwise_and ( qa_data, 7L ) # 7 is b000111, ie get the last 3 bits
        weight = np.where ( qa < 4, MAGIC**qa, 0. )
        return weight
    
    def _interpret_snowmask ( self, snow_data ):
        """Interpret the snow mask"""
        
        if self.no_snow:
            snow = np.where ( snow_data == 0, True, False ) # Snow free
        else:
            snow = np.where ( snow_data == 1, True, False ) # Snow 
        return snow
            
    def _interpret_landcover ( self, landcover ):
        """Interpret the landcover information"""
        # Landcover uses bits 04-07, so 8 + 16 + 32 + 64 = 120L
        # then shift back 3 positions
        mask = np.right_shift(np.bitwise_and ( landcover, 120L), 3)
        mask = np.where ( mask < 7, True, False ) # ignore deep ocean
        return mask
    
    def create_output ( self ):
        self.output_ptrs = {}
        gdal_opts = [ "COMPRESS=LZW", "INTERLEAVE=BAND", "TILED=YES" ]
        for band in self.bands:
            for doy in self.fnames_mcd43a1.iterkeys():
                output_fname = os.path.join ( root_dir, "%s.%s.%s.tif" % \
                            ( output, product, suffix ) )
                                                      
                if os.path.exists ( output_fname ):
                    print "Removing %s... " % output_fname, 
                    os.remove ( output_fname )
                    sys.stdout.flush()
                print "Creating %s" % output_fname,
                sys.stdout.flush()
                drv = gdal.GetDriverByName ( "GTiff" )
                output_prod = "%s_%s" % ( output, product )
                self.output_ptrs[output_prod] = drv.Create( output_fname, 2400, 2400, \
                    3, gdal.GDT_Float32, options=gdal_opts )
                
                print "... Created!"
    def stage1_prior ( self, band ):
        """Produce the stage 1 prior, which is simply a weighted average of the
        kernel weights. This"""
        
        
        for doy in self.fnames_mcd43a1.iterkeys():
            obs_fnames = [ 'HDF4_EOS:EOS_GRID:"%s":MOD_Grid_BRDF:BRDF_Albedo_Parameters_Band%d' % ( f, band ) \
                    for f in self.fnames_mcd43a1[doy] ]
            qa_fnames = [ 'HDF4_EOS:EOS_GRID:"%s":MOD_Grid_BRDF:BRDF_Albedo_Band_Quality' % ( f ) \
                    for f in self.fnames_mcd43a2[doy] ]
            snow_fnames = [ 'HDF4_EOS:EOS_GRID:"%s":MOD_Grid_BRDF:Snow_BRDF_Albedo' % ( f ) \
                    for f in self.fnames_mcd43a2[doy] ]
            land_fnames = [ 'HDF4_EOS:EOS_GRID:"%s":MOD_Grid_BRDF:BRDF_Albedo_Ancillary' % ( f ) \
                    for f in self.fnames_mcd43a2[doy] ]


            all_files = obs_fnames + qa_fnames
            first_time = True
            for (ds_config, this_X, this_Y, nx_valid, ny_valid, data_in )  \
                    in extract_chunk ( all_files ):
                n_years = len ( data_in )/4 # First lot will be BRDF parameters, 2nd will be QA,
                                            # 3 snow mask and fourth land mask
                brdf_data = data_in [ :n_years]
                
                self._interpret_qa ( qa_data )
                if first_time:
                    # set geotransforme etc for this glorious day
                    self.output_ptrs[output_prod].SetGeoTransform( ds_config['geoT'] )
                    self.output[output_prod].SetProjection( ds_config['proj'] )
                    first_time = False
                # data_in contains all the data. Half the samples bands are BRDF parameters        
            
if __name__ == "__main__":
    ga = GlobAlbedoPrior("h19v10", "/data/netapp_3/plewis/albedo/", "/data/netapp_3/plewis/albedo/prior", bands=[1,2] )
    