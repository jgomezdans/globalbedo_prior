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
    def __init__ ( self, tile, data_dir, output_dir, bands=[1,2,3,4,5,6,7] ):
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
        """
        self.tile = tile
        self.data_dir = data_dir
        self.output_dir = output_dir
        self.bands = bands
        # Now, get the list of filenames that we are going to use...
        self.get_modis_fnames ()
    
    def get_modis_fnames ( self ):
        """This method gets the MODIS filenames and stores them"""
        self.fnames_mcd43a1 = {}
        self.fnames_mcd43a2 = {}
        for doy in xrange ( 1, 367, 8 ):
            pattern = "MCD43A1.A????%03d.%s.005.*.hdf" % ( doy, self.tile )
            self.fnames_mcd43a1[doy] = [ f for f in locate( pattern, \
                root=self.data_dir ) ]
            pattern = "MCD43A2.A????%03d.%s.005.*.hdf" % ( doy, self.tile )
            self.fnames_mcd43a2[doy] = [ f for f in locate( pattern, \
                root=self.data_dir ) ]
            
if __name__ == "__main__":
    ga = GlobAlbedoPrior("h19v10", "/data/netapp_3/plewis/albedo/", "/data/netapp_3/plewis/albedo/prior", bands=[1,2] )
    