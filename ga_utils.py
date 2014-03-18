#!/usr/bin/env python

"""
Globalbedo prior utility functions
=====================================

A bunch of utility functions used in GA prior code to access
data etc.
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

# Authors etc
__author__ = "P Lewis & J Gomez-Dans (NCEO&UCL)"
__copyright__ = "(c) 2014"
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "J Gomez-Dans"
__email__ = "j.gomez-dans@ucl.ac.uk"
__status__ = "Development"

def locate( pattern, root=os.curdir ):
    """Locate all files matching supplied filename pattern in and below
    supplied root directory.
    Parameters
    ----------
    pattern: str
        A pattern in UNIX-speak e.g. "MCD43*.hdf"
    root: str
        The root directory where to start searching.
       """
    for path, dirs, files in os.walk( os.path.abspath( root ) ):
        for filename in fnmatch.filter( files, pattern ):
            yield os.path.join( path, filename )