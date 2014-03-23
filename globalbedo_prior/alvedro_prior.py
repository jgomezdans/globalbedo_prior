#!/usr/bin/env python
# -*- coding: utf-8 -*-

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
 
"""bootstrap: sample command line application"""
 
from .GlobAlbedoPrior import GlobAlbedoPrior
 
# Authors etc
__author__ = "P Lewis & J Gomez-Dans (NCEO&UCL)"
__copyright__ = "(c) 2014"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "J Gomez-Dans"
__email__ = "j.gomez-dans@ucl.ac.uk"
__status__ = "Development"

def main():
    print("%s, version %s." % (__doc__, __version__))
    parser = optparse.OptionParser(formatter=optparse.TitledHelpFormatter(), \
        type=str, usage=globals()['__doc__'])
    parser.add_option ('-t', '--tile', action='store', dest="tile", \
        type=str, help='Tile to process (e.g. h17v04)' )
    parser.add_option ('-d', '--datadir', action='store', dest="datadir", \
        type=str, help='Root directory where data are stored' )
    parser.add_option ('-o', '--outputdir', action='store', dest="outdir", \
        type=str, help='Output directory' )
    parser.add_option ('-b', '--bands', action="store", dest="bands", 
        type=str, nargs="*", help="Bands to process, separated by a space" )

    (options, args) = parser.parse_args()
    ga = GlobAlbedoPrior( options.tile, options.datadir, \
         options.outdir, bands=options.bands )
    ga.stage1_prior()
    ga.stage2_prior()
 
if __name__ == "__main__":
    main()            
        
if __name__ == "__main__":
    ga.stage1_prior()