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
    print("You need to have numpy installed!")

try:
    from osgeo import gdal
except ImportError:
    print("You need to have the GDAL Python bindings installed!")

# Authors etc
__author__ = "P Lewis & J Gomez-Dans (NCEO&UCL)"
__copyright__ = "(c) 2014"
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "J Gomez-Dans"
__email__ = "j.gomez-dans@ucl.ac.uk"
__status__ = "Development"


GDAL2NUMPY = {
    gdal.GDT_Byte: np.uint8,
    gdal.GDT_UInt16: np.uint16,
    gdal.GDT_Int16: np.int16,
    gdal.GDT_UInt32: np.uint32,
    gdal.GDT_Int32: np.int32,
    gdal.GDT_Float32: np.float32,
    gdal.GDT_Float64: np.float64,
    gdal.GDT_CInt16: np.complex64,
    gdal.GDT_CInt32: np.complex64,
    gdal.GDT_CFloat32: np.complex64,
    gdal.GDT_CFloat64: np.complex128,
}


def calculate_prior(brdf_data, mask):
    """Calculates the prior mean from the data & data mask
    Prior is tested, and looks OK, the variance is untested"""
    prior_mean = np.zeros((3, brdf_data.shape[-2], brdf_data.shape[-1]))
    prior_var = np.zeros((3, brdf_data.shape[-2], brdf_data.shape[-1]))
    for i in range(3):
        A = np.ma.array(
            brdf_data[:, i, :, :] * 0.0010,
            mask=np.logical_or(brdf_data[:, i, :, :] == 32767, mask == 0),
        )
        kw_mean = np.ma.average(A, axis=0, weights=mask)
        v1 = np.ma.sum(mask, axis=0)
        v2 = np.ma.sum(mask ** 2, axis=0)
        kw_weight = np.sum(mask * (A - kw_mean) ** 2, axis=0) * (
            v1 / (v1 * v1 - v2)
        )
        prior_mean[i, :, :] = kw_mean
        prior_var[i, :, :] = kw_weight

    return prior_mean, prior_var


def locate(pattern, root=os.curdir):
    """Locate all files matching supplied filename pattern in and below
    supplied root directory.
    Parameters
    ----------
    pattern: str
        A pattern in UNIX-speak e.g. "MCD43*.hdf"
    root: str
        The root directory where to start searching.
    """
    import fnmatch

    for path, dirs, files in os.walk(os.path.abspath(root)):
        for filename in fnmatch.filter(files, pattern):
            yield os.path.join(path, filename)


def extract_chunks(the_files, the_bands=None):
    """An iterator function that extracts a chunk from one or several 
     GDAL-friendly datafiles. The files are stroed in ``the_files``, 
     and if they are multiband files, you can extract a data cube of 
     the required size by specifying it in ``the_bands`` (otherwise,
     by default all bands are returned. 
     
     This iterator returns several bits of information: a general 
     "configuration dictionary" that includes information on the size,
     projection, geotransform, the size of the extracted window and 
     its location, and any relevant data etc.
     
     Parameters
     -----------
     the_files: list
         A list of filenames to read from. GDAL friendly types.
     the_bands: list or None
         A list or array of bands to retrieve from the files above
         
     Returns
     --------
     ds_config; dict
         Returns a dictionary with the following useful keys: ``nx``,
         ``ny``, ``nb``, ``geoT`` and ``proj``
     this_X: int
         Starting point in X of the extracted window
     this_Y: int
         Starting point in Y of the extracted window
     nx_valid: int
         Size of the returned array in X
     ny_valid: int
         Size of the returned array in Y
     data_in: list
         A list of numpy arrays with the read data
     """
    ds_config = {}
    gdal_ptrs = []
    datatypes = []
    for the_file in the_files:
        g = gdal.Open(the_file)
        gdal_ptrs.append(gdal.Open(the_file))
        datatypes.append(GDAL2NUMPY[g.GetRasterBand(1).DataType])

    block_size = g.GetRasterBand(1).GetBlockSize()
    nx = g.RasterXSize
    ny = g.RasterYSize
    if the_bands is None:
        the_bands = np.arange(g.RasterCount) + 1
    proj = g.GetProjectionRef()
    geoT = g.GetGeoTransform()
    ds_config["nx"] = nx
    ds_config["ny"] = ny
    ds_config["nb"] = g.RasterCount
    ds_config["geoT"] = geoT
    ds_config["proj"] = proj
    block_size = [block_size[0], block_size[1]]
    # block_size = [ 1200, 1200 ]
    # store these numbers in variables that may change later
    nx_valid = block_size[0]
    ny_valid = block_size[1]
    # find total x and y blocks to be read
    nx_blocks = (int)((nx + block_size[0] - 1) / block_size[0])
    ny_blocks = (int)((ny + block_size[1] - 1) / block_size[1])
    buf_size = block_size[0] * block_size[1]
    ################################################################
    # start looping through blocks of data
    ################################################################
    # loop through X-lines
    for X in range(nx_blocks):
        # change the block size of the final piece
        if X == nx_blocks - 1:
            nx_valid = nx - X * block_size[0]
            buf_size = nx_valid * ny_valid

        # find X offset
        this_X = X * block_size[0]

        # reset buffer size for start of Y loop
        ny_valid = block_size[1]
        buf_size = nx_valid * ny_valid

        # loop through Y lines
        for Y in range(ny_blocks):
            # change the block size of the final piece
            if Y == ny_blocks - 1:
                ny_valid = ny - Y * block_size[1]
                buf_size = nx_valid * ny_valid

            # find Y offset
            this_Y = Y * block_size[1]
            data_in = []
            for ig, ptr in enumerate(gdal_ptrs):
                # We can have different numbers of bands!
                the_bands = np.arange(1, ptr.RasterCount + 1)
                buf = ptr.ReadRaster(
                    this_X,
                    this_Y,
                    nx_valid,
                    ny_valid,
                    buf_xsize=nx_valid,
                    buf_ysize=ny_valid,
                    band_list=the_bands,
                )
                a = np.frombuffer(buf, dtype=datatypes[ig])
                data_in.append(
                    a.reshape((len(the_bands), ny_valid, nx_valid)).squeeze()
                )

            yield (ds_config, this_X, this_Y, nx_valid, ny_valid, data_in)
