globalbedo_prior
================

Code to produce the GlobAlbedo prior from MODIS data, using the MCD43A1 and MCD43A2 MODIS products.

The Python package can be installed using pip or easy_install (it's available from `the cheese shop <https://pypi.python.org/pypi/globalbedo_prior>`_. Once installed with e.g. ``pip install globalbedo_prior --user --upgrade``, a script called ``alvedro_prior`` should appear in your path (this depends on pip/easy_intall doing their work properly). This script can be used to produce a daily prior of BRDF kernel parameters derived from the MODIS products.

The generation of the prior is very simple, and is performed in two stages:

**Stage 1**
    For each pixel, the entire timeseries of 8-daily observations are averaged using a weight derived from the QA flags. This results in a mostly complete 8-daily kernel product, stored as GeoTIFF files. Note that in some regions with cloud problems, there can be empty pixels as no observations are available for the period of interest within the MODIS record. Note that we calculate both the mean and standar deviation of the kernel parameters.
    
**Stage 2**
    For each day of year and pixel, a simple Laplacian filter in time is used to interpolate temporally. The filter is quite peaky, and its weight decays to 0.5 8 days after the day of interest.
    

Usage
------

The usage using the ``alvedro_prior`` script is very simple: just stash the MCD43A1 and MCD43A2 files somewhere (no need for fancy directories or anything, although that helps you!), and select a tile. Then decide whether the output will be saved to, and execute a command like this:

.. code::

   nohup alvedro_prior --tile h17v04 --datadir <my_data_directory_root> --outdir <my_output_directory> &
   
The previous command will search for the MCD43A1/2 files below ``<my_data_directory_root>`` that relate to tile ``h17v04`` and save the Stage 1 and Stage 2 priors in ``<my_output_directory>``.


The data
---------

Stage 1 priors are written for each 8 day period in the year, for the three kernels and have two bands: the mean and the "weight" (or inverse of the variance). Stage 2 priors are given per kernel, and have 366 bands, each of them with the prior mean for that particular day. The uncertainty associated with the Stage 2 prior has not been calculated.