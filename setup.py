# -*- coding: utf-8 -*-
 
import re
from setuptools import setup
 
version = re.search(
    '^__version__\s*=\s*"(.*)"',
    open('globalbedo_prior/alvedro_prior.py').read(),
    re.M
    ).group(1)
 
setup(
    name = "globalbedo_prior",
    packages = ["globalbedo_prior"],
    entry_points = {
        "console_scripts": ['globalbedo_prior = globalbedo_prior.main:main']
        },
    version = version,
    description = "A Python implementation of the GlobAlbedo prior calculations.",
    author = "J Gomez-Dans",
    author_email = "j.gomez-dans@ucl.ac.uk",
    url = "http://github.com/jgomezdans/globalbedo_prior/",
 
    )