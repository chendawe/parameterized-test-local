#!/usr/bin/env python
"""
setup.py file for reverse chirp example pycbc waveform plugin package
"""

from setuptools import Extension, setup, Command
from setuptools import find_packages

VERSION = '0.0.0'

# The format for the entry_points is “capability”:”approximant_name = module_path:function_name”. The module path may include dots if the module is within a package or sub-package. The valid capbility is pycbc.waveform.td and pycbc.waveform.fd for time and frequency domain waveform models, respectively.

setup (
    name = 'dchi_modified_IMRPhenomXAS',
    version = VERSION,
    description = 'An IMRPhenomXAS-base non-GR waveform with PN coefficient modifiabilities.',
    long_description = open('descr.rst').read(),
    author = 'chendawww',
    author_email = '261364887@qq.com',
    # url = 'http://www.pycbc.org/',
    # download_url = 'https://github.com/gwastro/revchirp/tarball/v%s' % VERSION,
    keywords = ['pycbc', 'signal processing', 'gravitational waves'],
    install_requires = ['pycbc'],
    py_modules = ['dchi_modified_IMRPhenomXAS'],
    entry_points = {
        # "pycbc.waveform.td":"revchirp = revchirp:reverse_chirp_td",
                    "pycbc.waveform.fd":"dchi_modified_IMRPhenomXAS = dchi_modified_IMRPhenomXAS:get_dchi_modified_IMRPhenomXAS"},
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.6',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Scientific/Engineering :: Physics',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    ],
)