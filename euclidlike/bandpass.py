"""
@file bandpass.py

This file includes any routines needed to define the Euclid bandpass.
This module is heavily based on the roman bandpass.py file from the GalSim package.
https://github.com/GalSim-developers/GalSim/blob/releases/2.5/galsim/roman/roman_bandpass.py
The Euclid VIS bandpass is read in from the Euclid_VIS.vis.dat file which can be downloaded from
http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?mode=browse&gname=Euclid&gname2=VIS&asttype=.
"""

import numpy as np
import os
from galsim import Bandpass, LookupTable, galsim_warn
from importlib.resources import files


def getBandpasses(AB_zeropoint=True, default_thin_trunc=True, **kwargs):
    """
    Function to get the bandpass information for the Euclid VIS band.

    This routine reads in a file containing a list of wavelengths and
    transmission values for the Euclid VIS band. The file is located in the
    euclidlike.data directory. The routine then creates a Bandpass object
    using the LookupTable class from the GalSim package.

    Args:
    AB_zeropoint (bool) : If True, set the zeropoint of the bandpass to the AB magnitude system. [default: True]
    default_thin_trunc (bool) : If True, use the default thinning and truncation parameters. [default: True]
    kwargs : Additional keyword arguments to pass to either `Bandpass.thin` or `Bandpass.truncate`.

    TODO : Add the NISP bandpasses?
    """
    # Read in the bandpass file
    bandpass_file = files('euclidlike.data').joinpath('Euclid_VIS.vis.dat')
    bandpass = np.loadtxt(bandpass_file, dtype=float)
    wave = bandpass[:, 0]
    data = bandpass[:, 1]

    # In case we want to add NISP bandpasses
    data = np.atleast_2d(data)

    # Below is the original code from the GalSim package. 
    # I have modified it to read in the Euclid_VIS.vis.dat file.

    # Parse kwargs for truncation, thinning, etc., and check for nonsense.
    truncate_kwargs = ['blue_limit', 'red_limit', 'relative_throughput']
    thin_kwargs = ['rel_err', 'trim_zeros', 'preserve_range', 'fast_search']
    tmp_truncate_dict = {}
    tmp_thin_dict = {}
    if default_thin_trunc:
        if len(kwargs) > 0:
            galsim_warn('default_thin_trunc is true, but other arguments have been passed'
                        ' to getBandpasses().  Using the other arguments and ignoring'
                        ' default_thin_trunc.')
            default_thin_trunc = False
    if len(kwargs) > 0:
        for key in list(kwargs.keys()):
            if key in truncate_kwargs:
                tmp_truncate_dict[key] = kwargs.pop(key)
            if key in thin_kwargs:
                tmp_thin_dict[key] = kwargs.pop(key)
        if len(kwargs) != 0:
            raise TypeError("Unknown kwargs: %s"%(' '.join(kwargs.keys())))

    bandpass_dict = {}
    for index, bp_name in enumerate(['VIS']):
        # Create the bandpass object
        bp = Bandpass(LookupTable(wave, data[index]), wave_type='Angstrom')

        # Use any arguments related to truncation, thinning, etc.
        if len(tmp_truncate_dict) > 0 or default_thin_trunc:
            bp = bp.truncate(**tmp_truncate_dict)
        if len(tmp_thin_dict) > 0 or default_thin_trunc:
            bp = bp.thin(**tmp_thin_dict)

        # Set the zeropoint if requested by the user:
        if AB_zeropoint:
            bp = bp.withZeropoint('AB')

        bp.name = bp_name
        bandpass_dict[bp.name] = bp

    return bandpass_dict
    
