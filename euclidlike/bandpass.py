"""
@file bandpass.py

This file includes any routines needed to define the Euclid bandpass.
This module is heavily based on the roman bandpass.py file from the GalSim package.
https://github.com/GalSim-developers/GalSim/blob/releases/2.5/galsim/roman/roman_bandpass.py

The Euclid VIS bandpass is read in from the Euclid_VIS.vis.dat file which can be downloaded from
http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?mode=browse&gname=Euclid&gname2=VIS&asttype=.

The Euclid NISP bandpasses are read in from files downloaded from
https://euclid.esac.esa.int/msp/refdata/nisp/NISP-PHOTO-PASSBANDS-V1
"""

import numpy as np
import os
from galsim import Bandpass, LookupTable, galsim_warn
from importlib.resources import files


def getBandpasses(AB_zeropoint=True, default_thin_trunc=True, **kwargs):
    """
    Function to get the bandpass information for the Euclid VIS band the three Euclid NISP passbands.

    This routine reads in files containing a list of wavelengths and
    transmission values for the Euclid bands. The files are located in the
    euclidlike.data directory. The routine then creates a Bandpass object
    using the LookupTable class from the GalSim package.

    Args:
    AB_zeropoint (bool) : If True, set the zeropoint of the bandpass to the AB magnitude system. [default: True]
    default_thin_trunc (bool) : If True, use the default thinning and truncation parameters. [default: True]
    kwargs : Additional keyword arguments to pass to either `Bandpass.thin` or `Bandpass.truncate`.
    """
    # Read in the bandpass files, using a dict to distinguish the different filters
    wave = {}
    data = {}
    # Start with VIS
    bandpass_file = files('euclidlike.data').joinpath('Euclid_VIS.vis.dat')
    bandpass = np.loadtxt(bandpass_file, dtype=float)
    # Wavelengths in Angstroms
    wave['VIS'] = bandpass[:, 0]
    data['VIS'] = bandpass[:, 1]

    # Then do NISP
    nisp_bands = ['Y', 'H', 'J']
    all_bands = nisp_bands.copy()
    all_bands.append('VIS')
    for nisp_band in nisp_bands:
        bandpass_file = files('euclidlike.data').joinpath('NISP-PHOTO-PASSBANDS-V1-%s_throughput.dat.txt'%nisp_band)
        bandpass = np.loadtxt(bandpass_file, dtype=float)
        # These wavelengths are in nm but we want to be consistent for all bands, so multiply by 10
        # to get Angstroms
        wave[nisp_band] = bandpass[:,0]*10
        data[nisp_band] = bandpass[:,1]

    # Below is the original code from the GalSim package modified for the format of these Euclid
    # bandpass files. 

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
    for index, bp_name in enumerate(all_bands):
        # Create the bandpass object
        bp = Bandpass(LookupTable(wave[bp_name], data[bp_name]), wave_type='Angstrom')

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
    
