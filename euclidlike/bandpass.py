"""
@file bandpass.py

This file includes any routines needed to define the Euclid bandpasses.
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
from . import vis_red_limit, vis_blue_limit

def getBandpasses(AB_zeropoint=True, default_thin_trunc=True, full_bandpass=False, **kwargs):
    """
    Function to get the bandpass information for the Euclid VIS band and the three Euclid NISP passbands.

    This routine reads in files containing a list of wavelengths and
    transmission values for the Euclid bands. The files are located in the
    euclidlike.data directory. The routine then creates a Bandpass object
    using the LookupTable class from the GalSim package, and returns a dict with bandpasses for the
    keys.

    The bandpasses are publicly available from IPAC:
    http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?mode=browse&gname=Euclid&gname2=VIS&asttype=.
    https://euclid.esac.esa.int/msp/refdata/nisp/NISP-PHOTO-PASSBANDS-V1

    These are relatively old files that do not include the latest estimates of system response.
    They correspond to end-of-life estimates, with some expected degradation of the QE and filter
    transmission over time.  This can lead to flux estimates that are suppressed by 5-10% from
    beginning-of-life flux estimates.

    The VIS bandpass red and blue limits are set not by the transmission curve but by the range of
    wavelengths over which we have tabulated PSF images.  The wavelength range is read in from the
    instrument parameter file.

    Parameters:
        AB_zeropoint (bool) : If True, set the zeropoint of the bandpass to the AB magnitude system. [default: True]
        default_thin_trunc (bool) : If True, use the default thinning and truncation parameters. [default: True]
        full_bandpass (bool): if True, use the full bandpass without red/blue limits needed for PSF
                              calculations. [default: False]
        **kwargs : Additional keyword arguments to pass to either `Bandpass.thin` or `Bandpass.truncate`.
    
    Returns:
        A dict with bandpasses for the keys.
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

    # Then do NISP - make sure band names that are stored include 'NISP_'
    nisp_bands = ['Y', 'H', 'J']
    use_nisp_bands = []
    for band in nisp_bands:
        use_nisp_bands.append('NISP_%s'%band)
    for index, nisp_band in enumerate(nisp_bands):
        bandpass_file = files('euclidlike.data').joinpath('NISP-PHOTO-PASSBANDS-V1-%s_throughput.dat'%nisp_band)
        bandpass = np.loadtxt(bandpass_file, dtype=float)
        # These wavelengths are in nm but we want to be consistent for all bands, so multiply by 10
        # to get Angstroms
        wave[use_nisp_bands[index]] = bandpass[:,0]*10
        data[use_nisp_bands[index]] = bandpass[:,1]

    # make a list with all bands for later use.
    all_bands = use_nisp_bands.copy()
    all_bands.append('VIS')

        
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
        if bp_name == "VIS" and not full_bandpass:
            bp.blue_limit = vis_blue_limit
            bp.red_limit = vis_red_limit

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
    
