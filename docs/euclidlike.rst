The Euclid-like Module
################################

The ``euclidlike`` module contains telescope information and functionality needed for image simulations.
The demo script end_to_end_demo.py shows how to use many of the atrributes and functions described here.


Module-level Attributes
=======================

There are several of attributes of the ``euclidlike`` module which define some numerical
parameters related to a Euclid-like geometry.  Some of these parameters relate to the entire
wide-field imager.  Others, especially the return values of the functions to get the
PSF and WCS, are specific to each CCD and therefore are indexed based on the detector number. 
All detector-related arrays are 0-indexed, which might differ from the CCD indices for Euclid,
which run 1-1 to n_row-n_col.

gain
    The gain for all CCDs is expected to be the roughly the same.

pixel_scale
    The pixel scale in units of arcsec/pixel. 

diameter
    The telescope diameter in meters.

obscuration
    The linear obscuration of the telescope, expressed as a fraction of the diameter.

collecting_area
    The actual collecting area after accounting for obscuration, struts, etc. in
        units of cm^2.

long_exptime
    The typical exposure time for the longer exposures used for VIS images, in units of seconds.  
    The number that is stored is for a single dither.

short_exptime_nisp
    The typical exposure time for the shorter NISP imaging exposures used for VIS images, in units of seconds.  
    The number that is stored is for a single dither.

short_exptime_vis
    The typical exposure time for the shorter exposures with VIS taken in parallel with NISP imaging, in units of seconds.  
    The number that is stored is for a single dither.

n_dithers
    The number of dithers per filter.

n_ccd
    The number of CCDs in the focal plane.
    
n_ccd_row
    The number of CCDs in the each focal plane row.

n_ccd_col
    The number of CCDs in the each focal plane column.

n_pix_row
    Each CCD has n_pix_row total pixels in each row.
n_pix_col
    Each CCD has n_pix_col total pixels in each column.

pixel_scale_mm
    The physical pixel size, in units of mm.

plate_scale
    The plate scale, in units of arcsec / mm

read_noise
    Total readout noise, in units of e-. 

saturation
    Pixel saturation, in units of e-. 

det2ccd
    Mapping from DETID to CCDID.

ccd2det
    Mapping from CCDID to DETID.

min_sun_angle
    Minimum allowed angle from the telescope solar panels to the sun, in degrees.

max_sun_angle
    Maximum allowed angle from the telescope solar panels to the sun, in degrees.

vis_bands
    List of available VIS bands
nisp_bands
    List of available NISP bands

vis_blue_limit
    Bandpass blue limit  needed for consistency with the wavelength range covered by ur tabulated 
    PSF images, in nm.

vis_red_limit
    Bandpass red limit  needed for consistency with the wavelength range covered by ur tabulated 
    PSF images, in nm.

For example, to get the gain value, use euclidlike.gain.


Euclid-like Functions
===============

This module also contains the following routines:

`euclidlike.getBandpasses`
    A utility to get a dictionary containing galsim.Bandpass objects for each of
    the Euclid-like imaging bandpasses, which by default have AB zeropoints given using
    the GalSim zeropoint convention (see `getBandpasses` docstring for more details).

`euclidlike.getSkyLevel`
    A utility to find the expected sky level due to zodiacal light at a given
    position, in a given band.

`euclidlike.getZodiBackground`
    This helper routine is enables the calculation of the zodiacal light, in photons/m^2/arcsec^2/sec

`euclidlike.getPSF`
    A routine to get a chromatic representation of the PSF in a single CCD.
    PSFs are based on precomputed, oversampled (by 3x) PSF images on a grid in wavelength and
    focal plane position.

`euclidlike.getBrightPSF`
    Get a fake optical PSF for very bright objects.

`euclidlike.getWCS`
    This routine returns a dict containing a WCS for each of the Euclid CCDs.

`euclidlike.findCCD`
    This is a helper routine to calculate the minimum and maximum pixel values that should be
    considered within a CCD.



.. autofunction:: euclidlike.getBandpasses

.. autofunction:: euclidlike.getSkyLevel

.. autofunction:: euclidlike.getZodiBackground

.. autofunction:: euclidlike.getPSF

.. autofunction:: euclidlike.getBrightPSF

.. autofunction:: euclidlike.getWCS

.. autofunction:: euclidlike.findCCD
