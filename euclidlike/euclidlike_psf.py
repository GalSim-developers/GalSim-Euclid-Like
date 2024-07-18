import galsim
import numpy as np
import astropy.io.fits as pyfits
from importlib.resources import files

from . import n_ccd, n_pix_col, n_pix_row, pixel_scale, diameter, obscuration
from .bandpass import getBandpasses

"""
@file euclidlike_psf.py

Part of the Euclid-like simulation module. This file includes routines needed
to define a Euclid-like PSF.
"""


def get_euclid_wavelength_psf():
    """ This function get the oversampled PSF image as a function of wavelength
    from Lance as an input the Euclid-like simulation
    """
    # NOTE: We do not have PSF variation
    psf_file = files("euclidlike.data").joinpath("monopsfs_6_6.fits.gz")
    image_array = pyfits.getdata(psf_file)
    # get wavelengths directly from data
    wave_data = pyfits.getdata(psf_file, 1)
    # factor of 1e3 to convert from micrometer to nm
    wave_list = np.hstack(wave_data)*1e3
    # The following are the values for the data from Lance Miller
    nsample = len(wave_list)
    scale = pixel_scale/3  # images are oversampled by a factor of 3
    im_list = []
    for i in range(nsample):
        im_list.append(
            galsim.Image(image_array[i], scale=scale)
        )
    return wave_list, im_list


def getPSF(
        ccd, bandpass,
        ccd_pos=None, wcs=None,
        wavelength=None, gsparams=None,
        logger=None
):
    """Get a single PSF for Euclid like simulation.

    For applications that require very high accuracy in the modeling of the
    PSF, with very limited aliasing, you may want to lower the
    folding_threshold in the gsparams.  Otherwise very bright stars will show
    some reflections in the spider pattern and possibly some boxiness at the
    outskirts of the PSF.  Using ``gsparams =
    GSParams(folding_threshold=1.e-4)`` generally provides good results.

    Args:
    ccd (int):  Single value specifying the ccd for which the PSF should be
        loaded.
    bandpass (str): Single string specifying the bandpass to use when
        defining the pupil plane configuration and/or interpolation of
        chromatic PSFs.
    ccd_pos:  Single galsim.PositionD indicating the position within the ccd
        for which the PSF should be created. If None, the exact center of the
        ccd is chosen. [default: None]
    wcs:  The WCS to use to project the PSF into world coordinates. [default:
        galsim.PixelScale(euclidlike.pixel_scale)]
    wavelength (float):  An option to get an achromatic PSF for a single
        wavelength, for users who do not care about chromaticity of the PSF. If
        None, then the fully chromatic PSF is returned as an
        InterpolatedChromaticObject. Alternatively the user should supply
        either (a) a wavelength in nanometers, and they will get an
        InterpolatedImage object for that wavelength, or (b) a bandpass object,
        in which case they will get an InterpolatedImage objects defined at the
        effective wavelength of that bandpass. [default: None]
    gsparams:  An optional GSParams argument.  See the docstring for GSParams
        for details. [default: None]

    Returns:
        A single PSF object (either an InterpolatedChromaticObject or an
        InterpolatedImage depending on the inputs).

    """

    if ccd < 0 or ccd >= n_ccd:
        raise galsim.GalSimRangeError("Invalid ccd.", ccd, 0, n_ccd-1)

    assert bandpass == "VIS", "Only VIS band is supported"

    # ccd_pos: if None, then all should just be center of the ccd.
    if ccd_pos is None:
        ccd_pos = galsim.PositionD(x=n_pix_col/2, y=n_pix_row/2)

    if not isinstance(wavelength, (galsim.Bandpass, float, type(None))):
        raise TypeError(
            "wavelength should either be a Bandpass, float, or None."
        )

    # Now get psf model
    psf = _get_single_psf_obj(ccd, bandpass, ccd_pos, wavelength, gsparams)
    # Apply WCS.
    # The current version is in arcsec units, but oriented parallel to the
    # image coordinates. So to apply the right WCS, project to pixels using the
    # Euclid mean pixel_scale, then project back to world coordinates with the
    # provided wcs.
    if wcs is not None:
        scale = galsim.PixelScale(pixel_scale)
        psf = wcs.toWorld(scale.toImage(psf), image_pos=ccd_pos)

    return psf


def getBrightPSF(
        ccd_pos=None, wcs=None,
        wavelength=None, gsparams=None,
        logger=None
):
    """Get a fake achromatic optical PSF for really bright objects in Euclidlike simulations.
    Using getPSF() for really bright objects can cause some reflections in the spider pattern
    and some boxiness at the outskirts of the PSF.

    Args:
    ccd_pos:  Single galsim.PositionD indicating the position within the ccd
        for which the PSF should be created. If None, the exact center of the
        ccd is chosen. [default: None]
    wcs:  The WCS to use to project the PSF into world coordinates. [default:
        galsim.PixelScale(euclidlike.pixel_scale)]
    wavelength (float):  An option to get an achromatic PSF for a single
        wavelength, for users who do not care about chromaticity of the PSF. If
        None, then the fully chromatic PSF is returned as an
        InterpolatedChromaticObject. Alternatively the user should supply
        either (a) a wavelength in nanometers, and they will get an
        InterpolatedImage object for that wavelength, or (b) a bandpass object,
        in which case they will get an InterpolatedImage objects defined at the
        effective wavelength of that bandpass. [default: None]
    gsparams:  An optional GSParams argument.  See the docstring for GSParams
        for details. [default: None]

    Returns:
        A single PSF object (either an InterpolatedChromaticObject or an
        InterpolatedImage depending on the inputs).

    """
    if wavelength is None:
        euc_bp = getBandpasses()['VIS']
        euc_bp.red_limit = 910
        euc_bp.blue_limit = 540
        wavelength = euc_bp.effective_wavelength

    # ccd_pos: if None, then all should just be center of the ccd.
    if ccd_pos is None:
        ccd_pos = galsim.PositionD(x=n_pix_col/2, y=n_pix_row/2)

    if not isinstance(wavelength, (float, type(None))):
        raise TypeError(
            "wavelength should either be a float, or None."
        )

    # Now get psf model
    psf = _get_single_bright_psf_obj(wavelength, gsparams)
    # Apply WCS.
    # The current version is in arcsec units, but oriented parallel to the
    # image coordinates. So to apply the right WCS, project to pixels using the
    # Euclid mean pixel_scale, then project back to world coordinates with the
    # provided wcs.
    if wcs is not None:
        scale = galsim.PixelScale(pixel_scale)
        psf = wcs.toWorld(scale.toImage(psf), image_pos=ccd_pos)

    return psf


def _get_single_psf_obj(ccd, bandpass, ccd_pos, wavelength, gsparams = None):
    """
    Routine for making a single PSF.  This gets called by `getPSF` after it
    parses all the options that were passed in.  Users will not directly
    interact with this routine.
    """

    wave_list, im_list = get_euclid_wavelength_psf()
    # instantiate psf object from list of images and wavelengths
    psf_obj = galsim.InterpolatedChromaticObject.from_images(im_list, wave_list, gsparams = gsparams)
    if wavelength is not None:
        if isinstance(wavelength, galsim.Bandpass):
            wave = wavelength.effective_wavelength
        else:
            wave = wavelength
        psf_obj = psf_obj.evaluateAtWavelength(wave)

    return psf_obj

def _get_single_bright_psf_obj(wavelength, gsparams = None):
    """
    Routine for making a single PSF.  This gets called by `getPSF` after it
    parses all the options that were passed in.  Users will not directly
    interact with this routine.
    """

    psf_opt = galsim.OpticalPSF(
        diam=diameter,
        obscuration=obscuration,
        lam=wavelength,
        nstruts=6,
        strut_thick=0.1,
        gsparams=gsparams,
        pad_factor=10,
        strut_angle = galsim.Angle(np.pi/9, galsim.radians) ## rotated angle to match euclid psf struts rotation
        )
    psf_bright = psf_opt

    return psf_bright
