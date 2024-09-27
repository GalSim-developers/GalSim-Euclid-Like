import os
import galsim
import numpy as np
from galsim import roman
import astropy.io.fits as pyfits
from importlib.resources import files
from galsim.utilities import LRU_Cache
from . import n_ccd, n_pix_col, n_pix_row, pixel_scale, diameter, obscuration, det2ccd, collecting_area
from .bandpass import getBandpasses
from galsim.config import LoggerWrapper
import pathlib
import warnings

"""
@file euclidlike_psf.py

Part of the Euclid-like simulation module. This file includes routines needed
to define a Euclid-like PSF.  All PSF creation routines are quite approximate, and their limitations
are specified in the routine's docstring.
"""

meta_dir = files("euclidlike.data")
test_dir = files("euclidlike").joinpath('../tests')
# hard coded effective wavelength. Calculate from bandpass object using `bandpass.effective_wavelength`
effective_wave = {'VIS': 718.0867202226793 }
#read wavelengths sampled for PSF
wave_file = files("euclidlike.data").joinpath('psf_wavelengths.dat')
wave_list = np.genfromtxt(wave_file)

# if PSF data is not downloaded, the code defaults to use the PSF for a single quadrant of CCD = 7, stored in tests directory.
default_ccd=7

def _make_psf_list(psf_file):
    image_array = pyfits.getdata(psf_file)
    scale = pixel_scale/3  # images are oversampled by a factor of 3
    im_list = []
    nsample = len(image_array )
    # Provided PSF is normalized such that it includes the effects of obscuration; we want to remove it as obscuration is internally handled by Galsim by using the effective collecting area to compute the measured flux.
    # Diameter is in meters, collecting area in cm^2, need to convert diameter to cm
    m2cm_conv = 1e2
    obscuration_norm = collecting_area / ((m2cm_conv*diameter/2)**2*np.pi)
    for i in range(nsample):
        psf_arr = image_array[i]/obscuration_norm
        im_list.append(
            galsim.Image(psf_arr, scale=scale)
        )
    return im_list

def __get_quadrant_psf(ccd, bandpass, psf_dir):
    ccd_ID = det2ccd[ccd]
    col, row = int(ccd_ID[0]), int(ccd_ID[2])
    # get ccd quadrant IDs
    ul = str(col*2 - 1) + '_' + str(row*2 -1)
    ll = str(col*2 - 1) + '_' + str(row*2)
    ur = str(col*2) + '_' + str(row*2 - 1)
    lr = str(col*2 ) + '_' + str(row*2 )
    quadrants = [ll, ul, lr, ur]
    tags = ["ll", "ul", "lr", "ur"] # ll = lower left, ul = upper left, lr = lower right, ur = upper right
    tag_idx = []
    psf_images = {}
    for tag, CCD_quad in tuple(zip(tags,quadrants)):
        psf_file =psf_dir.joinpath("monopsfs_"+ CCD_quad + ".fits.gz")
        psf_images[tag] = _make_psf_list(psf_file)
    return psf_images

    
_get_quadrant_psf = LRU_Cache(__get_quadrant_psf)

def getPSF(
        ccd, bandpass,
        ccd_pos=None, wcs=None,
        wavelength=None, gsparams=None,
        logger=None, psf_dir = None
):
    """Get a single PSF for a Euclid-like simulation.

    These PSFs are based on precomputed, oversampled (by 3x) PSF images on a grid in wavelength and
    focal plane position.  These images were provided by Lance Miller and Chris Duncan.  They are
    meant to provide an approximately correct description of the PSF size and shape without the
    level of accuracy needed for precision tests of weak lensing shear inference.  The precomputed
    images are combined in a way that accounts for the bandpass and (if provided) SED.  The effects
    of detector linear charge diffusion, a simplified model for guiding errors, and wavefront error
    due to measured surface figure errors are incorporated in the precomputed images.  The telescope
    focus is set to a value from September 2023, and is therefore within a realistic range but does
    not reflect variation over time for any particular observation.

    Several important effects were omitted in these images, including the following:
    * the effect of the dichroic;
    * polarization aberrations; 
    * non-linear detector effects such as charge-transfer inefficiency, the brighter-fatter effect,
      and readout nonlinearity;
    * the effects of decontamination by ice.

    Functionally the lack of non-linear detector effects means the images produced should be thought
    of as post-instrument signature removal (assuming the ISR process is carried out perfectly).

    For applications that require very high accuracy in the modeling of the PSF, with very limited
    aliasing, you may want to lower the `folding_threshold` in the gsparams.  Otherwise very bright
    stars will show some reflections in the spider pattern and possibly some boxiness at the
    outskirts of the PSF due to the size of the precomputed images.  Using ``gsparams =
    GSParams(folding_threshold=1.e-4)`` generally provides good results.
    
    Note that before using, the oversampled PSF images used to create the
    PSF model need to be downloaded. This can be done using the terminal
    command `euclidlike_download_psf`. The images are sampled at the 4 quadrant
    centers of each CCD and at 17 discrete wavelengths.
    The `ccd` argument refers to the detector ID (integer between 0-35),
    not the focal plane position (in format column-row). The sampled
    PSF images are stored using the focal plane position format. Therefore,
    we convert the CCD detector ID to the appropiate focal plane position
    internally.

    In addition, the provided PSF images are normalized for obscuration, 
    vignetting and baffle effects. However, GalSim internally handles the 
    obscuration, so we remove this part of the normalization by dividing the
    PSF images by collecting_area / ((diameter/2)**2*np.pi).  As a result, the sum of the pixel values in the renormalized PSF images is very close to 1.
    

    Args:
    ccd (int):  Single value specifying the CCD for which the PSF
        should be loaded.
    bandpass (str): Single string specifying the bandpass to use when
        defining the pupil plane configuration and/or interpolation of
        chromatic PSFs.
    ccd_pos:  Single galsim.PositionD indicating the position within the CCD
        for which the PSF should be created. If None, the exact center of the
        CCD is chosen. [default: None]
    wcs:  The WCS to use to project the PSF into world coordinates. [default:
        galsim.PixelScale(euclid_like.roman.pixel_scale)]
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
    psf_dir (str): Directory where sampled PSF images can be accessed. If not
        given, look in ./data directory. [default: None] 

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
    
    if psf_dir is None:
        psf_dir = meta_dir.joinpath("monopsfs_euclidlike")
    elif isinstance(psf_dir, str):
        psf_dir = pathlib.Path(psf_dir)
    else:
        raise ValueError("psf_dir must be a string.")
        
    #check if psf directory exists 
    if not psf_dir.is_dir():
        warnings.warn(
            "Unable to use PSF images for full field of view because directory %s does not exist. "
            "Defaulting to use PSF from single quadrant in CCD = %d . All PSF images can be downloaded " 
            "by running the command `euclidlike_download_psf` in the terminal." % (psf_dir, default_ccd)
        )
        psf_dir = test_dir.joinpath("psfs")
        ccd = default_ccd
    logger = LoggerWrapper(logger)
    logger.debug('Loading PSF images from: ' + str(psf_dir))
    # Now get psf model
    psf = _get_single_psf_obj(ccd, bandpass, ccd_pos, wavelength, psf_dir, gsparams, logger)
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
       ccd,
    bandpass,
    ccd_pos=None,
    pupil_bin=4,
    wcs=None,
    n_waves=None,
    wavelength=None,
    gsparams=None,
    logger=None,
):
    """Get a fake optical PSF for very bright objects in Euclid-like simulations.
    Depending on the inputs, this routine returns a chromatic or achromatic PSF using the
    Euclid telescope diameter and Euclid-like aperture.
    
    Args:
    ccd (int):  Single value specifying the CCD for which the PSF should be
        loaded.
    bandpass (str): Single string specifying the bandpass to use when
        defining the pupil plane configuration and/or interpolation of
        chromatic PSFs.
    ccd_pos:  Single galsim.PositionD indicating the position within the CCD
        for which the PSF should be created. If None, the exact center of the
        CCD is chosen. [default: None]
    wcs:  The WCS to use to project the PSF into world coordinates. [default:
        galsim.PixelScale(euclidlike.pixel_scale)]
    n_waves (int): Number of wavelengths to use for setting up interpolation of the
        chromatic PSF objects, which can lead to much faster image
        rendering.  If None, then no interpolation is used. Note that
        users who want to interpolate can always set up the interpolation
        later on even if they do not do so when calling `getPSF`.
        [default: None]
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

    # ccd_pos: if None, then all should just be center of the ccd.
    if ccd_pos is None:
        ccd_pos = galsim.PositionD(x=n_pix_col/2, y=n_pix_row/2)

    if not isinstance(wavelength, (galsim.Bandpass, float, type(None))):
        raise TypeError(
            "wavelength should either be a Bandpass, float, or None."
        )

    # Now call _get_single_PSF().
    psf = _get_single_bright_psf_obj(
        ccd,
        bandpass,
        ccd_pos,
        pupil_bin,
        n_waves,
        wavelength,
        gsparams
    )

    # Apply WCS.
    # The current version is in arcsec units, but oriented parallel to the image coordinates.
    # So to apply the right WCS, project to pixels using the Roman mean pixel_scale, then
    # project back to world coordinates with the provided wcs.
    if wcs is not None:
        scale = galsim.PixelScale(pixel_scale)
        psf = wcs.toWorld(scale.toImage(psf), image_pos=ccd_pos)

    return psf

def _get_single_psf_obj(ccd, bandpass, ccd_pos, wavelength, psf_dir, gsparams, logger):
    """
    Routine for making a single PSF.  This gets called by `getPSF` after it
    parses all the options that were passed in.  Users will not directly
    interact with this routine.
    """

    psf_ims = _get_quadrant_psf(ccd,bandpass, psf_dir)
    # key: ll = lower left, lu = top left, ul = lower right, uu = top right
    # if position lies within the dividing line between quadrants, default
    # is to pick the quadrant center below and/or to the left of the boundary
    quad_row = 'l'
    quad_col = 'l'
    if ccd_pos.y > n_pix_row/2:
        quad_row = 'r'
    if ccd_pos.x > n_pix_col/2:
        quad_col = 'u'
    quad_pos = quad_col + quad_row 
    logger.debug('CCD position in quadrant ' + quad_pos)
    # instantiate psf object from list of images and wavelengths
    psf_obj = galsim.InterpolatedChromaticObject.from_images(psf_ims[quad_pos], wave_list, gsparams = gsparams)
    if wavelength is not None:
        if isinstance(wavelength, galsim.Bandpass):
            wave = wavelength.effective_wavelength
        else:
            wave = wavelength
        psf_obj = psf_obj.evaluateAtWavelength(wave)

    return psf_obj




def _get_single_bright_psf_obj(
    ccd,
    bandpass,
    ccd_pos,
    pupil_bin,
    n_waves,
    wavelength,
    gsparams
):
    """Routine for making a single Optical PSF, designed to handle very bright objects.   
       This gets called by `getBrightPSF` after it parses all the options that were passed in. 
       Users will not directly interact with this routine.
    """

    if wavelength is None:
        wave = effective_wave[bandpass]
    elif isinstance(wavelength, galsim.Bandpass):
        wave = wavelength = wavelength.effective_wavelength
    else:
        wave = wavelength

    # All parameters relevant to the aperture.  We may be able to use a cached version.
    aper = _make_aperture(pupil_bin, wave, gsparams)

    # Now set up the PSF, including the option to interpolate over waves
    if wavelength is None:
        PSF = galsim.ChromaticOpticalPSF(
            lam=wave,
            diam=diameter,
            aper=aper,
            gsparams=gsparams
        )
        if n_waves is not None:
            # To decide the range of wavelengths to use, check the bandpass.
            bp_dict = getBandpasses()
            bp = bp_dict[bandpass]
            PSF = PSF.interpolate(waves=np.linspace(bp.blue_limit, bp.red_limit, n_waves),
                                  oversample_fac=1.5)
    else:
        PSF = galsim.OpticalPSF(lam=wave, diam=diameter,
                         aper=aper, gsparams=gsparams)

    return PSF


def _make_aperture(pupil_bin, wave, gsparams):
    # Load the pupil plane image.
    pupil_plane_im_path = files("euclidlike.data").joinpath("euclid_pupil_plane.fits.gz")
    pupil_plane_im = galsim.fits.read(str(pupil_plane_im_path), read_header=True)

    pupil_plane_im = pupil_plane_im.bin(pupil_bin, pupil_bin)

    aper = galsim.Aperture(
        lam=wave, diam=diameter,
        obscuration=obscuration,
        pupil_plane_im=pupil_plane_im,
        gsparams=gsparams
    )
    return aper
