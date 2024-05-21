import galsim
import numpy as np
from galsim import roman

from . import n_ccd, n_pix_col, n_pix_row, pixel_scale

"""
@file euclidlike_psf.py

Part of the Euclid-like simulation module. This file includes routines needed
to define a Euclid-like PSF.
"""

fake_scale_default = 0.02
fake_npix_default = 257


def get_fake_wavelength_psf(
    scale, nsample, npix, wavelength_min=550, wavelength_max=950,
):
    """ This function get the oversampled PSF image as a function of wavelength
    as an input the Euclid-like simulation

    Args:
    scale (float):  the scale of the PSF image (needs oversampling)
    nsample (int):  number of samples in wavelength
    npix (int):  number of pixels of the PSF image, the area is npix x npix
    wavelength_min (float):  the minimum wave number to sample [units: nm]
    wavelength_max (float):  the maximum wave number to sample [units: nm]

    Returns:
    wl_array (ndarray):  wavelength array
    psfobjs (ndarray): PSF interpolated image for different wave lengths
    """
    psfobjs = []
    wl_array = np.linspace(wavelength_min, wavelength_max, nsample)
    for i, wl in enumerate(wl_array):
        psfobjs.append(
            galsim.InterpolatedImage(
                roman.getPSF(
                    8, "W146", wavelength=wl,
                ).drawImage(scale=scale, nx=npix, ny=npix, method="no_pixel")
            )
        )
    return wl_array, psfobjs


def getPSF(
        ccd, bandpass,
        ccd_pos=None, wcs=None,
        wavelength=None, gsparams=None,
        logger=None
):
    """Get a single PSF for Roman ST observations.

    The user must provide the ccd and bandpass; the latter is used when setting
    up the pupil plane configuration and when interpolating chromatic
    information, if requested.

    The PSF that is returned by default will be oriented with respect to the
    ccd coordinates, not world coordinates as is typical in GalSim.  The pupil
    plane has a fixed orientation with respect to the focal plane, so the PSF
    rotates with the telescope.  To obtain a PSF in world coordinates, which
    can be convolved with galaxies (that are normally described in world
    coordinates), you may pass in a ``wcs`` parameter to this function.  This
    will project the PSF into world coordinates according to that WCS before
    returning it.  Otherwise, the return value is equivalent to using
    ``wcs=galim.PixelScale(galsim.roman.pixel_scale)``.

    The calculation takes advantage of the fact that the diffraction limit and
    aberrations have a simple, understood wavelength-dependence.  (The Roman
    abberation data for Cycle 9 does in fact provide aberrations as a function
    of wavelength, but the deviation from the expected chromatic dependence is
    sub-percent so we neglect it here.)  For reference, the script used to
    parse the Zernikes given on the webpage and create the files in the GalSim
    repository can be found in ``devel/external/parse_roman_zernikes_1217.py``.
    The resulting chromatic object can be used to draw into any of the Roman
    bandpasses, though the pupil plane configuration will only be correct for
    those bands in the same range (i.e., long- or short-wavelength bands).

    For applications that require very high accuracy in the modeling of the
    PSF, with very limited aliasing, you may want to lower the
    folding_threshold in the gsparams.  Otherwise very bright stars will show
    some reflections in the spider pattern and possibly some boxiness at the
    outskirts of the PSF.  Using ``gsparams =
    GSParams(folding_threshold=2.e-3)`` generally provides good results even
    for very bright (e.g. mag=10) stars.

    Jitter and charge diffusion are, by default, not included.

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
        galsim.PixelScale(galsim.roman.pixel_scale)]
    wavelength (float):  An option to get an achromatic PSF for a single
        wavelength, for users who do not care about chromaticity of the PSF.
        If None, then the fully chromatic PSF is returned.  Alternatively the
        user should supply either (a) a wavelength in nanometers, and they will
        get achromatic OpticalPSF objects for that wavelength, or (b) a
        bandpass object, in which case they will get achromatic OpticalPSF
        objects defined at the effective wavelength of that bandpass.
        [default: None]
    gsparams:  An optional GSParams argument.  See the docstring for GSParams
        for details. [default: None]

    Returns:
        A single PSF object (either a ChromaticOpticalPSF or an OpticalPSF
        depending on the inputs).

    """

    if ccd <= 0 or ccd > n_ccd:
        raise galsim.GalSimRangeError("Invalid ccd.", ccd, 1, n_ccd)

    assert bandpass == "VIS", "Only VIS band is supported"

    # ccd_pos: if None, then all should just be center of the ccd.
    if ccd_pos is None:
        ccd_pos = galsim.PositionD(x=n_pix_col/2, y=n_pix_row/2)

    # Now get psf model
    psf = _get_single_psf_obj(ccd, bandpass, ccd_pos, wavelength, gsparams)
    # Apply WCS.
    # The current version is in arcsec units, but oriented parallel to the
    # image coordinates. So to apply the right WCS, project to pixels using the
    # Roman mean pixel_scale, then project back to world coordinates with the
    # provided wcs.
    if wcs is not None:
        scale = galsim.PixelScale(pixel_scale)
        psf = wcs.toWorld(scale.toImage(psf), image_pos=ccd_pos)

    return psf


def _get_single_psf_obj(ccd, bandpass, ccd_pos, wavelength, gsparams):
    """
    Routine for making a single PSF.  This gets called by `getPSF` after it
    parses all the options that were passed in.  Users will not directly
    interact with this routine.
    """

    # Now set up the PSF, including the option to interpolate over waves
    if wavelength is None:
        psf_obj = None
        pass
    else:
        wl_array, psfobjs = get_fake_wavelength_psf(
            scale=fake_scale_default, nsample=17, npix=fake_npix_default,
        )
        index = np.searchsorted(wl_array, wavelength)
        psf_obj = psfobjs[index]
    return psf_obj
