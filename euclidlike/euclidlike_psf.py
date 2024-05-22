import galsim
import numpy as np
from galsim import roman


def get_fake_wavelength_psf(
    scale, nsample, npix, wavelength_min=550, wavelength_max=950,
):
    """ This function get the oversampled PSF image as a function of wavelength
    as an input to the Euclid-like simulations.

    Args:
    scale (float):  the scale of the PSF image (needs oversampling)
    nsample (int):  number of samples in wavelength
    npix (int):  number of pixels of the PSF image, the area is npix x npix
    wavelength_min (float):  the minimum wave number to sample [units: nm]
    wavelength_max (float):  the maximum wave number to sample [units: nm]

    Returns:
    wl_array (ndarray):  wavelength array
    psfobjs (ndarray): PSF interpolated image for different wavelengths
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
