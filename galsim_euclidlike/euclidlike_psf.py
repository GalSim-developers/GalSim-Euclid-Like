import numpy as np
from galsim import roman


def get_fake_wavelength_psf(
    scale, nsample, ngrid, wavelength_min=550, wavelength_max=950,
):
    """ This function get the oversampled PSF image as a function of wavelength
    as an input the Euclid-like simulation

    Args:
    scale (float):  the scale of the PSF image (needs oversampling)
    nsample (int):  number of samples in wavelength
    ngrid (int):  number of grids (number of grids) of the PSF image

    Returns:
    wl_array (ndarray):  wavelength array
    psf_array (ndarray): PSF image array for different wave lengths

    """
    psf_array = np.zeros((nsample, ngrid, ngrid))
    wl_array = np.linspace(wavelength_min, wavelength_max, nsample)
    for i, wl in enumerate(wl_array):
        psf_array[i, :, :] = roman.getPSF(
            8, "W146", wavelength=wl,
        ).drawImage(scale=scale, nx=ngrid, ny=ngrid).array
    return wl_array, psf_array
