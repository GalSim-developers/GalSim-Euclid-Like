import numpy as np
import numba as nb

from astropy.time import Time

import galsim
from galsim.errors import GalSimConfigError

import ngmix

import euclidlike

MIN_IMG_SIZE = 51
MAX_IMG_SIZE = 501

nisp_gain = 2  # https://arxiv.org/pdf/2405.13496 Sect 4.3.7
nisp_pixel_scale = 0.3
nisp_dark_current = 0.02  # e.s-1.pix-1 https://arxiv.org/pdf/2405.13493 Sect 4.1.2
nisp_read_noise = 6.2  # e.pix-1 https://arxiv.org/pdf/2405.13496 Sect 4.3.5


def get_good_img_size(gmix, scale):

    _, _, sigma = gmix.get_e1e2sigma()
    size_pix = sigma / scale
    img_size = min(max(size_pix * 5, MIN_IMG_SIZE), MAX_IMG_SIZE)
    img_size = int(np.ceil(img_size))
    return img_size


@nb.njit
def get_flux(gmix, pixels):
    aper_flux = 0
    n_pixels = pixels.shape[0]
    for i in range(n_pixels):
        aper_flux += ngmix.gmix.gmix_nb.gmix_eval_pixel_fast(gmix, pixels[i])
    return aper_flux


@nb.njit
def circular_aper(N, r):
    """
    Draws a circle of radius r in a square image of size N x N.

    Parameters:
    N (int): The size of the square image (N x N).
    r (int): The radius of the circle.

    Returns:
    numpy.ndarray: A 2D numpy array representing the image with the circle.
    """
    # Create an empty N x N array (filled with zeros)
    image = np.zeros((N, N), dtype=np.int16)

    # Define the center of the image
    center = N // 2

    # Iterate over each pixel in the image
    for y in range(N):
        for x in range(N):
            # Calculate the distance from the center
            distance = np.sqrt((x - center)**2 + (y - center)**2)

            # If the distance is less than or equal to the radius, set the pixel to 1
            if distance <= r:
                image[y, x] = 1

    return image


def get_variance(bandpass, world_pos, mjd, exptime):

    filter = bandpass.name

    sky_lvl = euclidlike.getSkyLevel(
        bandpass,
        world_pos=world_pos,
        date=Time(mjd, format="mjd").datetime,
        exptime=exptime,
    )

    if filter == "VIS":
        return (
            sky_lvl * euclidlike.pixel_scale**2
            + euclidlike.read_noise**2
        ) / euclidlike.gain**2
    elif "NISP" in filter:
        return (
            sky_lvl * nisp_pixel_scale**2
            + nisp_read_noise**2
            + (nisp_dark_current * exptime)**2
        ) / nisp_gain**2


def get_flux_err(
    ra,
    dec,
    wcs,
    bandpass,
    mjd,
    exptime,
    aper_type,
    aper_size,
    Nexp,
    gal_pars=None,
    star_flux=None,
    model=None,
):
    if (
        (gal_pars is None and star_flux is None)
        | (gal_pars is not None and star_flux is not None)
    ):
        raise ValueError(
            "One of [gal_pars, star_flux] must be provided, not both."
        )
    if model is None:
        raise ValueError("Model must be provided in ['bd', 'gausse']")

    filter = bandpass.name
    wave = bandpass.effective_wavelength

    lam_over_diam = np.rad2deg(
        wave * 1e-9 / euclidlike.diameter
    ) * 3600
    psf_pars = [
        0, 0, 0, 0, 2 * (0.5 * lam_over_diam)**2, 1
    ]

    if gal_pars is not None:
        gal_gmix = ngmix.gmix.make_gmix_model(gal_pars, model)
        psf_gmix = ngmix.gmix.make_gmix_model(psf_pars, "gauss")
        obj_gmix = gal_gmix.convolve(psf_gmix)
    else:
        point_source_pars = psf_pars
        point_source_pars[-1] = star_flux
        obj_gmix = ngmix.gmix.make_gmix_model(point_source_pars, "gauss")

    world_pos = galsim.CelestialCoord(
        ra=ra * galsim.degrees, dec=dec * galsim.degrees
    )
    if filter == "VIS":
        jacobian = wcs.jacobian(world_pos=world_pos)
    else:
        jacobian = wcs.jacobian()
    pixel_scale = np.sqrt(jacobian.pixelArea())

    img_size = get_good_img_size(obj_gmix, pixel_scale)

    jacob = ngmix.Jacobian(
        row=(img_size - 1) / 2,
        col=(img_size - 1) / 2,
        wcs=jacobian,
    )

    gmix_img = obj_gmix.make_image(
        (img_size, img_size),
        jacobian=jacob,
        fast_exp=True,
    )

    noise_var = get_variance(bandpass, world_pos, mjd, exptime)
    if aper_type == "circular":
        if aper_size < 0:
            aperture_mask = np.ones((img_size, img_size))
        else:
            aperture_mask = circular_aper(img_size, aper_size / pixel_scale)
    else:
        raise NotImplementedError
    weight_img = (
        aperture_mask
        * np.ones_like(gmix_img) / ((noise_var / Nexp + gmix_img))
    )
    pixels = ngmix.pixels.make_pixels(
        gmix_img,
        weight_img,
        jacob,
    )

    flux = get_flux(
        obj_gmix.get_data(),
        pixels,
    )

    snr = np.sqrt(
        ngmix.gmix.gmix_nb.get_model_s2n_sum(
            obj_gmix.get_data(),
            pixels,
        )
    )

    flux_err = flux / snr

    return flux, flux_err, snr
