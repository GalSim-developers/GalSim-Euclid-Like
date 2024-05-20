from galsim_euclidlike import euclidlike_psf


def test_psf_fake():
    scale = 0.01
    nsample = 10
    ngrid = 257
    wl_array, psf_array = euclidlike_psf.get_fake_wavelength_psf(
        scale, nsample, ngrid,
    )
    assert psf_array.shape == (nsample, ngrid, ngrid)
    assert wl_array.shape == (nsample, )
    return
