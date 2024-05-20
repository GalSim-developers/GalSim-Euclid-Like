from galsim_euclidlike import euclidlike_psf


def test_psf_fake():
    scale = 0.01
    nsample = 10
    npix = 257
    wl_array, psfobjs = euclidlike_psf.get_fake_wavelength_psf(
        scale, nsample, npix,
    )
    assert len(psfobjs) == nsample
    assert psfobjs[0].image.array.shape == (npix, npix)
    assert wl_array.shape == (nsample, )
    return
