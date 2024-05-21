from euclidlike import get_fake_wavelength_psf, getPSF


def test_psf_fake():
    scale = 0.01
    nsample = 10
    npix = 257
    wl_array, psfobjs = get_fake_wavelength_psf(
        scale, nsample, npix,
    )
    assert len(psfobjs) == nsample
    assert psfobjs[0].image.array.shape == (npix, npix)
    assert wl_array.shape == (nsample, )
    return

def test_get_psf_function():
    fake_scale_default = 0.02
    fake_npix_default = 257
    psfobj = getPSF(ccd=1, bandpass="VIS", wavelength=814)
    assert psfobj.image.array.shape == (fake_npix_default, fake_npix_default)
    assert psfobj.image.scale == fake_scale_default
    return

