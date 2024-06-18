import galsim
import euclidlike
import numpy as np

from euclidlike import (
    get_fake_wavelength_psf, get_euclid_wavelength_psf, getPSF
)


def test_psf_fake():
    scale = 0.01
    nsample = 10
    npix = 257
    wl_array, psfobjs = get_fake_wavelength_psf(
        scale, nsample, npix,
    )
    assert len(psfobjs) == nsample
    assert psfobjs[0].array.shape == (npix, npix)
    assert wl_array.shape == (nsample, )
    return


def test_psf_euclid():
    wl_array, psfobjs = get_euclid_wavelength_psf()
    assert len(psfobjs) == 17
    assert psfobjs[0].array.shape == (480, 480)
    assert wl_array.shape == (17, )
    for i in range(len(psfobjs)):
        assert isinstance(
            psfobjs[i], galsim.Image,
        )  # check objects in list are Galsim Image objects
    return


def test_get_psf_function():
    euc_bp = euclidlike.getBandpasses()['VIS']
    euc_bp.red_limit = 910
    euc_bp.blue_limit = 540
    scale = euclidlike.pixel_scale/3

    wl_array, psfobjs = get_euclid_wavelength_psf()
    wave = wl_array[9]  # random wavelength from oversampled images
    trueobj = psfobjs[9]
    psfobj = getPSF(ccd=1, bandpass="VIS", wavelength=wave)
    # check it returns identical oversampled psf
    np.testing.assert_allclose(psfobj.image.array, trueobj.array, atol = 0,
        err_msg = 'getPSF() fails to initialize input images correctly')

    #check passing bandpass
    psfobj_bp = getPSF(ccd=1, bandpass="VIS", wavelength=euc_bp)
    psfobj_wl = getPSF(ccd=1, bandpass="VIS", wavelength=euc_bp.effective_wavelength)
    np.testing.assert_allclose(psfobj_bp.image.array, psfobj_wl.image.array, atol = 0,
        err_msg = 'getPSF() fails to reproduce image if bandpass is the input wavelength')

    #check full PSF with delta function SED at desired wavelength returns identical image
    psfobj = getPSF(ccd=1, bandpass="VIS")

    # creating narrow SED at desired wavelength
    x = np.linspace(euc_bp.blue_limit, euc_bp.red_limit, 500)
    x[285] =  wave #index 285 roughly corresponds to this wavelength, manually setting
    f = np.full(len(x), 1e-30)
    f[285] = 100
    lk_table = galsim.LookupTable(x = x, f = f)
    deltaf_sed = galsim.SED(lk_table, wave_type = 'nm', flux_type = 'fphotons')
    star = galsim.Gaussian(fwhm=1.e-8) * deltaf_sed

    # drawing objects
    trueobj = getPSF(ccd=1, bandpass="VIS", wavelength=wave)
    psf_obj = galsim.Convolve(psfobj, star)
    true_obj = galsim.Convolve(trueobj, star)
    np.testing.assert_allclose(psfobj.ims[9].array, trueobj.image.array, atol=0,
        err_msg='getPSF() with wavelength=None fails to initialize input images correctly')
    psf_img = psf_obj.drawImage(euc_bp, scale=scale)
    im_interp = psf_img.copy()
    true_img = true_obj.drawImage(euc_bp, scale=scale)
    # images identicial within 0.01% of total flux
    np.testing.assert_allclose(
        psf_img.array, true_img.array, atol=1e-4*np.sum(true_img.array),
        err_msg='getPSF() does replicate image with very narrow SED centered at desired wavelength')
    return

if __name__ == "__main__":
    testfns = [v for k, v in vars().items() if k[:5] == 'test_' and callable(v)]
    for testfn in testfns:
        testfn()
