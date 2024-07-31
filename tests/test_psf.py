import galsim
import euclidlike
import numpy as np

from euclidlike import (
     get_euclid_wavelength_psf, getPSF, getBrightPSF
)



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

def test_get_bright_psf_function():
    euc_bp = euclidlike.getBandpasses()['VIS']
    euc_bp.red_limit = 910
    euc_bp.blue_limit = 540
    # check optical PSF wavelength at bandpass effective wavelength
    bright_psf = getBrightPSF(1, "VIS", wavelength = 800.)
    np.testing.assert_allclose(
        bright_psf._lam, 800., atol=0,
        err_msg='getBrightPSF() fails to properly initialize OpticalPSF with user-input wavelength')
    
    
    # check size of chromatic bright PSF is reasonable in comparison with euclid-like PSF
    bright_psf = getBrightPSF(ccd=1, bandpass="VIS")
    normal_psf = getPSF(ccd=1, bandpass="VIS")
    scale = euclidlike.pixel_scale
    sed = galsim.SED('vega.txt', 'nm', 'flambda')
    sed = galsim.SED(galsim.LookupTable([100, 2600], [1,1], interpolant='linear'),
                                  wave_type='nm', flux_type='fphotons')
    star = galsim.Convolve(normal_psf, galsim.DeltaFunction(flux = 1e2)*sed)
    star_bright = galsim.Convolve(bright_psf, galsim.DeltaFunction(flux = 1e2)*sed)
    img =  star.drawImage(euc_bp, nx=160, ny=160, scale= scale, method = 'auto')
    img_bright =  star_bright.drawImage(euc_bp, nx=160, ny=160, scale= scale, method = 'auto')
    img_mom = galsim.hsm.FindAdaptiveMom(img)
    img_bright_mom = galsim.hsm.FindAdaptiveMom(img_bright)
    np.testing.assert_allclose(img_bright_mom.moments_sigma, img_mom.moments_sigma, rtol = 0.01,
            err_msg = 'Bright chromatic PSF has size at least 1% larger than getPSF()')
    
    # check size of achromatic bright PSF is reasonable in comparison with euclid-like PSF
    bright_psf = getBrightPSF(ccd=1, bandpass="VIS", wavelength = euc_bp.effective_wavelength)
    star_bright = galsim.Convolve(bright_psf, galsim.DeltaFunction(flux = 1e2)*sed)
    img_bright =  star_bright.drawImage(euc_bp, nx=160, ny=160, scale= scale, method = 'auto')
    img_bright_mom = galsim.hsm.FindAdaptiveMom(img_bright)
    np.testing.assert_allclose(img_bright_mom.moments_sigma, img_mom.moments_sigma, rtol = 0.01,
            err_msg = 'Bright ahcromatic PSF, at effective wavelength, has size at least 1% larger than getPSF()')
    
    # Check nwaves implementation works correctly
    n_waves = 3
    psf_int = getBrightPSF(ccd=1, bandpass="VIS", n_waves=n_waves)
    obj_int = psf_int.evaluateAtWavelength(euc_bp.effective_wavelength  )
    im_int = obj_int.drawImage(scale=scale)
    # Check that evaluation at a single wavelength is consistent with previous results.
    psf_achrom = getBrightPSF(ccd=1, bandpass="VIS", wavelength = euc_bp.effective_wavelength)
    im_achrom = psf_achrom.drawImage(image= im_int.copy(), scale=scale)
    # From roman.getPSF() tests: These images should agree well, but not perfectly.
    #Check for agreement at the level of 1e-3
    diff_im = (im_int.array-im_achrom.array)
    np.testing.assert_array_almost_equal(
        diff_im, np.zeros_like(diff_im), decimal=3,
        err_msg='PSF at a given wavelength and interpolated chromatic one evaluated at that '
        'wavelength disagree.')

if __name__ == "__main__":
    testfns = [v for k, v in vars().items() if k[:5] == 'test_' and callable(v)]
    for testfn in testfns:
        testfn()
