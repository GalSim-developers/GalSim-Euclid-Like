import galsim
import euclidlike
import numpy as np
import os
from importlib.resources import files
from euclidlike import (
      getPSF, getBrightPSF
)
import astropy.io.fits as pyfits


script_dir = os.path.dirname(__file__)
psf_dir = os.path.join(script_dir, 'psfs') # directory for psfs, monopsfs_6_6.fits.gz is real, the rest are symlinks to this one for testin

# we want to test the PSF at some random focal plane positions. We use these two ccd's for testing purposes
test_ccd=7
test_adj_ccd=8

# set tolerance levels for different tests
flux_atol = 1e-4
size_rtol = 0.05

def test_get_psf_function():
    euc_bp = euclidlike.getBandpasses()['VIS']
    scale = euclidlike.pixel_scale/3
    
    wave_file = files("euclidlike.data").joinpath('psf_wavelengths.dat')
    wave_list = np.genfromtxt(wave_file)
    psf_file = os.path.join(psf_dir, "monopsfs_6_6.fits.gz")
    psfobjs = euclidlike.euclidlike_psf._make_psf_list(psf_file)
    wave = wave_list[9]  # random wavelength from oversampled images
    trueobj = psfobjs[9]
    psfobj = getPSF(ccd=test_ccd, bandpass="VIS", wavelength=wave, psf_dir = psf_dir)
    # check it returns identical oversampled psf
    np.testing.assert_allclose(psfobj.image.array, trueobj.array, atol = 0,
        err_msg = 'getPSF() fails to initialize input images correctly')

    # check sum of PSF image pixels after normalizing by obscuration is within 1% of 1.
    # It won't add to exactly 1 since PSF includes vignetting and baffle effects.
    for i in range(len(psfobjs)):
          np.testing.assert_allclose(np.sum(psfobjs[i].array), 1.0, atol = 0.01,
              err_msg = 'Sum of PSF pixels != 1 by more than 1%') 
      
    # check different ccd works
    psfobj = getPSF(ccd=test_adj_ccd, bandpass="VIS", psf_dir = psf_dir)

    #check 4 quadrants positions on ccd
    col_move = [20, -20,-20, 20 ]
    row_move = [20, -20,20, -20 ]
    for i in range(4): 
        ccd_pos = galsim.PositionD(x=euclidlike.n_pix_col/2 + col_move[i], y=euclidlike.n_pix_row/2 + row_move[i])
        psf_ur = getPSF(ccd=test_ccd, bandpass="VIS",ccd_pos = ccd_pos, psf_dir = psf_dir)

    #check passing bandpass
    psfobj_bp = getPSF(ccd=test_ccd, bandpass="VIS", wavelength=euc_bp, psf_dir = psf_dir)
    psfobj_wl = getPSF(ccd=test_ccd, bandpass="VIS", wavelength=euc_bp.effective_wavelength, psf_dir = psf_dir)
    np.testing.assert_allclose(psfobj_bp.image.array, psfobj_wl.image.array, atol = 0,
        err_msg = 'getPSF() fails to reproduce image if bandpass is the input wavelength')

    #check full PSF with delta function SED at desired wavelength returns identical image
    psfobj = getPSF(ccd=test_ccd, bandpass="VIS", psf_dir = psf_dir)

    # creating narrow SED at desired wavelength
    x = np.linspace(euc_bp.blue_limit, euc_bp.red_limit, 500)
    x[285] =  wave #index 285 roughly corresponds to this wavelength, manually setting
    f = np.full(len(x), 1e-30)
    f[285] = 100
    lk_table = galsim.LookupTable(x = x, f = f)
    deltaf_sed = galsim.SED(lk_table, wave_type = 'nm', flux_type = 'fphotons')
    star = galsim.Gaussian(fwhm=1.e-8) * deltaf_sed

    # drawing objects
    trueobj = getPSF(ccd=test_ccd, bandpass="VIS", wavelength=wave, psf_dir = psf_dir)
    psf_obj = galsim.Convolve(psfobj, star)
    true_obj = galsim.Convolve(trueobj, star)
    np.testing.assert_allclose(psfobj.ims[9].array, trueobj.image.array, atol=0,
        err_msg='getPSF() with wavelength=None fails to initialize input images correctly')
    psf_img = psf_obj.drawImage(euc_bp, scale=scale)
    im_interp = psf_img.copy()
    true_img = true_obj.drawImage(euc_bp, scale=scale)
    # images identicial within 0.01% of total flux
    np.testing.assert_allclose(
        psf_img.array, true_img.array, atol=flux_atol*np.sum(true_img.array),
        err_msg='getPSF() does replicate image with very narrow SED centered at desired wavelength')
    return

def test_get_bright_psf_function():
    euc_bp = euclidlike.getBandpasses()['VIS']
    # check optical PSF wavelength at bandpass effective wavelength
    bright_psf = getBrightPSF(1, "VIS", wavelength = 800.)
    np.testing.assert_allclose(
        bright_psf._lam, 800., atol=0,
        err_msg='getBrightPSF() fails to properly initialize OpticalPSF with user-input wavelength')
    
    
    # check size of chromatic bright PSF is reasonable in comparison with euclid-like PSF
    bright_psf = getBrightPSF(ccd=test_ccd, bandpass="VIS")
    normal_psf = getPSF(ccd=test_ccd, bandpass="VIS", psf_dir = psf_dir)
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
    np.testing.assert_allclose(img_bright_mom.moments_sigma, img_mom.moments_sigma, rtol = size_rtol,
            err_msg = 'Bright chromatic PSF has size at least %f%% larger than getPSF()'%(100* size_rtol))
    
    # check size of achromatic bright PSF is reasonable in comparison with euclid-like PSF
    bright_psf = getBrightPSF(ccd=1, bandpass="VIS", wavelength = euc_bp.effective_wavelength)
    star_bright = galsim.Convolve(bright_psf, galsim.DeltaFunction(flux = 1e2)*sed)
    img_bright =  star_bright.drawImage(euc_bp, nx=160, ny=160, scale= scale, method = 'auto')
    img_bright_mom = galsim.hsm.FindAdaptiveMom(img_bright)
    np.testing.assert_allclose(img_bright_mom.moments_sigma, img_mom.moments_sigma, rtol = size_rtol,
            err_msg = 'Bright ahcromatic PSF, at effective wavelength, has size at least %f%% larger than getPSF()'%(100* size_rtol))
    
    # Check nwaves implementation works correctly
    n_waves = 3
    psf_int = getBrightPSF(ccd=test_ccd, bandpass="VIS", n_waves=n_waves)
    obj_int = psf_int.evaluateAtWavelength(euc_bp.effective_wavelength  )
    im_int = obj_int.drawImage(scale=scale)
    # Check that evaluation at a single wavelength is consistent with previous results.
    psf_achrom = getBrightPSF(ccd=test_ccd, bandpass="VIS", wavelength = euc_bp.effective_wavelength)
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
