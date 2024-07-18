import numpy as np
import galsim
import euclidlike

def test_euclid_bandpass():
    # Obtain the bandpasses with AB_zeropoint set
    bp = euclidlike.getBandpasses(AB_zeropoint=True)

    # Check existence - did we get all the bands?
    for band in euclidlike.vis_bands:
        assert bp[band]
    for band in euclidlike.nisp_bands:
        assert bp[band]
    
    # Check if the zeropoints have been set correctly
    AB_spec = lambda x: (3631e-23)
    AB_sed = galsim.SED(spec=AB_spec, wave_type='nm', flux_type='fnu')
    for filter_name, filter_ in bp.items():
        mag = AB_sed.calculateMagnitude(bandpass=filter_)
        np.testing.assert_almost_equal(mag,0.0,decimal=6,
            err_msg="Zeropoint not set accurately enough for bandpass filter "+filter_name)

    # Check if the zeropoints have not been set
    nozp_bp = euclidlike.getBandpasses(AB_zeropoint=False)
    for key in nozp_bp:
        assert nozp_bp[key].zeropoint is None
    return

if __name__ == "__main__":
    testfns = [v for k, v in vars().items() if k[:5] == 'test_' and callable(v)]
    for testfn in testfns:
        testfn()