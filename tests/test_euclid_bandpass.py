import numpy as np
import galsim
from euclidlike import bandpass

def test_euclid_bandpass():
    # Obtain the bandpasses with AB_zeropoint set
    bp = bandpass.getBandpasses(AB_zeropoint=True)

    # Check if the zeropoints have been set correctly
    AB_spec = lambda x: (3631e-23)
    AB_sed = galsim.SED(spec=AB_spec, wave_type='nm', flux_type='fnu')
    for filter_name, filter_ in bp.items():
        mag = AB_sed.calculateMagnitude(bandpass=filter_)
        np.testing.assert_almost_equal(mag,0.0,decimal=6,
            err_msg="Zeropoint not set accurately enough for bandpass filter "+filter_name)

    # Check if the zeropoints have not been set
    nozp_bp = galsim.roman.getBandpasses(AB_zeropoint=False)
    for key in nozp_bp:
        assert nozp_bp[key].zeropoint is None
    return
