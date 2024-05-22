import galsim
import galsim.roman
import euclidlike
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np


def test_sky():
    """Verify numerical accuracy compared to reference code"""
    # select a few positions - Euclid deep field centers
    pos1 = SkyCoord('17:58:55.9 +66:01:03.7', unit=(u.hourangle, u.deg))
    pos2 = SkyCoord('04:04:57.84 -48:25:22.8', unit=(u.hourangle, u.deg))
    pos3 = SkyCoord('03:31:43.6 -28:05:18.6', unit=(u.hourangle, u.deg), frame='icrs')

    # compare our implementation against GalSim's -- they should not agree perfectly because we do
    # exact calculations and GalSim is interpolating.  But <0.1% agreement is expected.
    for pos in [pos1, pos2, pos3]:
        galsim_pos = galsim.CelestialCoord(ra = pos.ra.deg*galsim.degrees,
                                           dec=pos.dec.deg*galsim.degrees)
        sky1 = euclidlike.getSkyLevel(galsim.roman.getBandpasses()['H158'], galsim_pos)
        sky2 = galsim.roman.getSkyLevel(galsim.roman.getBandpasses()['H158'], world_pos=galsim_pos,
                                 exptime=euclidlike.long_exptime)
        # correct galsim.roman outputs to Euclid-like scenario
        sky2 *= euclidlike.collecting_area/galsim.roman.collecting_area*euclidlike.gain
        fracdiff = (sky2-sky1)/sky1
        np.testing.assert_almost_equal(fracdiff, 0, 3)
    return

def test_sky_noarg():
    """Confirm it runs if no input position is given"""
    foo = euclidlike.getSkyLevel(euclidlike.getBandpasses()['VIS'])
    return

def test_sky_units():
    """Confirm it raises an exception if the units are wrong in bandpass"""
    foo = euclidlike.getBandpasses()['VIS']
    foo.wave_type = "fake_value"
    np.testing.assert_raises(ValueError, euclidlike.getSkyLevel, foo)
    return

if __name__ == "__main__":
    testfns = [v for k, v in vars().items() if k[:5] == 'test_' and callable(v)]
    for testfn in testfns:
        testfn()
