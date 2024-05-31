import numpy as np
from unittest import mock
from numpy.testing import assert_raises

import galsim

import euclidlike


def test_euclidlike_wcs():
    """Test the Euclid WCS
    """

    import datetime

    # Date for testing
    date = datetime.datetime(2025, 1, 12)
    # Position of a Euclid deep field for testing
    world_pos = galsim.CelestialCoord(
        ra=269.73*galsim.degrees,
        dec=66.02*galsim.degrees
    )

    # Check whether we're allowed to look at certain positions on certain dates.
    # Let's choose RA=90 degrees, dec=10 degrees.
    # We know that it's best to look about 90 degrees from the Sun.  So on the vernal and autumnal
    # equinox, this should be a great place to look, but not midway in between.
    pos = galsim.CelestialCoord(90.*galsim.degrees, 10.*galsim.degrees)
    import datetime
    assert euclidlike.allowedPos(pos, datetime.date(2025, 3, 20))
    assert euclidlike.allowedPos(pos, datetime.date(2025, 9, 20))
    assert not euclidlike.allowedPos(pos, datetime.date(2025, 6, 20))
    assert euclidlike.bestPA(pos, datetime.date(2025, 6, 20)) is None

    # Finally make sure it does something reasonable for the observatory position angle.
    # When the sun is at (0,0), and we look at (90,0), then +Z points towards the Sun and +Y points
    # North, giving a PA of 0 degrees.
    pos = galsim.CelestialCoord(90.*galsim.degrees, 0.*galsim.degrees)
    test_date = datetime.datetime(2025, 3, 20, 9, 2)
    pa = euclidlike.bestPA(pos, test_date)
    np.testing.assert_almost_equal(pa.rad, 0., decimal=3)
    # Now make it look at the same RA as the sun but quite different declination. It wants +Z
    # pointing North toward Sun, so we'll get a -90 degree angle for the PA.
    pos = galsim.CelestialCoord(0.*galsim.degrees, -89.*galsim.degrees)
    pa = euclidlike.bestPA(pos, test_date)
    np.testing.assert_almost_equal(pa.rad, -np.pi/2, decimal=3)

    sun_pos = galsim.CelestialCoord(0*galsim.degrees, 0*galsim.degrees)
    sun_pa = euclidlike.bestPA(sun_pos, test_date)
    assert sun_pa is None

    with assert_raises(TypeError):
        euclidlike.getWCS(world_pos=galsim.PositionD(300, 400))
    with assert_raises(galsim.GalSimError):
        euclidlike.getWCS(world_pos=sun_pos, date=test_date)
    with assert_raises(TypeError):
        euclidlike.getWCS(world_pos=pos, PA=33.)
    with assert_raises(galsim.GalSimRangeError):
        euclidlike.getWCS(world_pos=pos, CCDs=[-1, 1])
    with assert_raises(galsim.GalSimRangeError):
        euclidlike.getWCS(world_pos=pos, CCDs=[0, 40])

    south_pole = galsim.CelestialCoord(0*galsim.degrees, -90*galsim.degrees)
    wcs = euclidlike.getWCS(world_pos=south_pole, CCDs=0)
    with assert_raises(TypeError):
        euclidlike.findCCD(wcs_dict=None, world_pos=pos)
    with assert_raises(TypeError):
        euclidlike.findCCD(wcs_dict=wcs, world_pos=galsim.PositionD(300, 400))
    with mock.patch('euclidlike.euclidlike_wcs.pv_filename', 'pv_coeffs.dat'):
        with assert_raises(OSError):
            euclidlike.getWCS(world_pos=world_pos, date=date)


if __name__ == "__main__":
    test_euclidlike_wcs()
