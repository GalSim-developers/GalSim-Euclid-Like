import galsim
import euclidlike
import numpy as np


def test_optics():
    """
    Verify that the optics constants in the param file is self-consistent
    """
    
    diameter = euclidlike.instrument_params.diameter
    focal_length = euclidlike.instrument_params.focal_length
    fratio = euclidlike.instrument_params.fratio

    # Make sure the fratio in the params file is consistent with the focal length
    # and diameter of the main mirror. 
    np.testing.assert_almost_equal(diameter/focal_length, 1/fratio, 5)
    
    # Make sure the physical pixel size is consistent with its angular size.
    pixel_size_physical_meter = euclidlike.instrument_params.pixel_scale_mm / 1e3
    pixel_size_rad = pixel_size_physical_meter/focal_length
    pixel_size_arcsec = pixel_size_rad * 180.0 / np.pi *3600
    
    np.testing.assert_almost_equal(pixel_size_arcsec, euclidlike.instrument_params.pixel_scale, 2)
    
if __name__ == "__main__":
    test_optics()