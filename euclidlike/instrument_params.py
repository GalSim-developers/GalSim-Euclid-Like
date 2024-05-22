"""
Basic information GalSim needs to produce Euclid-like simulations.

When editing this file, you should also modify the imports in __init__.py
"""
import os
import numpy as np

# Set these quantities
gain = 3.4  # e-/ADU, https://www.euclid-ec.org/public/mission/vis/
pixel_scale = 0.1  # arcsec / pixel
diameter = 1.2  # meters
obscuration = 0.42 # meters, inferred from diameter and collecting area
focal_length = 24.5 # meters
fratio = focal_length/diameter # dimension-less, inferred from diameter and focal_length
collecting_area = 9926 # cm^2, from https://www.euclid-ec.org/science/overview/
long_exptime = 565  # s (for the longer exposures taken when NISP is measuring spectra)
short_exptime = 100  # s (for the shorter exposures taken when NISP is imaging)
read_noise = 4.4  # e-, https://www.euclid-ec.org/public/mission/vis/
n_dithers = 4
n_ccd = 36
n_pix_row = 4132
n_pix_col = 4096

# Physical pixel size
pixel_scale_mm = 0.012  # mm

plate_scale = 8.33  # arcesec / mm

# Mapping DETID to CCDID
det2ccd = {
    0: '1-1',
    1: '1-2',
    2: '1-3',
    20: '1-4',
    19: '1-5',
    18: '1-6',
    3: '2-1',
    4: '2-2',
    5: '2-3',
    23: '2-4',
    22: '2-5',
    21: '2-6',
    6: '3-1',
    7: '3-2',
    8: '3-3',
    26: '3-4',
    25: '3-5',
    24: '3-6',
    9: '4-1',
    10: '4-2',
    11: '4-3',
    29: '4-4',
    28: '4-5',
    27: '4-6',
    12: '5-1',
    13: '5-2',
    14: '5-3',
    32: '5-4',
    31: '5-5',
    30: '5-6',
    15: '6-1',
    16: '6-2',
    17: '6-3',
    35: '6-4',
    34: '6-5',
    33: '6-6',
}

# Maxinum allowed angle from the telecope solar panels to the sun in degrees.
min_sun_angle = 3.
max_sun_angle = 20.

# Items to potentially do later; part of the galsim.roman setup that currently has no correspondence
# here.
#   dark_current
#   nonlinearity_beta
#   reciprocity_alpha - probably not relevant for VIS
#   persistence_coefficients, persistence_fermi_parameters
#   thermal_backgrounds: not important for VIS, only for NISP
#   stray_light_fraction
#   IPC information - probably not relevant for VIS
#
# Currently we are including charge diffusion and jitter as part of the PSF rather than including it
# separately, so there is no parameter for setting them (e.g. a Gaussian sigma) in
# this module.  However, pointing errors / jitter are quantified in the VIS paper.
#
# The impact of the dichroic and of charge transfer inefficiency are not included.  Same for
# saturation (200000 e-). 
