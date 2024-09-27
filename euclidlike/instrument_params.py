"""
Basic information GalSim needs to produce Euclid-like simulations.

When editing this file, you should also modify the imports in __init__.py
"""

import galsim

# Set these quantities
gain = 3.4  # e-/ADU, https://www.euclid-ec.org/public/mission/vis/
pixel_scale = 0.1  # arcsec / pixel
diameter = 1.2  # meters
obscuration = 0.42  # meters, inferred from diameter and collecting area
focal_length = 24.5  # meters
fratio = focal_length / diameter  # dimension-less, inferred from diameter and focal_length
collecting_area = 9926  # cm^2, from https://www.euclid-ec.org/science/overview/
long_exptime = 566  # s (for the longer exposures used for VIS images)
short_exptime_nisp = 112  # s (for the shorter NISP imaging exposures)
short_exptime_vis = 95 # s (for the shorter exposures with VIS taken in parallel with NISP imaging)
read_noise = 4.4  # e-, https://www.euclid-ec.org/public/mission/vis/
saturation = 200000 # e-, from https://www.euclid-ec.org/public/mission/vis/#:~:text=VIS%20pixels%20are%20saturated%20at,IE(AB)%3D16.0.
n_dithers = 4
n_ccd = 36
n_ccd_row = 6
n_ccd_col = 6
n_pix_row = 4132
n_pix_col = 4096

# Physical pixel size
pixel_scale_mm = 0.012  # mm

plate_scale = 8.33  # arcsec / mm
# Mapping DETID to CCDID
det2ccd = {
    0: "1-1",
    1: "1-2",
    2: "1-3",
    20: "1-4",
    19: "1-5",
    18: "1-6",
    3: "2-1",
    4: "2-2",
    5: "2-3",
    23: "2-4",
    22: "2-5",
    21: "2-6",
    6: "3-1",
    7: "3-2",
    8: "3-3",
    26: "3-4",
    25: "3-5",
    24: "3-6",
    9: "4-1",
    10: "4-2",
    11: "4-3",
    29: "4-4",
    28: "4-5",
    27: "4-6",
    12: "5-1",
    13: "5-2",
    14: "5-3",
    32: "5-4",
    31: "5-5",
    30: "5-6",
    15: "6-1",
    16: "6-2",
    17: "6-3",
    35: "6-4",
    34: "6-5",
    33: "6-6",
}
ccd2det = {
    "1-1": 0,
    "1-2": 1,
    "1-3": 2,
    "1-4": 20,
    "1-5": 19,
    "1-6": 18,
    "2-1": 3,
    "2-2": 4,
    "2-3": 5,
    "2-4": 23,
    "2-5": 22,
    "2-6": 21,
    "3-1": 6,
    "3-2": 7,
    "3-3": 8,
    "3-4": 26,
    "3-5": 25,
    "3-6": 24,
    "4-1": 9,
    "4-2": 10,
    "4-3": 11,
    "4-4": 29,
    "4-5": 28,
    "4-6": 27,
    "5-1": 12,
    "5-2": 13,
    "5-3": 14,
    "5-4": 32,
    "5-5": 31,
    "5-6": 30,
    "6-1": 15,
    "6-2": 16,
    "6-3": 17,
    "6-4": 35,
    "6-5": 34,
    "6-6": 33,
}

# Maximum allowed angle from the telescope solar panels to the sun in degrees.
min_sun_angle = 3.0 * galsim.degrees
max_sun_angle = 20.0 * galsim.degrees

# Define variables that distinguish between VIS and NISP bandpasses
vis_bands = ['VIS']
nisp_bands = ['NISP_Y', 'NISP_J', 'NISP_H']

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
