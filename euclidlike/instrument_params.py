"""
Basic information GalSim needs to produce Euclid-like simulations.

When editing this file, you should also modify the imports in __init__.py
"""
import os
import numpy as np

# Set these quantities
gain = 3.4 # e-/ADU, https://www.euclid-ec.org/public/mission/vis/
pixel_scale = 0.1  # arcsec / pixel
diameter = 1.2  # meters
obscuration = 0.42 # meters, inferred from diameter and collecting area
focal_length = 24.5 # meters
fratio = 20.42 # dimension-less, inferred from diameter and focal_length
collecting_area = 9926 # cm^2, from https://www.euclid-ec.org/science/overview/
long_exptime = 565  # s (for the longer exposures taken when NISP is measuring spectra)
short_exptime = 100 # s (for the shorter exposures taken when NISP is imaging)
read_noise = 4.4 # e-, https://www.euclid-ec.org/public/mission/vis/
n_dithers = 4
n_ccd = 36
n_pix_row = 4132
n_pix_col = 4096

# Physical pixel size
pixel_scale_mm = 0.012 # mm

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
