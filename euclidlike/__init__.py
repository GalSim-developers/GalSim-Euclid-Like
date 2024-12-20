from .instrument_params import gain, pixel_scale, diameter, obscuration, collecting_area
from .instrument_params import focal_length, fratio
from .instrument_params import vis_bands, nisp_bands
from .instrument_params import long_exptime, short_exptime_vis, nisp_exptime_eff, read_noise, n_dithers
from .instrument_params import n_ccd, n_pix_row, n_pix_col, pixel_scale_mm, det2ccd
from .instrument_params import vis_red_limit, vis_blue_limit
from . import instrument_params
from .bandpass import getBandpasses
from .euclidlike_wcs import getWCS, findCCD, allowedPos, bestPA, convertCenter

from .euclidlike_psf import (
   getPSF, getBrightPSF
)
from .backgrounds import getSkyLevel, getZodiBackground

from ._version import __version__, __version_info__
version = __version__
