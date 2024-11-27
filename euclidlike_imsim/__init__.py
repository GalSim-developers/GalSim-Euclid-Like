try:
    from lsst.utils.threads import disable_implicit_threading

    disable_implicit_threading()
except:
    pass
from .obseq import *
from .psf import *
from .ccd import *
from .stamp import *
from .wcs import *
from .skycat import *
from .photonOps import *
from .bandpass import *

# from .detector_physics import *
from ._version import __version__, __version_info__
version = __version__
