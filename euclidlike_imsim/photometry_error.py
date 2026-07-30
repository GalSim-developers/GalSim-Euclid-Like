from email.mime import base

import numpy as np
import galsim
from galsim.config import RegisterValueType, InputLoader, RegisterInputType
from galsim.errors import GalSimConfigError
from photerr import EuclidWideErrorModel


def get_euclid_band(obs_kind):
    """
    Get the Euclid band from the obs_kind.

    Parameters
    ----------
    obs_kind : str
        The obs_kind

    Returns
    -------
    band : str
        The Euclid band corresponding to the obs_kind. In the format PhotErr expects.
    """
    if obs_kind == "VIS_LONG" or obs_kind == "VIS_SHORT":
        band = "VIS"
    elif obs_kind == "NISP_Y":
        band = "Y"
    elif obs_kind == "NISP_J":
        band = "J"
    elif obs_kind == "NISP_H":
        band = "H"
    else:
        raise ValueError(f"Invalid obs_kind: {obs_kind}")
    return band


def get_effective_g_and_hlr(
    hlr_bulge,
    hlr_disk,
    g1_bulge,
    g2_bulge,
    g1_disk,
    g2_disk,
    BTT,
):
    """
    Calculate the effective ellipticity for a bulge+disk galaxy.

    Parameters
    ----------
    hlr_bulge : float
        Half-light radius of the bulge component.
    hlr_disk : float
        Half-light radius of the disk component.
    g1_bulge : float
        First component of the ellipticity of the bulge component.
    g2_bulge : float
        Second component of the ellipticity of the bulge component.
    g1_disk : float
        First component of the ellipticity of the disk component.
    g2_disk : float
        Second component of the ellipticity of the disk component.
    BTT : float
        Bulge-to-total ratio.

    Returns
    -------
    g1_eff : float
        First component of the effective ellipticity.
    g2_eff : float
        Second component of the effective ellipticity.
    r_eff : float
        Effective half-light radius.
    """
    g_bulge = g1_bulge + 1j * g2_bulge
    g_disk = g1_disk + 1j * g2_disk

    w_bulge = BTT * hlr_bulge**2
    w_disk = (1.0 - BTT) * hlr_disk**2

    geff = (w_bulge * g_bulge + w_disk * g_disk) / (w_bulge + w_disk)
    r_eff = np.sqrt(BTT * hlr_bulge**2 + (1.0 - BTT) * hlr_disk**2)

    return geff.real, geff.imag, r_eff


def get_minor_major_axis(hlr, g1, g2):
    """
    Calculate the minor and major axes of a galaxy.

    Parameters
    ----------
    hlr : float
        Half-light radius (in any system of units).
    g1 : float
        First component of the ellipticity.
    g2 : float
        Second component of the ellipticity.

    Returns
    -------
    a : float
        Major axis (same units as input hlr).
    b : float
        Minor axis.
    """

    g = np.hypot(g1, g2)

    q = (1 - g) / (1 + g)

    a = hlr / np.sqrt(q)
    b = hlr * np.sqrt(q)

    return a, b


def prepare_inputs(
    minor_, major_, g1_, g2_, hlr_, g1_disk_, g2_disk_, hlr_disk_, g1_bulge_, g2_bulge_, hlr_bulge_, BTT_
):
    """
    Prepare the inputs for the photometric error calculation.

    Parameters
    ----------
    minor_ : float or None
        Minor axis.
    major_ : float or None
        Major axis.
    g1_ : float or None
        First component of the ellipticity.
    g2_ : float or None
        Second component of the ellipticity.
    hlr_ : float or None
        Half-light radius.
    g1_disk_ : float or None
        First component of the ellipticity of the disk component.
    g2_disk_ : float or None
        Second component of the ellipticity of the disk component.
    hlr_disk_ : float or None
        Half-light radius of the disk component.
    g1_bulge_ : float or None
        First component of the ellipticity of the bulge component.
    g2_bulge_ : float or None
        Second component of the ellipticity of the bulge component.
    hlr_bulge_ : float or None
        Half-light radius of the bulge component.
    BTT_ : float or None
        Bulge-to-total ratio.
    """

    check = [
        (np.isnan(minor_)),
        (np.isnan(major_)),
        (np.isnan(g1_)),
        (np.isnan(g2_)),
        (np.isnan(hlr_)),
        (np.isnan(g1_disk_)),
        (np.isnan(g2_disk_)),
        (np.isnan(hlr_disk_)),
        (np.isnan(g1_bulge_)),
        (np.isnan(g2_bulge_)),
        (np.isnan(hlr_bulge_)),
        (np.isnan(BTT_)),
    ]
    if any(check):
        # This should handle the case where we get this for let's say stars and those quantities do not exist in the catalogue
        return None, None

    # check minor / major
    check = [(minor_ is None), (major_ is None)]
    if not any(check) and all(check):
        raise ValueError("Both minor and major axes must be provided or both must be None.")
    if minor_ is not None:
        minor = np.array([minor_])
        major = np.array([major_])
        return minor, major

    # check g1 / g2 / hlr
    check = [(g1_ is None), (g2_ is None), (hlr_ is None)]
    if not any(check) and all(check):
        raise ValueError("g1, g2, and hlr must all be provided or all must be None.")
    if g1_ is not None:
        g1 = np.array([g1_])
        g2 = np.array([g2_])
        hlr = np.array([hlr_])
        major, minor = get_minor_major_axis(hlr, g1, g2)
        return minor, major

    # check g1_disk / g2_disk / hlr_disk / g1_bulge / g2_bulge / hlr_bulge / BTT
    check = [
        (g1_disk_ is None),
        (g2_disk_ is None),
        (hlr_disk_ is None),
        (g1_bulge_ is None),
        (g2_bulge_ is None),
        (hlr_bulge_ is None),
        (BTT_ is None),
    ]
    if not any(check) and all(check):
        raise ValueError(
            "g1_disk, g2_disk, hlr_disk, g1_bulge, g2_bulge, hlr_bulge, and BTT must all be provided or all must be None."
        )
    if g1_disk_ is not None:
        g1_disk = np.array([g1_disk_])
        g2_disk = np.array([g2_disk_])
        hlr_disk = np.array([hlr_disk_])
        g1_bulge = np.array([g1_bulge_])
        g2_bulge = np.array([g2_bulge_])
        hlr_bulge = np.array([hlr_bulge_])
        BTT = np.array([BTT_])
        g1_eff, g2_eff, hlr_eff = get_effective_g_and_hlr(
            hlr_bulge, hlr_disk, g1_bulge, g2_bulge, g1_disk, g2_disk, BTT
        )
        major, minor = get_minor_major_axis(hlr_eff, g1_eff, g2_eff)
        return minor, major

    return None, None


def PhotErr(config, base, value_type):
    """
    Return the photometric error for a given magnitude and band (obs_kind).
    """

    try:
        photerr_config = galsim.config.GetInputObj("photerr_config", config, base, "PhotErr")
    except GalSimConfigError:
        photerr_config = PhotErrConfig()

    req = {"mag": float, "obs_kind": str}
    opt = {
        "minor": float,
        "major": float,
        "g1": float,
        "g2": float,
        "hlr": float,
        "g1_disk": float,
        "g2_disk": float,
        "hlr_disk": float,
        "g1_bulge": float,
        "g2_bulge": float,
        "hlr_bulge": float,
        "BTT": float,
        "obj_type": str,
    }
    params, safe = galsim.config.GetAllParams(config, base, req=req, opt=opt)
    mag_ = params["mag"]
    minor_ = params.get("minor", None)
    major_ = params.get("major", None)
    g1_ = params.get("g1", None)
    g2_ = params.get("g2", None)
    hlr_ = params.get("hlr", None)
    g1_disk_ = params.get("g1_disk", None)
    g2_disk_ = params.get("g2_disk", None)
    hlr_disk_ = params.get("hlr_disk", None)
    g1_bulge_ = params.get("g1_bulge", None)
    g2_bulge_ = params.get("g2_bulge", None)
    hlr_bulge_ = params.get("hlr_bulge", None)
    BTT_ = params.get("BTT", None)
    obj_type = params.get("obj_type", None)
    obs_kind = params["obs_kind"]

    # Get Euclid band from obs_kind
    band = get_euclid_band(obs_kind)

    mag = np.atleast_2d(mag_)
    minor, major = prepare_inputs(
        minor_,
        major_,
        g1_,
        g2_,
        hlr_,
        g1_disk_,
        g2_disk_,
        hlr_disk_,
        g1_bulge_,
        g2_bulge_,
        hlr_bulge_,
        BTT_,
    )

    euclid_err_model = EuclidWideErrorModel(**photerr_config.get_config(obj_type))

    if minor is None or major is None:
        euclid_err_model = EuclidWideErrorModel(extendedType="point")

    obs_mag, obs_mag_err, obs_flux = euclid_err_model._get_obs_and_errs(
        mag - 2.5 * np.log10(photerr_config.gain[band]),
        major,
        minor,
        [band],
        np.random.default_rng(base["obj_num_seed"]),
    )

    return obs_mag_err[0][0], safe


class PhotErrLoader(InputLoader):
    def getKwargs(self, config, base, logger):
        req = {}
        opt = {
            "m5": dict,
            "theta": dict,
            "extendedType": (str, dict),
            "gain": dict,
        }

        kwargs, safe = galsim.config.GetAllParams(config, base, req=req, opt=opt)

        if isinstance(kwargs.get("gain", None), dict):
            opt = {
                "VIS": float,
                "Y": float,
                "J": float,
                "H": float,
            }
            gains, safe = galsim.config.GetAllParams(kwargs["gain"], base, req={}, opt=opt)
            kwargs["gain"] = gains

        return kwargs, True


class PhotErrConfig:
    def __init__(self, m5=None, theta=None, extendedType=None, gain=None):
        self.m5 = m5
        self.theta = theta
        if isinstance(extendedType, str):
            self.extendedType = {"": extendedType}
        elif isinstance(extendedType, dict):
            self.extendedType = extendedType

        self._set_gain(gain)

    def get_config(self, obj_type=None):
        config = {
            "m5": self.m5,
            "theta": self.theta,
        }
        if obj_type is not None:
            config["extendedType"] = self.extendedType[obj_type]
        else:
            config["extendedType"] = self.extendedType[""]
        return config

    def _set_gain(self, gain):
        self.gain = {
            "VIS": 1.0,
            "Y": 1.0,
            "J": 1.0,
            "H": 1.0,
        }
        if isinstance(gain, dict):
            self.gain.update(gain)


RegisterInputType("photerr_config", PhotErrLoader(PhotErrConfig))
RegisterValueType(
    "PhotErr",
    PhotErr,
    [float, int, None],
)
