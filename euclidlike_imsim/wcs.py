from astropy.time import Time
import galsim
import euclidlike
from galsim.config import WCSBuilder, RegisterWCSType
from galsim.angle import Angle
from galsim.celestial import CelestialCoord


class EuclidlikeWCS(WCSBuilder):

    def buildWCS(self, config, base, logger):

        req = {
            "CCD": int,
            "ra": Angle,
            "dec": Angle,
            "pa": Angle,
            "mjd": float,
        }
        opt = {
            "saa": Angle,
            "min_sun_angle": float,
            "max_sun_angle": float,
            "force_cvz": bool,
        }

        kwargs, safe = galsim.config.GetAllParams(
            config,
            base,
            req=req,
            opt=opt,
        )
        if "min_sun_angle" in kwargs:
            euclidlike.instrument_params.min_sun_angle = \
                kwargs["min_sun_angle"] * galsim.degrees
            euclidlike.euclidlike_wcs.min_sun_angle = \
                kwargs["min_sun_angle"] * galsim.degrees
        if "max_sun_angle" in kwargs:
            euclidlike.instrument_params.max_sun_angle = \
                kwargs["max_sun_angle"] * galsim.degrees
            euclidlike.euclidlike_wcs.max_sun_angle = \
                kwargs["max_sun_angle"] * galsim.degrees
        if "saa" in kwargs:
            SAA = kwargs["saa"]
        else:
            SAA = None
        pointing = CelestialCoord(ra=kwargs["ra"], dec=kwargs["dec"])
        wcs = euclidlike.getWCS(
            world_pos=pointing,
            PA=kwargs["pa"],
            date=Time(kwargs["mjd"], format="mjd").datetime,
            CCDs=kwargs["CCD"],
            PA_is_FPA=True,
            SAA=SAA,
        )[kwargs["CCD"]]
        return wcs


RegisterWCSType("EuclidlikeWCS", EuclidlikeWCS())
