from astropy.time import Time
import pandas as pd
import galsim
import galsim.config
from galsim.angle import Angle
from galsim.config import InputLoader, RegisterValueType, RegisterInputType
from galsim.errors import GalSimConfigValueError


OBS_KIND = [
    "VIS_LONG",
    "VIS_SHORT",
    "NISP_J",
    "NISP_H",
    "NISP_Y",
]

class ObSeqDataLoader(object):
    """Read the exposure information from the observation sequence."""

    _req_params = {"file_name": str, "visit": int, "obs_kind": str, "CCD": int}

    def __init__(self, file_name, visit, obs_kind, CCD, logger=None):
        self.logger = galsim.config.LoggerWrapper(logger)
        self.file_name = file_name
        self.visit = visit
        self.ccd = CCD

        if obs_kind not in OBS_KIND:
            raise GalSimConfigValueError(
                "Invalid obs_kind.", obs_kind, OBS_KIND
            )
        self.obs_kind = obs_kind

        # try:
        self.read_obseq()
        # except:
        #     # Read visit info from the config file.
        #     self.logger.warning('Reading visit info from config file.')

    def read_obseq(self):
        """Read visit info from the obseq file."""
        if self.file_name is None:
            raise ValueError(
                "No obseq filename provided, trying to build from config information."
            )
        if self.visit is None:
            raise ValueError(
                "The visit must be set when reading visit info from an obseq file."
            )

        self.logger.warning(
            "Reading info from obseq file %s for visit %s", self.file_name, self.visit
        )

        ob = pd.read_pickle(self.file_name).loc[self.visit]

        self.ob = {}
        for obs_kind in OBS_KIND:
            _ob = {}
            _ob["visit"] = self.visit
            _ob["ccd"] = self.ccd
            _ob["ra"] = ob.loc[obs_kind]["ra"] * galsim.degrees
            _ob["dec"] = ob.loc[obs_kind]["dec"] * galsim.degrees
            _ob["pa"] = ob.loc[obs_kind]["pa"] * galsim.degrees
            _ob["saa"] = ob.loc[obs_kind]["saa"] * galsim.degrees
            _ob["date"] = Time(ob.loc[obs_kind]["date"], format="mjd").datetime
            _ob["mjd"] = ob.loc[obs_kind]["date"]
            _ob["filter"] = ob.loc[obs_kind]["filter"]
            _ob["exptime"] = ob.loc[obs_kind]["exptime"]
            self.ob[obs_kind] = _ob

    def get(self, field, default=None, obs_kind=None):

        if obs_kind is None:
            obs_kind = self.obs_kind
        else:
            if obs_kind not in OBS_KIND:
                raise KeyError(
                    f"OpsimData obs_kind {obs_kind} not present in ob, "
                    f"must be in {OBS_KIND}."
                )
        ob = self.ob[obs_kind]

        if field not in ob and default is None:
            raise KeyError("OpsimData field %s not present in ob" % field)

        return ob.get(field, default)


def ObSeqData(config, base, value_type):
    """Returns the obseq data for a pointing."""
    pointing = galsim.config.GetInputObj("obseq_data", config, base, "OpSeqDataLoader")

    req = {"field": str}
    opt = {"obs_kind": str}

    kwargs, safe = galsim.config.GetAllParams(config, base, req=req, opt=opt)
    field = kwargs["field"]
    obs_kind = kwargs.get("obs_kind", None)

    val = value_type(pointing.get(field, obs_kind=obs_kind))
    return val, safe


RegisterInputType(
    "obseq_data", InputLoader(ObSeqDataLoader, file_scope=True, takes_logger=True)
)
RegisterValueType(
    "ObSeqData", ObSeqData, [float, int, str, Angle], input_type="obseq_data"
)
