"""
Interface to obtain objects from skyCatalogs.
"""

import os
import numpy as np
import galsim
from galsim.config import (
    InputLoader,
    RegisterInputType,
    RegisterValueType,
    RegisterObjectType,
)

import euclidlike


class SkyCatalogInterface:
    """Interface to skyCatalogs package."""

    _trivial_sed = galsim.SED(
        galsim.LookupTable([100, 2600], [1, 1], interpolant="linear"),
        wave_type="nm",
        flux_type="fphotons",
    )

    def __init__(
        self,
        file_name,
        exptime,
        wcs=None,
        mjd=None,
        bandpass=None,
        xsize=None,
        ysize=None,
        obj_types=None,
        edge_pix=100,
        max_flux=None,
        logger=None,
    ):
        """
        Parameters
        ----------
        file_name : str
            Name of skyCatalogs yaml config file.
        wcs : galsim.WCS
            WCS of the image to render.
        mjd : float
            MJD of the midpoint of the exposure.
        exptime : float
            Exposure time.
        xsize : int
            Size in pixels of CCD in x-direction.
        ysize : int
            Size in pixels of CCD in y-direction.
        obj_types : list-like [None]
            List or tuple of object types to render, e.g., ('star', 'galaxy').
            If None, then consider all object types.
        edge_pix : float [100]
            Size in pixels of the buffer region around nominal image
            to consider objects.
        logger : logging.Logger [None]
            Logger object.
        """
        self.file_name = file_name
        self.wcs = wcs
        self.mjd = mjd
        self.exptime = exptime
        self.bandpass = bandpass
        if xsize is not None:
            self.xsize = xsize
        else:
            self.xsize = euclidlike.n_pix_col
        if ysize is not None:
            self.ysize = ysize
        else:
            self.ysize = euclidlike.n_pix_row
        self.obj_types = obj_types
        self.edge_pix = edge_pix
        self.logger = galsim.config.LoggerWrapper(logger)

        if obj_types is not None:
            self.logger.warning(f"Object types restricted to {obj_types}")
        self.ccd_center = wcs.toWorld(
            galsim.PositionD(self.xsize / 2.0, self.ysize / 2.0)
        )
        self._objects = None

    @property
    def objects(self):
        from skycatalogs import skyCatalogs

        if self._objects is None:

            # Select objects from polygonal region bounded by CCD edges
            corners = (
                (-self.edge_pix, -self.edge_pix),
                (self.xsize + self.edge_pix, -self.edge_pix),
                (self.xsize + self.edge_pix, self.ysize + self.edge_pix),
                (-self.edge_pix, self.ysize + self.edge_pix),
            )
            vertices = []
            for x, y in corners:
                sky_coord = self.wcs.toWorld(galsim.PositionD(x, y))
                vertices.append(
                    (sky_coord.ra / galsim.degrees, sky_coord.dec / galsim.degrees)
                )
            region = skyCatalogs.PolygonalRegion(vertices)
            sky_cat = skyCatalogs.open_catalog(self.file_name)
            self._objects = sky_cat.get_objects_by_region(
                region, obj_type_set=self.obj_types, mjd=self.mjd
            )
            if not self._objects:
                self.logger.warning("No objects found on image.")
            else:
                self._build_dtype_dict()
        return self._objects

    def _build_dtype_dict(self):
        self._dtype_dict = {}
        obj_types = []
        for coll in self._objects.get_collections():
            objects_type = coll._object_type_unique
            if objects_type in obj_types:
                continue
            col_names = list(coll.native_columns)
            for col_name in col_names:
                try:
                    # Some columns cannot be read in snana
                    np_type = coll.get_native_attribute(col_name).dtype.type()
                except Exception as e:
                    self.logger.warning(
                        f"The column {col_name} could not be read from skyCatalog."
                    )
                    continue
                if np_type is None:
                    py_type = str
                else:
                    py_type = type(np_type.astype(object))
                self._dtype_dict[col_name] = py_type

    def get_ccd_center(self):
        """
        Return the CCD center.
        """
        return self.ccd_center

    def getNObjects(self):
        """
        Return the number of GSObjects to render
        """
        return len(self.objects)

    def getApproxNObjects(self):
        """
        Return the approximate number of GSObjects to render, as set in
        the class initializer.
        """
        return self.getNObjects()

    def getWorldPos(self, index):
        """
        Return the sky coordinates of the skyCatalog object
        corresponding to the specified index.

        Parameters
        ----------
        index : int
            Index of the (object_index, subcomponent) combination.

        Returns
        -------
        galsim.CelestialCoord
        """
        skycat_obj = self.objects[index]
        ra, dec = skycat_obj.ra, skycat_obj.dec
        return galsim.CelestialCoord(ra * galsim.degrees, dec * galsim.degrees)

    def getFlux(self, index, filter=None, mjd=None, exptime=None):

        if filter is None:
            filter = self.bandpass.name
        if mjd is None:
            mjd = self.mjd
        if exptime is None:
            exptime = self.exptime

        skycat_obj = self.objects[index]
        # We cache the SEDs for potential later use
        self._seds = skycat_obj.get_observer_sed_components()
        for i, sed in enumerate(self._seds.values()):
            if i == 0:
                sed_sum = sed
            else:
                sed_sum += sed
        raw_flux = skycat_obj.get_euclid_flux(
            filter,
            sed_sum,
            mjd=mjd,
            cache=False
        )
        if hasattr(skycat_obj, "get_wl_params"):
            _, _, mu = skycat_obj.get_wl_params()
        else:
            mu = 1.
        # raw_flux = skycat_obj.get_euclid_flux(filter, mjd=mjd, cache=False)
        flux = raw_flux * mu * exptime * euclidlike.collecting_area

        return flux
    def getValue(self, index, field):

        skycat_obj = self.objects[index]

        if field not in self._dtype_dict:
            # We cannot raise an error because one could have a field for snana
            # in the config and we don't want to crash because there are no SN
            # in this particular image. We then default to False which might not
            # be the right type for the required column but we have no way of knowing
            # the correct type if the column do not exist.
            self.logger.warning(f"The field {field} was not found in skyCatalog.")
            return None
        elif field not in skycat_obj.native_columns:
            if self._dtype_dict[field] is int:
                # There are no "special value" for integer so we default to
                # hopefully something completely off
                return -9999
            elif self._dtype_dict[field] is float:
                return np.nan
            elif self._dtype_dict[field] is str:
                return None
        else:
            return skycat_obj.get_native_attribute(field)

    def getObj(self, index, gsparams=None, rng=None, exptime=30):
        """
        Return the galsim object for the skyCatalog object
        corresponding to the specified index.  If the skyCatalog
        object is a galaxy, the returned galsim object will be
        a galsim.Sum.

        Parameters
        ----------
        index : int
            Index of the object in the self.objects catalog.

        Returns
        -------
        galsim.GSObject
        """
        if not self.objects:
            raise RuntimeError("Trying to get an object from an empty sky catalog")

        faint = False
        skycat_obj = self.objects[index]
        gsobjs = skycat_obj.get_gsobject_components(gsparams)

        # Compute the flux or get the cached value.
        flux = self.getFlux(index)
        if np.isnan(flux):
            return None

        # Set up simple SED if too faint
        if flux < 40:
            faint = True
        if not faint:
            seds = skycat_obj.get_observer_sed_components(mjd=self.mjd)

        gs_obj_list = []
        for component in gsobjs:
            if faint:
                gsobjs[component] = gsobjs[component].evaluateAtWavelength(
                    self.bandpass
                )
                gs_obj_list.append(
                    gsobjs[component]
                    * self._trivial_sed
                    * self.exptime
                    * euclidlike.collecting_area
                )
            else:
                if component in seds:
                    gs_obj_list.append(
                        gsobjs[component]
                        * seds[component]
                        * self.exptime
                        * euclidlike.collecting_area
                    )

        if not gs_obj_list:
            return None

        if len(gs_obj_list) == 1:
            gs_object = gs_obj_list[0]
        else:
            gs_object = galsim.Add(gs_obj_list)

        # Give the object the right flux
        gs_object.flux = flux
        gs_object.withFlux(gs_object.flux, self.bandpass)

        # Get the object type
        if (skycat_obj.object_type == "diffsky_galaxy") | (
            skycat_obj.object_type == "galaxy"
        ):
            gs_object.object_type = "galaxy"
        if skycat_obj.object_type == "star":
            gs_object.object_type = "star"
        if skycat_obj.object_type == "snana":
            gs_object.object_type = "transient"

        return gs_object


class SkyCatalogLoader(InputLoader):
    """
    Class to load SkyCatalogInterface object.
    """

    def getKwargs(self, config, base, logger):
        req = {"file_name": str, "exptime": float}
        opt = {
            "edge_pix": float,
            "obj_types": list,
            "mjd": float,
        }
        kwargs, safe = galsim.config.GetAllParams(config, base, req=req, opt=opt)
        wcs = galsim.config.BuildWCS(base["image"], "wcs", base, logger=logger)
        kwargs["wcs"] = wcs
        kwargs["logger"] = logger

        if "bandpass" not in config:
            base["bandpass"] = galsim.config.BuildBandpass(
                base["image"], "bandpass", base, logger=logger
            )[0]

        kwargs["bandpass"] = base["bandpass"]
        # Sky catalog object lists are created per CCD, so they are
        # not safe to reuse.
        safe = False
        return kwargs, safe


def SkyCatObj(config, base, ignore, gsparams, logger):
    """
    Build an object according to info in the sky catalog.
    """
    skycat = galsim.config.GetInputObj("sky_catalog", config, base, "SkyCatObj")

    # Ensure that this sky catalog matches the CCD being simulated by
    # comparing center locations on the sky.
    world_center = base["world_center"]
    ccd_center = skycat.get_ccd_center()
    sep = ccd_center.distanceTo(base["world_center"]) / galsim.arcsec
    # Centers must agree to within at least 1 arcsec:
    if sep > 1.0:
        message = (
            "skyCatalogs selection and CCD center do not agree: \n"
            "skycat.ccd_center: "
            f"{ccd_center.ra/galsim.degrees:.5f}, "
            f"{ccd_center.dec/galsim.degrees:.5f}\n"
            "world_center: "
            f"{world_center.ra/galsim.degrees:.5f}, "
            f"{world_center.dec/galsim.degrees:.5f} \n"
            f"Separation: {sep:.2e} arcsec"
        )
        raise RuntimeError(message)

    # Setup the indexing sequence if it hasn't been specified.  The
    # normal thing with a catalog is to just use each object in order,
    # so we don't require the user to specify that by hand.  We can do
    # it for them.
    galsim.config.SetDefaultIndex(config, skycat.getNObjects())

    req = {"index": int}
    opt = {"num": int}
    kwargs, safe = galsim.config.GetAllParams(config, base, req=req, opt=opt)
    index = kwargs["index"]

    rng = galsim.config.GetRNG(config, base, logger, "SkyCatObj")

    obj = skycat.getObj(index, gsparams=gsparams, rng=rng)
    base["object_id"] = skycat.objects[index].id

    return obj, safe


def SkyCatWorldPos(config, base, value_type):
    """Return a value from the object part of the skyCatalog"""
    skycat = galsim.config.GetInputObj("sky_catalog", config, base, "SkyCatWorldPos")

    # Setup the indexing sequence if it hasn't been specified.  The
    # normal thing with a catalog is to just use each object in order,
    # so we don't require the user to specify that by hand.  We can do
    # it for them.
    galsim.config.SetDefaultIndex(config, skycat.getNObjects())

    req = {"index": int}
    opt = {"num": int}
    kwargs, safe = galsim.config.GetAllParams(config, base, req=req, opt=opt)
    index = kwargs["index"]

    pos = skycat.getWorldPos(index)
    return pos, safe


def SkyCatValue(config, base, value_type):

    skycat = galsim.config.GetInputObj("sky_catalog", config, base, "SkyCatValue")

    # Setup the indexing sequence if it hasn't been specified.  The
    # normal thing with a catalog is to just use each object in order,
    # so we don't require the user to specify that by hand.  We can do
    # it for them.
    galsim.config.SetDefaultIndex(config, skycat.getNObjects())

    req = {"field": str, "index": int}
    opt = {"obs_kind": str}
    params, safe = galsim.config.GetAllParams(config, base, req=req, opt=opt)
    field = params["field"]
    index = params["index"]
    obs_kind = params.get("obs_kind", None)

    if field == "flux":
        if obs_kind is None:
            val = skycat.getFlux(index)
        else:
            pointing = galsim.config.GetInputObj("obseq_data", config, base, "OpSeqDataLoader")
            filter = pointing.get("filter", obs_kind=obs_kind)
            exptime = pointing.get("exptime", obs_kind=obs_kind)
            mjd = pointing.get("mjd", obs_kind=obs_kind)
            val = skycat.getFlux(index, filter=filter, exptime=exptime, mjd=mjd)
    else:
        val = skycat.getValue(index, field)

    return val, safe



RegisterInputType("sky_catalog", SkyCatalogLoader(SkyCatalogInterface, has_nobj=True))
RegisterObjectType("SkyCatObj", SkyCatObj, input_type="sky_catalog")
RegisterValueType(
    "SkyCatWorldPos", SkyCatWorldPos, [galsim.CelestialCoord], input_type="sky_catalog"
)

# Here we have to provide None as a type otherwise Galsim complains but I don't know why..
RegisterValueType(
    "SkyCatValue", SkyCatValue, [float, int, str, None]  # , input_type="sky_catalog"
)
