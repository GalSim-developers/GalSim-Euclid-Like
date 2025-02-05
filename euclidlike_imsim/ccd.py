import galsim
import galsim.config
from galsim.config import RegisterImageType
from galsim.config.image_scattered import ScatteredImageBuilder
from galsim.errors import GalSimConfigValueError
from galsim.image import Image
from astropy.time import Time
import numpy as np

import euclidlike
from euclidlike.instrument_params import gain, saturation
from .noise import cfg_noise_key, parse_noise_config, get_noise


class EuclidlikeCCDImageBuilder(ScatteredImageBuilder):
    def setup(self, config, base, image_num, obj_num, ignore, logger):
        """Do the initialization and setup for building the image.

        This figures out the size that the image will be, but doesn't actually build it yet.

        Parameters:
            config:     The configuration dict for the image field.
            base:       The base configuration dict.
            image_num:  The current image number.
            obj_num:    The first object number in the image.
            ignore:     A list of parameters that are allowed to be in config that we can
                        ignore here. i.e. it won't be an error if these parameters are present.
            logger:     If given, a logger object to log progress.

        Returns:
            xsize, ysize
        """
        logger.debug(
            "image %d: Building EuclidlikeCCD: image, obj = %d,%d",
            image_num,
            image_num,
            obj_num,
        )

        self.nobjects = self.getNObj(config, base, image_num, logger=logger)
        logger.debug("image %d: nobj = %d", image_num, self.nobjects)

        # These are allowed for Scattered, but we don't use them here.
        extra_ignore = [
            "image_pos",
            "world_pos",
            "stamp_size",
            "stamp_xsize",
            "stamp_ysize",
            "nobjects",
        ]
        req = {"CCD": int, "filter": str, "mjd": float, "exptime": float}
        opt = {
            "draw_method": str,
        }
        opt.update({key: bool for key in cfg_noise_key})
        params = galsim.config.GetAllParams(
            config, base, req=req, opt=opt, ignore=ignore + extra_ignore
        )[0]

        self.ccd = params["CCD"]
        base["CCD"] = self.ccd
        self.filter = params["filter"]
        self.mjd = params["mjd"]
        self.exptime = params["exptime"]

        self.cfg_noise = parse_noise_config(params)

        # If draw_method isn't in image field, it may be in stamp.  Check.
        self.draw_method = params.get(
            "draw_method", base.get("stamp", {}).get("draw_method", "auto")
        )

        # If user hasn't overridden the bandpass to use, get the standard one.
        if "bandpass" not in config:
            base["bandpass"] = galsim.config.BuildBandpass(
                base["image"], "bandpass", base, logger=logger
            )

        return euclidlike.n_pix_col, euclidlike.n_pix_row

    def buildImage(self, config, base, image_num, obj_num, logger):
        """Build an Image containing multiple objects placed at arbitrary locations.

        Parameters:
            config:     The configuration dict for the image field.
            base:       The base configuration dict.
            image_num:  The current image number.
            obj_num:    The first object number in the image.
            logger:     If given, a logger object to log progress.

        Returns:
            the final image and the current noise variance in the image as a tuple
        """
        full_xsize = base["image_xsize"]
        full_ysize = base["image_ysize"]
        wcs = base["wcs"]

        full_image = Image(full_xsize, full_ysize, dtype=float)
        full_image.setOrigin(base["image_origin"])
        full_image.wcs = wcs
        full_image.setZero()

        full_image.header = galsim.FitsHeader()
        full_image.header["EXPTIME"] = self.exptime
        full_image.header["MJD-OBS"] = self.mjd
        full_image.header["DATE-OBS"] = str(
            Time(self.mjd, format="mjd").datetime
        )
        full_image.header["FILTER"] = self.filter
        full_image.header["GAIN"] = gain
        full_image.header["ZPTMAG"] = 2.5 * np.log10(
            self.exptime * euclidlike.collecting_area
        )

        base["current_image"] = full_image

        if "image_pos" in config and "world_pos" in config:
            raise GalSimConfigValueError(
                "Both image_pos and world_pos specified for Scattered image.",
                (config["image_pos"], config["world_pos"]),
            )

        if "image_pos" not in config and "world_pos" not in config:
            xmin = base["image_origin"].x
            xmax = xmin + full_xsize - 1
            ymin = base["image_origin"].y
            ymax = ymin + full_ysize - 1
            config["image_pos"] = {
                "type": "XY",
                "x": {"type": "Random", "min": xmin, "max": xmax},
                "y": {"type": "Random", "min": ymin, "max": ymax},
            }

        nbatch = self.nobjects // 1000 + 1
        for batch in range(nbatch):
            start_obj_num = self.nobjects * batch // nbatch
            end_obj_num = self.nobjects * (batch + 1) // nbatch
            nobj_batch = end_obj_num - start_obj_num
            if nbatch > 1:
                logger.warning(
                    "Start batch %d/%d with %d objects [%d, %d)",
                    batch + 1,
                    nbatch,
                    nobj_batch,
                    start_obj_num,
                    end_obj_num,
                )
            stamps, current_vars = galsim.config.BuildStamps(
                nobj_batch,
                base,
                logger=logger,
                obj_num=start_obj_num,
                do_noise=False,
            )
            base["index_key"] = "image_num"

            for k in range(nobj_batch):
                # This is our signal that the object was skipped.
                if stamps[k] is None:
                    continue
                bounds = stamps[k].bounds & full_image.bounds
                if not bounds.isDefined():  # pragma: no cover
                    # These noramlly show up as stamp==None, but technically it is possible
                    # to get a stamp that is off the main image, so check for that here to
                    # avoid an error.  But this isn't covered in the imsim test suite.
                    continue

                logger.debug(
                    "image %d: full bounds = %s",
                    image_num,
                    str(full_image.bounds),
                )
                logger.debug(
                    "image %d: stamp %d bounds = %s",
                    image_num,
                    k + start_obj_num,
                    str(stamps[k].bounds),
                )
                logger.debug("image %d: Overlap = %s", image_num, str(bounds))
                full_image[bounds] += stamps[k][bounds]
            stamps = None

        return full_image, None

    def addNoise(
        self, image, config, base, image_num, obj_num, current_var, logger
    ):
        """Add the final noise to a Scattered image

        Parameters:
            image:          The image onto which to add the noise.
            config:         The configuration dict for the image field.
            base:           The base configuration dict.
            image_num:      The current image number.
            obj_num:        The first object number in the image.
            current_var:    The current noise variance in each postage stamps.
            logger:         If given, a logger object to log progress.
        """
        # check ignore noise
        if not self.cfg_noise["use_noise"]:
            return

        if "noise_image" not in base.keys():
            get_noise(self.cfg_noise, config, base, logger)

        # We first have to apply gain and quantize the image
        image /= gain
        image.quantize()

        image += base["noise_image"]

        # Apply saturation
        saturation_ADU = np.round(saturation / gain)
        mask_saturated = image.array > saturation_ADU
        base["saturated_mask"] = mask_saturated
        image.array[mask_saturated] = saturation_ADU

        if self.cfg_noise["sky_subtract"]:
            image -= base["sky_image"]


# Register this as a valid type
RegisterImageType("euclidlike_ccd", EuclidlikeCCDImageBuilder())
