import numpy as np

from astropy.time import Time

import galsim
from galsim.errors import GalSimConfigError

import euclidlike
from euclidlike.instrument_params import gain, read_noise

# Some place holder stuff, will need to update later in euclidlike
# This is just to avoid changing the code and don't have "errors"
stray_light_fraction = 0.
dark_current = 0.

cfg_noise_key = [
    "ignore_noise",
    "stray_light",
    "thermal_background",
    "reciprocity_failure",
    "dark_current",
    "nonlinearity",
    "ipc",
    "read_noise",
    "sky_subtract",
]


def parse_noise_config(params):

    cfg_noise = {}
    for key in cfg_noise_key:
        cfg_noise[key] = params.get(key, False)

    cfg_noise["exptime"] = params["exptime"]
    cfg_noise["mjd"] = params["mjd"]

    return cfg_noise


def get_noise(cfg_noise, cfg_image, base, logger):

    noise_img = base["current_image"].copy()
    noise_img.fill(0)

    wcs = base["wcs"]
    bp = base["bandpass"]
    rng = galsim.config.GetRNG(cfg_image, base)
    logger.info(
        "image %d: Start EuclidlikeCCD detector effects", base.get("image_num", 0)
    )

    # Things that will eventually be subtracted (if sky_subtract) will have their expectation
    # value added to sky_image.  So technically, this includes things that aren't just sky.
    # E.g. includes dark_current and thermal backgrounds.
    sky_image = noise_img.copy()
    sky_level = euclidlike.getSkyLevel(
        bp,
        world_pos=wcs.toWorld(noise_img.true_center),
        date=Time(cfg_noise["mjd"], format="mjd").datetime,
    )
    logger.debug("Adding sky_level = %s", sky_level)
    if cfg_noise["stray_light"]:
        logger.debug("Stray light fraction = %s", stray_light_fraction)
        sky_level *= 1.0 + stray_light_fraction
    wcs.makeSkyImage(sky_image, sky_level)

    # The image up to here is an expectation value.
    # Realize it as an integer number of photons.
    poisson_noise = galsim.noise.PoissonNoise(rng)
    logger.debug("Adding poisson noise to sky photons")
    sky_image1 = sky_image.copy()
    sky_image1.addNoise(poisson_noise)
    noise_img.quantize()
    # the image won't necessarily be integers.
    noise_img += sky_image1

    # Apply the detector effects here.  Not all of these are "noise" per se, but they
    # happen interspersed with various noise effects, so apply them all in this step.

    if cfg_noise["dark_current"]:
        dc = dark_current * cfg_noise["exptime"]
        logger.debug("Adding dark current: %s", dc)
        sky_image += dc
        dark_noise = galsim.noise.DeviateNoise(
            galsim.random.PoissonDeviate(rng, dc)
        )
        noise_img.addNoise(dark_noise)

    if cfg_noise["read_noise"]:
        logger.debug("Adding read noise %s", read_noise)
        noise_img.addNoise(galsim.GaussianNoise(rng, sigma=read_noise))

    logger.debug("Applying gain %s", gain)
    noise_img /= gain

    # Make integer ADU now.
    noise_img.quantize()

    sky_image /= gain
    sky_image.quantize()

    base["noise_image"] = noise_img
    base["sky_image"] = sky_image


class NoiseImageBuilder(galsim.config.ExtraOutputBuilder):

    def initialize(self, data, scratch, config, base, logger):
        """Do any initial setup for this builder at the start of a new output file.

        The base class implementation saves two work space items into self.data and self.scratch
        that can be used to safely communicate across multiple processes.

        Parameters:
            data:       An empty list of length nimages to use as work space.
            scratch:    An empty dict that can be used as work space.
            config:     The configuration field for this output object.
            base:       The base configuration dict.
            logger:     If given, a logger object to log progress. [default: None]
        """
        req = {"CCD": int, "filter": str, "mjd": float, "exptime": float}
        opt = {key: bool for key in cfg_noise_key}
        ignore = galsim.config.image.image_ignore
        extra_ignore = [
            "image_pos",
            "world_pos",
            "stamp_size",
            "stamp_xsize",
            "stamp_ysize",
            "nobjects",
        ]
        params = galsim.config.GetAllParams(
            base["image"], base, req=req, opt=opt, ignore=ignore+extra_ignore,
        )[0]
        self.cfg_noise = parse_noise_config(params)

        self._check_input()

        self.data = data
        self.scratch = scratch
        self.final_data = None

    def _check_input(self):
        if self.cfg_noise["ignore_noise"]:
            raise GalSimConfigError(
                "You cannot ignore the noise and request the noise image at the same time."
                " Either active the noise or remove the output noise image."
            )

    def processImage(self, index, obj_nums, config, base, logger):
        """Perform any necessary processing at the end of each image construction.

        This function will be called after each full image is built.

        Compute the noise for the current image and add it to the image.
        It will also create an independent noise image.
        The code optionally subtract the background if requested.

        Parameters:
            index:      The index in self.data to use for this image.  This isn't the image_num
                        (which can be accessed at base['image_num'] if needed), but rather
                        an index that starts at 0 for the first image being worked on and
                        goes up to nimages-1.
            obj_nums:   The object numbers that were used for this image.
            config:     The configuration field for this output object.
            base:       The base configuration dict.
            logger:     If given, a logger object to log progress. [default: None]
        """
        if "noise_image" not in base.keys():
            get_noise(self.cfg_noise, config, base, logger)

        noise_image = base["noise_image"].copy()

        if self.cfg_noise["sky_subtract"]:
            noise_image -= base["sky_image"]

        self.data[index] = noise_image


class SkyImageBuilder(NoiseImageBuilder):

    def _check_input(self):
        if self.cfg_noise["ignore_noise"]:
            raise GalSimConfigError(
                "You cannot ignore the noise and request the sky image at the same time."
                " Either activate the noise or remove the output sky image."
            )

    def processImage(self, index, obj_nums, config, base, logger):
        """Perform any necessary processing at the end of each image construction.

        This function will be called after each full image is built.

        Compute the sky background and return it in an image.

        Parameters:
            index:      The index in self.data to use for this image.  This isn't the image_num
                        (which can be accessed at base['image_num'] if needed), but rather
                        an index that starts at 0 for the first image being worked on and
                        goes up to nimages-1.
            obj_nums:   The object numbers that were used for this image.
            config:     The configuration field for this output object.
            base:       The base configuration dict.
            logger:     If given, a logger object to log progress. [default: None]
        """
        if "noise_image" not in base.keys():
            get_noise(self.cfg_noise, config, base, logger)

        self.data[index] = base["sky_image"].copy()


class WeightImageBuilder(NoiseImageBuilder):

    def _check_input(self):
        if self.cfg_noise["ignore_noise"]:
            raise GalSimConfigError(
                "You cannot ignore the noise and request the weight image at the same time."
                " Either activate the noise or remove the output sky image."
            )

    def processImage(self, index, obj_nums, config, base, logger):
        """Perform any necessary processing at the end of each image construction.

        This function will be called after each full image is built.

        Compute the weight map from the noise image and return it in an image.

        Parameters:
            index:      The index in self.data to use for this image.  This isn't the image_num
                        (which can be accessed at base['image_num'] if needed), but rather
                        an index that starts at 0 for the first image being worked on and
                        goes up to nimages-1.
            obj_nums:   The object numbers that were used for this image.
            config:     The configuration field for this output object.
            base:       The base configuration dict.
            logger:     If given, a logger object to log progress. [default: None]
        """
        if "noise_image" not in base.keys():
            get_noise(self.cfg_noise, config, base, logger)

        noise_var = np.var(base["noise_image"].array)

        weight_image = galsim.ImageF(base["image_bounds"])
        weight_image += noise_var
        weight_image.invertSelf()

        self.data[index] = weight_image


galsim.config.RegisterExtraOutput('noise_image', NoiseImageBuilder())
galsim.config.RegisterExtraOutput('sky_image', SkyImageBuilder())
galsim.config.RegisterExtraOutput('weight_image', WeightImageBuilder())
