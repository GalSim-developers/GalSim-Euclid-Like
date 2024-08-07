import galsim
import euclidlike
import galsim.config
from galsim.config import RegisterObjectType, RegisterInputType, OpticalPSF, InputLoader

import euclidlike


class EuclidlikePSF(object):
    """Class building needed Euclidlike PSFs."""

    def __init__(
        self,
        CCD=None,
        WCS=None,
        n_waves=None,
        bpass=None,
        extra_aberrations=None,
        logger=None,
    ):

        logger = galsim.config.LoggerWrapper(logger)

        if n_waves == -1:
            if bpass.name == "W146":
                n_waves = 10
            else:
                n_waves = 5

        corners = [
            galsim.PositionD(1, 1),
            galsim.PositionD(1, euclidlike.n_pix_col),
            galsim.PositionD(euclidlike.n_pix_row, 1),
            galsim.PositionD(euclidlike.n_pix_row, euclidlike.n_pix_col),
        ]
        cc = galsim.PositionD(euclidlike.n_pix_row / 2, euclidlike.n_pix_col / 2)
        tags = ["ll", "lu", "ul", "uu"]
        self.PSF = {}
        pupil_bin = 8
        self.PSF[pupil_bin] = {}
        for tag, CCD_pos in tuple(zip(tags, corners)):
            self.PSF[pupil_bin][tag] = self._psf_call(
                CCD, bpass, CCD_pos, WCS, pupil_bin, n_waves, logger, extra_aberrations
            )
        for pupil_bin in [4, 2, "achromatic"]:
            self.PSF[pupil_bin] = self._psf_call(
                CCD, bpass, cc, WCS, pupil_bin, n_waves, logger, extra_aberrations
            )

    def _parse_pupil_bin(self, pupil_bin):
        if pupil_bin == "achromatic":
            return 8
        else:
            return pupil_bin

    def _psf_call(
        self, CCD, bpass, CCD_pos, WCS, pupil_bin, n_waves, logger, extra_aberrations
    ):

        if pupil_bin == 8:
            psf = euclidlike.getPSF(
                CCD,
                bpass.name,
                ccd_pos=CCD_pos,
                wcs=WCS,
                # pupil_bin=pupil_bin,
                # n_waves=n_waves,
                logger=logger,
                # Don't set wavelength for this one.
                # We want this to be chromatic for photon shooting.
                # wavelength          = bpass.effective_wavelength,
                # extra_aberrations=extra_aberrations,
            )
        else:
            psf = euclidlike.getBrightPSF(
                CCD,
                bpass.name,
                ccd_pos=CCD_pos,
                wcs=WCS,
                pupil_bin=self._parse_pupil_bin(pupil_bin),
                n_waves=n_waves,
                logger=logger,
                # Note: setting wavelength makes it achromatic.
                # We only use pupil_bin = 2,4 for FFT objects.
                wavelength=bpass.effective_wavelength,
            )

            # psf = euclidlike.getPSF(
            #     CCD,
            #     bpass.name,
            #     ccd_pos=CCD_pos,
            #     wcs=WCS,
            #     # pupil_bin=self._parse_pupil_bin(pupil_bin),
            #     # n_waves=n_waves,
            #     logger=logger,
            #     # Note: setting wavelength makes it achromatic.
            #     # We only use pupil_bin = 2,4 for FFT objects.
            #     wavelength=bpass.effective_wavelength,
            #     # extra_aberrations=extra_aberrations,
            # )
        if pupil_bin == 4:
            return psf.withGSParams(maximum_fft_size=16384, folding_threshold=1e-3)
        elif pupil_bin == 2:
            return psf.withGSParams(maximum_fft_size=16384, folding_threshold=1e-4)
        else:
            return psf.withGSParams(maximum_fft_size=16384)

    def getPSF(self, pupil_bin, pos):
        """
        Return a PSF to be convolved with sources.

        @param [in] what pupil binning to request.
        """

        # temporary
        # psf = self.PSF[pupil_bin]['ll']
        # if ((pos.x-roman.n_pix)**2+(pos.y-roman.n_pix)**2)<((pos.x-1)**2+(pos.y-1)**2):
        #     psf = self.PSF[pupil_bin]['uu']
        # if ((pos.x-1)**2+(pos.y-roman.n_pix)**2)<((pos.x-roman.n_pix)**2+(pos.y-roman.n_pix)**2):
        #     psf = self.PSF[pupil_bin]['lu']
        # if ((pos.x-roman.n_pix)**2+(pos.y-1)**2)<((pos.x-1)**2+(pos.y-roman.n_pix)**2):
        #     psf = self.PSF[pupil_bin]['ul']
        # if ((pos.x-roman.n_pix/2)**2+(pos.y-roman.n_pix/2)**2)<((pos.x-roman.n_pix)**2+(pos.y-1)**2):
        #     psf = self.PSF[pupil_bin]['cc']
        # return psf

        psf = self.PSF[pupil_bin]
        if pupil_bin != 8:
            return psf

        wll = (euclidlike.n_pix_col - pos.x) * (euclidlike.n_pix_row - pos.y)
        wlu = (euclidlike.n_pix_col - pos.x) * (pos.y - 1)
        wul = (pos.x - 1) * (euclidlike.n_pix_row - pos.y)
        wuu = (pos.x - 1) * (pos.y - 1)
        return (
            wll * psf["ll"] + wlu * psf["lu"] + wul * psf["ul"] + wuu * psf["uu"]
        ) / ((euclidlike.n_pix_row - 1) * (euclidlike.n_pix_col - 1))


class PSFLoader(InputLoader):
    """PSF loader."""

    def __init__(self):
        # Override some defaults in the base init.
        super().__init__(init_func=EuclidlikePSF, takes_logger=True, use_proxy=False)

    def getKwargs(self, config, base, logger):
        logger.debug("Get kwargs for PSF")

        req = {}
        opt = {
            "n_waves": int,
        }
        ignore = ["extra_aberrations"]

        # If CCD is in base, then don't require it in the config file.
        # (Presumably because using Euclidlike image type, which sets it there for convenience.)
        if "CCD" in base:
            opt["CCD"] = int
        else:
            req["CCD"] = int

        kwargs, safe = galsim.config.GetAllParams(
            config, base, req=req, opt=opt, ignore=ignore
        )

        # If not given in kwargs, then it must have been in base, so this is ok.
        if "CCD" not in kwargs:
            kwargs["CCD"] = base["CCD"]

        kwargs["extra_aberrations"] = galsim.config.ParseAberrations(
            "extra_aberrations", config, base, "EuclidlikePSF"
        )
        kwargs["WCS"] = galsim.config.BuildWCS(
            base["image"], "wcs", base, logger=logger
        )
        kwargs["bpass"] = galsim.config.BuildBandpass(
            base["image"], "bandpass", base, logger
        )[0]

        logger.debug("kwargs = %s", kwargs)

        return kwargs, False


# Register this as a valid type
RegisterInputType("euclidlike_psf", PSFLoader())
# RegisterObjectType('euclidlike_psf', BuildEuclidlikePSF, input_type='euclidlikepsf_loader')
