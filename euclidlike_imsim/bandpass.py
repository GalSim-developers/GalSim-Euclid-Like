import euclidlike
from galsim.config import BandpassBuilder, RegisterBandpassType, GetAllParams


class EuclidlikeBandpassBuilder(BandpassBuilder):
    """A class for loading a Bandpass from a file

    FileBandpass expected the following parameter:

        name (str)          The name of the Euclid filter to get. (required)
    """
    def buildBandpass(self, config, base, logger):
        """Build the Bandpass based on the specifications in the config dict.

        Parameters:
            config:     The configuration dict for the bandpass type.
            base:       The base configuration dict.
            logger:     If provided, a logger for logging debug statements.

        Returns:
            the constructed Bandpass object.
        """
        req = {'name': str}
        kwargs, safe = GetAllParams(config, base, req=req)

        name = kwargs['name']
        # Hard set the limit due to PSF definition
        bandpass = euclidlike.getBandpasses()[name]
        if name == "VIS":
            bandpass.blue_limit = 540
            bandpass.red_limit = 910

        return bandpass, safe


RegisterBandpassType('EuclidlikeBandpass', EuclidlikeBandpassBuilder())
RegisterBandpassType('EuclidlikeBandpassTrimmed', EuclidlikeBandpassBuilder())
