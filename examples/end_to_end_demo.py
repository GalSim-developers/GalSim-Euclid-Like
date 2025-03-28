"""
@file end_to_end_demo.py

This script is an adaptation of the demo13 script in Galsim for Euclid. The script is intended to
produce a relatively realistic scene of galaxies and stars as will be observed by the
Euclid Space Telescope, including the Euclid-like PSF, WCS, and various detector effects.

The demo13.py script can be found here:

https://github.com/GalSim-developers/GalSim/blob/releases/2.5/examples/demo13.py

The script uses the Euclid-like module to set up the approximate filters, PSF, and WCS for the
Euclid Space Telescope. It also uses the COSMOSCatalog class to read in the COSMOS catalog
of galaxy properties.

The current package misses the following features:

1) Non-linearity: Charge-dependent gain in converting from units of electrons to ADU. Non-linearity
   in some form is also relevant for CCDs in addition to NIR detectors.
2) And any non-linear effects that are specific to the Euclid detectors such as
   charge-transfer inefficiency.

It also uses chromatic photon shooting, which is generally a more efficient way to simulate
scenes with many faint galaxies.  The default FFT method for drawing chromatically is fairly
slow, since it needs to integrate the image over the bandpass.  With photon shooting, the
photons are assigned wavelengths according to the SED of the galaxy, and then each photon has
the appropriate application of the chromatic PSF according to the wavelength.

To select a different set of Euclid filters, you may use the `filters` option on the command line.
E.g. `python end_to_end_demo.py --filters=VIS`
"""

import argparse
import numpy as np
import sys
import os
import logging
import galsim
import euclidlike
import datetime


def parse_args(argv):
    parser = argparse.ArgumentParser(prog='end_to_end_demo', add_help=True)
    parser.add_argument('-f', '--filters', nargs='+', action='store', default="VIS",
                        help='Which filters to simulate')
    parser.add_argument('-o', '--outpath', type=str, default='output',
                        help='Which directory to put the output files')
    parser.add_argument('-n', '--nobj', type=int, default=1000,
                        help='How many objects to draw')
    parser.add_argument('--seed', type=int, default=12345,
                        help='Initial seed for random numbers')
    parser.add_argument('-c', '--ccd', type=int, default=8, choices=range(0,35),
                        help='Which ccd to simulate (default is arbitrarily CCD 8)')
    parser.add_argument('-t', '--test', action='store_const', default=False, const=True,
                        help='Whether to use the smaller test sample, rather than the full COSMOS samples')
    parser.add_argument('-v', '--verbosity', type=int, default=2, choices=range(0, 4),
                        help='Verbosity level')

    args = parser.parse_args(argv)
    return args


def main(argv):

    args = parse_args(argv)
    use_filters = args.filters
    outpath = args.outpath
    nobj = args.nobj
    seed = args.seed
    seed1 = galsim.BaseDeviate(seed).raw()
    use_CCD = args.ccd

    # Make output directory if not already present.
    if not os.path.isdir(outpath):
        os.mkdir(outpath)

    # Use a logger to output some information about the run.
    logging.basicConfig(format="%(message)s", stream=sys.stdout)
    logger = logging.getLogger("demo13")
    logging_levels = {
        0: logging.CRITICAL,
        1: logging.WARNING,
        2: logging.INFO,
        3: logging.DEBUG,
    }
    level = logging_levels[args.verbosity]
    logger.setLevel(level)

    # Read in the Euclid filters, setting an AB zeropoint appropriate for this telescope given its
    # collecting area and (since we didn't use any keyword arguments to modify this) using the typical
    # exposure time for Euclid images.  By default, this routine truncates the parts of the
    # bandpasses that are near 0 at the edges, and thins them by the default amount.
    euclidlike_filters = euclidlike.getBandpasses(AB_zeropoint=True)
    logger.debug('Read in Euclid imaging filters.')

    # Get the names of the ones we will use here.
    filters = [filter_name for filter_name in euclidlike_filters if filter_name in use_filters]
    logger.debug('Using filters: %s', filters)

    # We'll use this one for our flux normalization of stars, so we'll need this regardless of
    # which bandpass we are simulating.
    vis_bandpass = euclidlike_filters['VIS']
    vis_bandpass.red_limit = 910
    vis_bandpass.blue_limit = 540

    # Note: This example uses both the I<23.5 and I<25.2 COSMOS catalogs to try to better span
    #       a range from bigger, bright galaxies to fainter ones.  We also dilate and magnify
    #       the bright galaxies to make them a bit more visually compelling in this example.
    #
    #       You should download both of these COSMOS catalogs if you haven't yet.
    #
    #       We recommend specifying the directory for the download, rather than let it use the
    #       default directory, since that will be in the GalSim python share directory, which will
    #       be overwritten whenever you reinstall GalSim.  This command sets up a symlink from that
    #       location to a directory in your home directory.  (Feel free to use any other convenient
    #       directory of course, depending on your situation.)  You still need to rerun this command
    #       after reinstalls of GalSim, but it will just need to update the link, not actually
    #       re-download anything.
    #
    #           galsim_download_cosmos -s 23.5 -d ~/share
    #           galsim_download_cosmos -s 25.2 -d ~/share
    #
    # The area and exposure time here rescale the fluxes to be appropriate for the Euclid collecting
    # area and exposure time, rather than the default HST collecting area and 1 second exposures.
    #
    # If you really want to use the smaller test sample, you can use --sample test, but there
    # are only 100 galaxies there, so most galaxies will be repeated several times.

    if args.test:
        cat1 = galsim.COSMOSCatalog('real_galaxy_catalog_23.5_example.fits', dir='data',
                                    area=euclidlike.collecting_area, exptime=euclidlike.long_exptime)
        cat2 = cat1
    else:
        cat1 = galsim.COSMOSCatalog(sample='25.2', area=euclidlike.collecting_area, exptime=euclidlike.long_exptime)
        cat2 = galsim.COSMOSCatalog(sample='23.5', area=euclidlike.collecting_area, exptime=euclidlike.long_exptime)

    logger.info('Read in %d galaxies from I<25.2 catalog' % cat1.nobjects)
    logger.info('Read in %d galaxies from I<23.5 catalog' % cat2.nobjects)

    # For the stars, we'll use the vega SED, since that's the only stellar SED we have in the
    # GalSim share directory.  Which means all the stars will be pretty blue.
    vega_sed = galsim.SED('vega.txt', 'nm', 'flambda')

    # Pick a plausible observation (Center of the Euclid Deep Field North (EDF-N))
    ra_targ = galsim.Angle.from_hms('17:45:40')
    dec_targ = galsim.Angle.from_dms('66:0:0')
    targ_pos = galsim.CelestialCoord(ra=ra_targ, dec=dec_targ)

    # Get the WCS for an observation at this position.
    date = datetime.datetime(2025, 5, 16)

    # The output of this routine is a dict of WCS objects, one for each CCD. We then take the WCS
    # for the CCD that we are using.
    wcs_dict = euclidlike.getWCS(world_pos=targ_pos, CCDs=use_CCD, date=date)
    wcs = wcs_dict[use_CCD]

    # Now start looping through the filters to draw.
    for ifilter, filter_name in enumerate(filters):

        logger.info('Beginning work for {0}.'.format(filter_name))

        # GalSim uses the term Bandpass for the class that defines the throughput across a
        # filter bandpass, partly because "filter" is a reserved word in python.  So we follow
        # that convention here as well.
        bandpass = euclidlike_filters[filter_name]
        bandpass.red_limit = 910
        bandpass.blue_limit = 540

        # Create the PSF
        # We are ignoring the position-dependence of the PSF within each CCD, just using the PSF
        # at the center of the sensor.
        logger.info('Building PSF for CCD %d, filter %s.'%(use_CCD, filter_name))
        psf = euclidlike.getPSF(use_CCD, filter_name, wcs=wcs)
        lam = 700     # nm
        diam = 1.2    # meters
        lam_over_diam = (lam * 1.e-9) / diam  # radians
        lam_over_diam *= 206265  # Convert to arcsec
        #psf = galsim.Airy(lam_over_diam)

        # Set up the full image for the galaxies
        full_image = galsim.ImageF(euclidlike.n_pix_col, euclidlike.n_pix_row, wcs=wcs)

        # Also separately build up the sky image, which we need to get the noise right,
        # even though we'll subtract off the expectation of the sky image.
        sky_image = galsim.ImageF(euclidlike.n_pix_col, euclidlike.n_pix_row, wcs=wcs)

        # We have one rng for image-level stuff, and two others for the stars and galaxies.
        # There are simpler ways to do this in a python script (e.g. probably only need 2
        # rngs, not 3), but this way of setting it up matches the way the config file initializes
        # the random number generators.
        # Also, note that the second seed given in the config file, doesn't get the
        # BaseDeviate(...).raw() treatment.  Only the first item, which parses as an int.
        # When a random_seed config item is already a dict, GalSim leaves it as is.
        image_rng = galsim.UniformDeviate(seed1 + ifilter * nobj)

        # Start with the flux from the sky. This is a little easier to do first before adding
        # the light from the objects, since we will have to apply Poisson noise to the sky flux
        # manually, but the photon shooting will automatically include Poisson noise for the
        # objects.

        # First we get the amount of zodaical light for a position corresponding to the center of
        # this CCD.  The results are provided in units of e-/arcsec^2, using the default Euclid
        # exposure time since we did not explicitly specify one.  Then we multiply this by a factor
        # >1 to account for the amount of stray light that is expected.  If we do not provide a date
        # for the observation, then it will assume that it's the vernal equinox (sun at (0,0) in
        # ecliptic coordinates) in 2025.
        CCD_cent_pos = wcs.toWorld(sky_image.true_center)
        sky_level = euclidlike.getSkyLevel(bandpass, world_pos=CCD_cent_pos)

        # Currently the stray light fraction is not defined in the Euclid-like module
        # sky_level *= (1.0 + euclidlike.stray_light_fraction)

        # Note that makeSkyImage() takes a bit of time. If you do not care about the variable pixel
        # scale, you could simply compute an approximate sky level in e-/pix by multiplying
        # sky_level by euclidlike.pixel_scale**2, and add that to sky_image.
        wcs.makeSkyImage(sky_image, sky_level)

        # The other background is the expected thermal backgrounds in this band.
        # These are provided in e-/pix/s, so we have to multiply by the exposure time.
        # The thermal background is not important for VIS, only for NISP and is
        # currently not defined in the Euclid-like module.
        # sky_image += euclidlike.thermal_backgrounds[filter_name]*euclidlike.long_exptime

        # Draw the galaxies and stars into the image.
        # We want (most of) the object properties to be the same for all the filters.
        # E.g. the position, orintation, etc. should match up for all the observations.
        # To make this happen, we start an rng from the same seed each time.
        for i_obj in range(nobj):
            logger.info('Drawing image for object {} in band {}'.format(i_obj, filter_name))

            # The rng for object parameters should be the same for each filter to make sure
            # we get the same parameters, position, rotation in each color.
            # Note that this one follows the explicit sequence given in the config file,
            # starting with 12345, not BaseDeviate(12345).raw().
            obj_rng = galsim.UniformDeviate(seed + 1 + 10**6 + i_obj)
            # The rng for photon shooting should be different for each filter.
            phot_rng = galsim.UniformDeviate(seed1 + 1 + i_obj + ifilter*nobj)

            # We'll deal with this below.  The config processing calculates this before the
            # position, so to get the same answers, we need to do the same here.
            p = obj_rng()

            # Pick a random position in the image to draw it.
            # If we had a real galaxy catalog with positions in terms of RA, Dec we could use
            # wcs.toImage() to find where those objects should be in terms of (x, y).
            # Note that we could use wcs.toWorld() to get the (RA, Dec) for these (x, y) positions
            # if we wanted that information, but we don't need it.
            x = obj_rng() * euclidlike.n_pix_row
            y = obj_rng() * euclidlike.n_pix_col
            image_pos = galsim.PositionD(x, y)
            logger.debug('Position = %s', image_pos)

            # Now decide which of our three kinds of objects we want to draw:
            # 80% faint galaxy
            # 10% star
            # 10% bright galaxy
            if p < 0.8:
                # Faint galaxy
                logger.debug('Faint galaxy')

                # Select a random galaxy from the catalog.
                obj = cat1.makeGalaxy(chromatic=True, gal_type='parametric', rng=obj_rng)
                logger.debug('galaxy index = %s', obj.index)

                # Rotate the galaxy randomly
                theta = obj_rng() * 2 * np.pi * galsim.radians
                logger.debug('theta = %s', theta)
                obj = obj.rotate(theta)

            elif p < 0.9:
                # Star
                logger.debug('Star')

                # Use a log-normal distribution for the stellar fluxes.  This gives a few very
                # bright stars, but not too many.
                # cf. https://en.wikipedia.org/wiki/Log-normal_distribution
                # mu_x, sigma_x are the target mean, std.dev.
                # mu,sigma are the parameters of the function.
                mu_x = 1.e5
                sigma_x = 2.e5
                mu = np.log(mu_x**2 / (mu_x**2+sigma_x**2)**0.5)
                sigma = (np.log(1 + sigma_x**2/mu_x**2))**0.5
                gd = galsim.GaussianDeviate(obj_rng, mean=mu, sigma=sigma)
                flux = np.exp(gd())
                logger.debug('flux = %s', flux)

                # Normalize the SED to have this flux in the VIS band.
                sed = vega_sed.withFlux(flux, bandpass)

                obj = galsim.DeltaFunction() * sed

            else:
                # Bright galaxy
                logger.debug('Bright galaxy')

                obj = cat2.makeGalaxy(chromatic=True, gal_type='parametric', rng=obj_rng)
                logger.debug('galaxy index = %s', obj.index)

                # Scale up the area by a factor of 2, and the flux by a factor of 4.
                # This is not necessarily physical, but it is intended to add some more big,
                # bright galaxies to the scene to make the final image a bit more interesting.
                obj = obj.dilate(2) * 4

                # Rotate the galaxy randomly
                theta = obj_rng() * 2 * np.pi * galsim.radians
                logger.debug('theta = %s',theta)
                obj = obj.rotate(theta)

            # Convolve the (chromatic) object with the (chromatic) PSF.
            final = galsim.Convolve(obj, psf)
            stamp = final.drawImage(bandpass, center=image_pos, wcs=wcs.local(image_pos),
                                    method='phot', rng=phot_rng)

            # Find the overlapping bounds between the large image and the individual stamp.
            bounds = stamp.bounds & full_image.bounds

            # Add this to the corresponding location in the large image.
            full_image[bounds] += stamp[bounds]

        # Now we're done with the per-object drawing for this image.  The rest will be done for the
        # entire image at once.
        logger.info('All objects have been drawn for filter %s.', filter_name)
        logger.info('Adding the noise and detector non-idealities.')

        # At this point in the image generation process, an integer number of photons gets
        # detected.  Because of how GalSim's photon shooting works for InterpolatedImage
        # (used implicitly in the PSF implementation), the image has non-integral values at this
        # point.  So the first thing we do is quantize that to an integer number of photons.
        full_image.quantize()

        # Add the sky image.  Note: the galaxies already have Poisson noise because we are photon
        # shooting, but the sky image does not.  We want to preserve the expectation value of the
        # sky image (to subtract it off below), so we need a copy, which we can add noise to.
        poisson_noise = galsim.PoissonNoise(image_rng)
        sky_image_realized = sky_image.copy()
        sky_image_realized.addNoise(poisson_noise)
        full_image += sky_image_realized

        # Now that all sources of signal (from astronomical objects and background) have been added
        # to the image, we can start adding noise and detector effects. 
        # Here we step through the process and explain these in a bit
        # more detail without using that utility.

        # The subsequent steps account for the non-ideality of the detectors.

        # 1) Adding dark current to the image:
        # Even when the detector is unexposed to any radiation, the electron-hole pairs that
        # are generated within the depletion region due to finite temperature are swept by the
        # high electric field at the junction of the photodiode. This small reverse bias
        # leakage current is referred to as 'dark current'. It is specified by the average
        # number of electrons reaching the detectors per unit time and has an associated
        # Poisson noise since it is a random event.
        # dark_current = euclidlike.dark_current*euclidlike.exptime
        # dark_noise = galsim.DeviateNoise(galsim.PoissonDeviate(image_rng, dark_current))
        # full_image.addNoise(dark_noise)
        # sky_image += dark_current # (also want to subtract this expectation value along with sky)

        # NOTE: Sky level and dark current might appear like a constant background that can be
        # simply subtracted. However, these contribute to the shot noise and matter for the
        # non-linear effects that follow. Hence, these must be included at this stage of the
        # image generation process. We subtract these backgrounds in the end.

        # 2) Applying a quadratic non-linearity:
        # In order to convert the units from electrons to ADU, we must use the gain factor. The gain
        # has a weak dependency on the charge present in each pixel. This dependency is accounted
        # for by changing the pixel values (in electrons) and applying a constant nominal gain
        # later, which is unity in our demo.

        # Apply the Euclid nonlinearity routine, which knows all about the nonlinearity expected in
        # the Euclid detectors.
        # This routine is currently not defined in the Euclid-like module.
        # euclidlike.applyNonlinearity(full_image)

        # Note that users who wish to apply some other nonlinearity function (perhaps for other NISP
        # detectors, or for CCDs) can use the more general nonlinearity routine, which uses the
        # following syntax:
        # full_image.applyNonlinearity(NLfunc=NLfunc)
        # with NLfunc being a callable function that specifies how the output image pixel values
        # should relate to the input ones.
        # logger.debug('Applied nonlinearity to {0}-band image'.format(filter_name))

        # 3) Adding read noise:
        # Read noise is the noise due to the on-chip amplifier that converts the charge into an
        # analog voltage.  We already applied the Poisson noise due to the sky level, so read noise
        # should just be added as Gaussian noise:
        read_noise = galsim.GaussianNoise(image_rng, sigma=euclidlike.read_noise)
        full_image.addNoise(read_noise)
        logger.debug('Added readnoise to {0}-band image'.format(filter_name))

        # We divide by the gain to convert from e- to ADU. For now, there is just a single number
        # for each CCD, which is equal to 3.4
        full_image /= euclidlike.gain
        sky_image /= euclidlike.gain

        # Finally, the analog-to-digital converter reads in an integer value.
        full_image.quantize()
        # sky_image.quantize() # Quantizing the background will actually lead to some issues.
        # Note that the image type after this step is still a float.  If we want to actually
        # get integer values, we can do new_img = galsim.Image(full_image, dtype=int)

        # Here we add quantization noise to the image. This prevent to have problems due
        # to the rounding. For example, given the very low amplitude of the sky background,
        # we can have spatial variations of only ~1 ADU. This will be impossible to 
        # properly pick up by tools like SExtractor. Adding this noise prevents this problem
        # and does not change the signal in the image.
        quantization_noise = full_image.copy()
        quantization_noise.fill(0)
        quantization_noise.addNoise(galsim.DeviateNoise(galsim.UniformDeviate(image_rng)))
        full_image += quantization_noise
        full_image -= 0.5

        # Since many people are used to viewing background-subtracted images, we provide a
        # version with the background subtracted (also rounding that to an int).
        full_image -= sky_image

        logger.debug('Subtracted background for {0}-band image'.format(filter_name))
        # Write the final image to a file.
        out_filename = os.path.join(outpath, 'end_to_end_demo_{0}.fits'.format(filter_name))
        full_image.write(out_filename)

        logger.info('Completed {0}-band image.'.format(filter_name))


if __name__ == "__main__":
    main(sys.argv[1:])
