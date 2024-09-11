"""
@file euclidlike_wcs.py

This file includes any routines needed to define and use the Euclid-like WCS.
Current version is based on the focal plane description detailed in 
Scaramella et al. (Fig. 2 and Table 1).
The distortion coefficients were derived from the ERO release.

Scaramella et al.: https://arxiv.org/abs/2108.01201
"""
import numpy as np
import os
import coord
import datetime

from importlib.resources import files

from galsim import GSFitsWCS, FitsHeader
from galsim import PositionD
from galsim import BoundsI
from galsim import GalSimRangeError, GalSimError

from .instrument_params import (
    pixel_scale_mm,
    plate_scale,
    n_ccd,
    n_ccd_row,
    n_ccd_col,
    n_pix_row,
    n_pix_col,
    det2ccd,
    ccd2det,
    min_sun_angle,
    max_sun_angle,
)


# Basic Euclid reference info.
n_pv = 10  # Number of PV coefficients used, where arrays are n_pv x 2 in dimension

# Version-related information. At the moment the version are based on Scaramella et al. 2021
tel_name = "Euclid"
instr_name = "VIS"
optics_design_ver = "S21"

# Information about center points of the CCDs in the VIS focal plane coordinate system
# coordinates.
# Units are mm.
infile = files('euclidlike.data').joinpath('ccd_data.dat')
ccd_data = np.genfromtxt(
    infile,
    dtype=[
        ("CCD_ID", "<U3"),
        ("xc", np.float64),
        ("yc", np.float64),
        ("uc", np.float64),
        ("vc", np.float64)
    ]
)
ccd_xc_mm = ccd_data["xc"]
ccd_yc_mm = ccd_data["yc"]
ccd_crval_u_deg = ccd_data["uc"]
ccd_crval_v_deg = ccd_data["vc"]
# Nominal center of FPA from the payload axis in this coordinate system, in mm and as an angle
# (neglecting distortions - to be included later).
fpa_xc_mm = 0.0
fpa_yc_mm = 0.859/plate_scale
xc_fpa = 0.*coord.radians
yc_fpa = np.deg2rad(0.859)*coord.radians

# Gaps between CCDs
x_gap_as = 12.7
y_gap_as = 64.4
x_gap_mm = x_gap_as/plate_scale
y_gap_mm = y_gap_as/plate_scale

# The next array contains rotation offsets of individual CCD Y axis relative to FPA.
ccd_rot = np.zeros_like(ccd_xc_mm)

# Rotation of WFI local axes relative to payload axes: this is expressed as a CCW rotation
# relative to observatory +Z direction.
theta_fpa = 0.*coord.degrees

# File with PV coefficients.
pv_filename = files('euclidlike.data').joinpath('pv_coeffs.dat')


def getWCS(world_pos, PA=None, date=None, CCDs=None, PA_is_FPA=False, SAA=None):
    """
    This routine returns a dict containing a WCS for each of the Euclid CCDs.
    The Euclid CCDs are labeled 0-35, so these numbers are used as the keys in
    the dict.  Alternatively the user can request a subset of the CCDs using
    the ``CCDs`` option.  The basic instrument parameters used to create the
    WCS correspond to those in Early Release Observations (ERO) data and are 
    likely to change in the future.

    The user must specify a position for observation, at which the center of
    the focal plane array will point. This must be supplied as a CelestialCoord
    ``world_pos``.  In general, only certain positions are observable on
    certain dates, and for a given position there is an optimal position angle
    for the observatory (with the solar panels pointed as directly towards the
    sun as possible).  Users who are knowledgable about these details may
    choose to supply a position angle as ``PA``, either for the observatory or
    for the focal plane (using ``PA_is_FPA`` to indicate this).  But otherwise,
    the routine will simply choose the optimal position angle for a given date.

    Parameters:
        world_pos:      A `galsim.CelestialCoord` indicating the position to
                        observe at the center of the focal plane array (FPA).
                        Note that if the given position is not observable on
                        the given date, then the routine will raise an
                        exception.
        PA:             A `galsim.Angle` representing the position angle of the
                        observatory +Y axis, unless ``PA_is_FPA=True``, in
                        which case it's the position angle of the FPA.  For
                        users to do not care about this, then leaving this as
                        None will result in the routine using the supplied
                        ``date`` and ``world_pos`` to select the optimal
                        orientation for the observatory.  Note that if a user
                        supplies a ``PA`` value, the routine does not check
                        whether this orientation is actually allowed.
                        [default: None]
        date:           The date of the observation, as a python datetime
                        object.  If None, then the vernal equinox in 2025 will
                        be used.  [default: None]
        PA_is_FPA:      If True, then the position angle that was provided was
                        the PA of the focal plane array, not the observatory.
                        [default: False]
        CCDs:           A single number or iterable giving the CCDs for which
                        the WCS should be obtained.  If None, then the WCS is
                        calculated for all CCDs. [default: None]
        SAA:            A `galsim.Angle` representing the Solar Aspect Angle
                        of the telescope for the observation. If not provided,
                        it will be computed internally.

    Returns:
        A dict of WCS objects for each CCD.
    """
    # First just parse the input quantities.
    date, CCDs, pa_fpa, pa_obsy = _parse_WCS_inputs(
        world_pos, PA, date, PA_is_FPA, CCDs, SAA,
    )

    # Note, this routine reads in the coeffs.  We don't use them until later, but read them in for
    # all CCDs at once.
    PVs = _parse_pv_file(pv_filename)

    # Loop over CCDs:
    wcs_dict = {}
    for i_ccd in CCDs:
        # Set up the header.
        header = []
        # Populate some necessary variables in the FITS header that are always the same, regardless of
        # input and CCD number.
        _populate_required_fields(header)

        # And populate some things that just depend on the overall locations or other input, not on
        # the CCD.
        header.extend([
            ("RA",  world_pos.ra / coord.degrees, "reconstructed pointing RA"),
            ("DEC", world_pos.dec / coord.degrees, "reconstructed pointing DEC"),
            ("PA", pa_fpa / coord.degrees, "reconstructed position angle"),
            ('RA_COMM', world_pos.ra / coord.degrees, "commanding pointing RA"),
            ('DEC_COMM', world_pos.dec / coord.degrees, "commanding pointing DEC"),
            # ('PA_OBSY', pa_obsy / coord.degrees, "position angle of observatory Y axis (deg)"),
            ('PA_COMM', pa_fpa / coord.degrees, "commanding position angle"),
            ('DETID', i_ccd, "CCD-ID field from science TM (0-35)"),
            ('CCDID', det2ccd[i_ccd], "e.g. Detector ID, e.g. '0-0', '1-1' ... '6-6'"),
        ])

        # Get position of CCD center given the center of the FPA and the
        # orientation angle of the focal plane.
        crval, u, v = _get_ccd_center_pos(i_ccd, world_pos, pa_fpa)

        # Compute the position angle of the local pixel Y axis.
        # This requires projecting local North onto the detector axes.
        # Start by adding any CCD-unique rotation relative to FPA axes:s
        sca_tp_rot = pa_fpa + ccd_rot[i_ccd]*coord.degrees

        # Go some reasonable distance from crval in the +y direction.  Say, 1 degree.
        plus_y = world_pos.deproject(u, v + 1*coord.degrees, projection='gnomonic')
        # Find the angle between this point, crval and due north.
        north = coord.CelestialCoord(0.*coord.degrees, 90.*coord.degrees)
        pa_sca = sca_tp_rot - crval.angleBetween(plus_y, north)

        # Rotate by pa_fpa.
        cos_pa_sca = np.cos(pa_sca)
        sin_pa_sca = np.sin(pa_sca)

        header.extend([
            ('CRVAL1', crval.ra / coord.degrees),
            ('CRVAL2', crval.dec / coord.degrees),
            ('CD1_1', cos_pa_sca * PVs["CD1_1"][i_ccd] - sin_pa_sca * PVs["CD1_2"][i_ccd]),
            ('CD1_2', sin_pa_sca * PVs["CD1_1"][i_ccd] + cos_pa_sca * PVs["CD1_2"][i_ccd]),
            ('CD2_1', cos_pa_sca * PVs["CD2_1"][i_ccd] - sin_pa_sca * PVs["CD2_2"][i_ccd]),
            ('CD2_2', sin_pa_sca * PVs["CD2_1"][i_ccd] + cos_pa_sca * PVs["CD2_2"][i_ccd]),
        ])
        for key in PVs.dtype.names:
            if "PV" not in key:
                continue
            header.append((key, PVs[key][i_ccd]))

        header = FitsHeader(header)
        wcs = GSFitsWCS(header=header)
        # Store the original header as an attribute of the WCS. This ensures that we have all the
        # extra keywords for whenever an image with this WCS is written to file.
        wcs.header = header
        wcs_dict[i_ccd] = wcs

    return wcs_dict


def convertCenter(world_pos, CCD, PA=None, date=None, SAA=None, PA_is_FPA=False, tol=0.5*coord.arcsec):
    """
    This is a simple helper routine that takes an input position ``world_pos`` that is meant to
    correspond to the position of the center of an CCD, and tells where the center of the focal
    plane array should be.  The goal is to provide a position that can be used as an input to
    getWCS(), which wants the center of the focal plane array.

    The results of the calculation are deterministic if given a fixed position angle (PA).  If it's
    not given one, it will try to determine the best one for this location and date, like getWCS()
    does.

    Because of distortions varying across the focal plane, this routine has to iteratively correct
    its initial result based on empirical tests.  The ``tol`` kwarg can be used to adjust how
    careful it will be, but it always does at least one iteration.

    Parameters:
        world_pos:  A galsim.CelestialCoord indicating the position to observe at the center of the
                    given CCD.  Note that if the given position is not observable on
                    the given date, then the routine will raise an exception.
        CCD:        A single number giving the CCD for which the center should be located at
                    ``world_pos``.
        PA:         galsim.Angle representing the position angle of the observatory +Y axis, unless
                    ``PA_is_FPA=True``, in which case it's the position angle of the FPA.  For
                    users to do not care about this, then leaving this as None will result in the
                    routine using the supplied ``date`` and ``world_pos`` to select the optimal
                    orientation for the observatory.  Note that if a user supplies a ``PA`` value,
                    the routine does not check whether this orientation is actually allowed.
                    [default: None]
        SAA:        A `galsim.Angle` representing the Solar Aspect Angle
                    of the telescope for the observation. If not provided,
                    it will be computed internally.
        date:       The date of the observation, as a python datetime object.  If None, then the
                    vernal equinox in 2025 will be used.  [default: None]
        PA_is_FPA:  If True, then the position angle that was provided was the PA of the focal
                    plane array, not the observatory. [default: False]
        tol:        Tolerance for errors due to distortions, as a galsim.Angle.
                    [default: 0.5*galsim.arcsec]

    Returns:
        A CelestialCoord object indicating the center of the focal plane array.
    """
    if not isinstance(CCD, int):
        raise TypeError("Must pass in an int corresponding to the CCD")
    if not isinstance(tol, coord.Angle):
        raise TypeError("tol must be a galsim.Angle")
    use_CCD = CCD
    # Parse inputs appropriately.
    _, _, pa_fpa, _ = _parse_WCS_inputs(
        world_pos, PA, date, PA_is_FPA, [CCD], SAA,
    )

    # Now pretend world_pos was the FPA center and we want to find the location of this CCD:
    _, u, v = _get_ccd_center_pos(use_CCD, world_pos, pa_fpa)
    # The (u, v) values give an offset, and we can invert this.
    fpa_cent = world_pos.deproject(-u, -v, projection='gnomonic')
    # This is only approximately correct, especially for detectors that are far from the center of
    # the FPA, because of distortions etc.  We can do an iterative correction.
    # For the default value of 'tol', typically just 1-2 iterations are needed.
    shift_val = 1000.0  # arcsec
    while shift_val > tol/coord.arcsec:
        test_wcs = getWCS(fpa_cent, PA, date, use_CCD, PA_is_FPA)[use_CCD]
        im_cent_pos = PositionD(n_pix_row/2, n_pix_col/2)
        test_sca_pos = test_wcs.toWorld(im_cent_pos)
        delta_ra = np.cos(world_pos.dec)*(world_pos.ra-test_sca_pos.ra)
        delta_dec = world_pos.dec-test_sca_pos.dec
        shift_val = np.abs(world_pos.distanceTo(test_sca_pos)/coord.arcsec)
        fpa_cent = coord.CelestialCoord(fpa_cent.ra + delta_ra, fpa_cent.dec + delta_dec)

    return fpa_cent


def findCCD(wcs_dict, world_pos, include_border=False):
    """
    This is a subroutine to take a dict of WCS (one per CCD) from euclidlike.getWCS() and query
    which CCD a particular real-world coordinate would be located on.  The position (``world_pos``)
    should be specified as a galsim.CelestialCoord.  If the position is not located on any of the
    CCDs, the result will be None.  Note that if ``wcs_dict`` does not include all CCDs in it, then
    it's possible the position might lie on one of the CCDs that was not included.

    Depending on what the user wants to do with the results, they may wish to use the
    ``include_border`` keyword.  This keyword determines whether or not to include an additional
    border corresponding to half of the gaps between CCDs.  For example, if a user is drawing a
    single image they may wish to only know whether a given position falls onto a CCD, and if so,
    which one (ignoring everything in the gaps).  In contrast, a user who plans to make a sequence
    of dithered images might find it most useful to know whether the position is either on a CCD or
    close enough that in a small dither sequence it might appear on the CCD at some point.  Use of
    ``include_border`` switches between these scenarios.

    Parameters:
        wcs_dict:        The dict of WCS's output from euclidlike.getWCS().
        world_pos:       A galsim.CelestialCoord indicating the sky position of interest.
        include_border:  If True, then include the half-border around CCD to cover the gap
                         between each sensor. [default: False]

    Returns:
        an integer value of the CCD on which the position falls, or None if the position is not
        on any CCD.

    """
    # Sanity check args.
    if not isinstance(wcs_dict, dict):
        raise TypeError("wcs_dict should be a dict containing WCS output by euclidlike.getWCS.")

    if not isinstance(world_pos, coord.CelestialCoord):
        raise TypeError("Position on the sky must be given as a galsim.CelestialCoord.")

    # Set up the minimum and maximum pixel values, depending on whether or not to include the
    # border.  We put it immediately into a galsim.BoundsI(), since the routine returns xmin, xmax,
    # ymin, ymax:
    xmin, xmax, ymin, ymax = _calculate_minmax_pix(include_border)
    bounds_list = [
        BoundsI(x1, x2, y1, y2)
        for x1, x2, y1, y2 in zip(xmin, xmax, ymin, ymax)
    ]

    ccd = None
    for i_ccd in wcs_dict:
        wcs = wcs_dict[i_ccd]
        image_pos = wcs.toImage(world_pos)
        if bounds_list[i_ccd].includes(image_pos):
            ccd = i_ccd
            break

    return ccd


def _calculate_minmax_pix(include_border=False):
    """
    This is a helper routine to calculate the minimum and maximum pixel values that should be
    considered within a CCD, possibly including the complexities of including 1/2 of the gap
    between CCDs. In that case it depends on the detailed geometry of the Euclid focal plane.

    Parameters:
        include_border:     A boolean value that determines whether to include 1/2 of the gap
                            between CCDs as part of the CCD itself.  [default: False]

    Returns:
        a tuple of NumPy arrays for the minimum x pixel value, maximum x pixel value, minimum y
        pixel value, and maximum y pixel value for each CCD.
    """
    # First, set up the default (no border).
    # The minimum and maximum pixel values are (1, n_pix).
    min_x_pix = np.ones(n_ccd).astype(int)
    max_x_pix = min_x_pix + n_pix_col - 1
    min_y_pix = min_x_pix.copy()
    max_y_pix = min_y_pix + n_pix_row - 1

    # Then, calculate the half-gaps, grouping together CCDs whenever possible.
    if include_border:
        # The gaps are set following https://arxiv.org/abs/2108.01201 table 1
        # At the moment the ccd edges that are also at the edge of the focal
        # plane are handle separately. But, the border is like any other
        # border.
        half_border_x_pix = int(0.5*x_gap_mm/pixel_scale_mm)
        half_border_y_pix = int(0.5*y_gap_mm/pixel_scale_mm)
        val_x_edge = half_border_x_pix
        val_y_edge = half_border_y_pix

        for i in range(1, n_ccd_row+1):
            for j in range(1, n_ccd_col+1):
                i_ccd = ccd2det[f"{i}-{j}"]
                if i == 1:
                    min_x_pix[i_ccd] -= val_x_edge
                    max_x_pix[i_ccd] += half_border_x_pix
                elif i == 6:
                    min_x_pix[i_ccd] -= half_border_x_pix
                    max_x_pix[i_ccd] += val_x_edge
                else:
                    min_x_pix[i_ccd] -= half_border_x_pix
                    max_x_pix[i_ccd] += half_border_x_pix

                if j == 1:
                    min_y_pix[i_ccd] -= half_border_y_pix
                    max_y_pix[i_ccd] += val_y_edge
                elif j == 6:
                    min_y_pix[i_ccd] -= val_y_edge
                    max_y_pix[i_ccd] += half_border_y_pix
                else:
                    min_y_pix[i_ccd] -= half_border_y_pix
                    max_y_pix[i_ccd] += half_border_y_pix

    return min_x_pix, max_x_pix, min_y_pix, max_y_pix


def _populate_required_fields(header):
    """
    Utility routine to do populate some of the basic fields for the WCS headers for Euclid that
    don't require any interesting calculation.
    """
    header.extend([
        ('EQUINOX', 2000.0, "equinox of celestial coordinate system"),
        ('WCSAXES', 2, "number of World Coordinate System axes"),
        ('WCSNAME', 'viswcs_'+optics_design_ver),
        ('CRPIX1', n_pix_row/2, "x-coordinate of reference pixel"),
        ('CRPIX2', n_pix_col/2, "y-coordinate of reference pixel"),
        ('CTYPE1', "RA---TPV"),
        ('CTYPE2', "DEC--TPV"),
        ('SIMPLE', 'True'),
        ('BITPIX', 16),  # NOTE: BITPIX=8 in the real Euclid headers
        ('NAXIS', 0),
        ('EXTEND', 'True'),
        ('BZERO', 0),
        ('BSCALE', 1),
        ('TELESCOP', tel_name, "telescope used to acquire data"),
        ('INSTRUME', instr_name, "identifier for instrument used to acquire data"),
    ])


def _parse_pv_file(file):
    """
    Utility routine to parse the file with the CD and PV coefficients and hand back an array to be used
    for later calculations.
    """
    if not os.path.exists(file):
        raise OSError(
            "Cannot find file that should have Euclid PVs coefficients: %s" % file
        )

    # Parse the file, generated by make_sip_file.py in devel/roman directory.
    data = np.loadtxt(file, usecols=range(1, 39))

    PVs = np.zeros(
        n_ccd,
        dtype=[
            (f"CD{i}_{j}", np.float64)
            for i in range(1, 2+1)
            for j in range(1, 2+1)
        ] +
        [
            (f"PV{i}_{j}", np.float64)
            for i in range(1, 2+1)
            for j in range(1, n_pv+1)
        ]
    )
    for k in range(4):
        i = int(data[k, 0])
        j = int(data[k, 1])
        for i_ccd in range(n_ccd):
            PVs[f"CD{i}_{j}"][i_ccd] = data[k, 2 + i_ccd]
    for k in range(4, n_pv*2 + 4):
        i = int(data[k, 0])
        j = int(data[k, 1])
        for i_ccd in range(n_ccd):
            PVs[f"PV{i}_{j}"][i_ccd] = data[k, 2 + i_ccd]

    return PVs


def _get_ccd_center_pos(i_ccd, world_pos, pa_fpa):
    """
    This helper routine calculates the center position for a given CCD ``ccd`` given the position of
    the center of the focal plane array ``world_pos`` and an orientation angle for the observation.
    It is used by getWCS() and other routines.
    """
    # Go from the tangent plane position of the CCD center, to the actual celestial coordinate,
    # using `world_pos` as the center point of the tangent plane projection.  This celestial
    # coordinate for the CCD center is `crval`, which goes into the WCS as CRVAL1, CRVAL2.
    cos_pa = np.cos(pa_fpa)
    sin_pa = np.sin(pa_fpa)
    u = ccd_crval_v_deg[i_ccd] * cos_pa - ccd_crval_u_deg[i_ccd] * sin_pa
    v = ccd_crval_v_deg[i_ccd] * sin_pa + ccd_crval_u_deg[i_ccd] * cos_pa
    u = u * coord.degrees
    v = v * coord.degrees
    crval = world_pos.deproject(u, v, projection='gnomonic')
    return crval, u, v


def _parse_CCDs(CCDs):
    # This is a helper routine to parse the input CCDs (single number or iterable) and put it into a
    # convenient format.  It is used in roman_wcs.py.
    #
    # Check which CCDs are to be done.  Default is all (and they are 0-indexed).
    all_CCDs = np.arange(0, n_ccd, 1)
    # Later we will use the list of selected CCDs to decide which ones we're actually going to do
    # the calculations for.  For now, just check for invalid numbers.
    if CCDs is not None:
        # Make sure CCDs is iterable.
        if not hasattr(CCDs, '__iter__'):
            CCDs = [CCDs]
        # Then check for reasonable values.
        if min(CCDs) < 0 or max(CCDs) > n_ccd:
            raise GalSimRangeError("Invalid CCD.", CCDs, 0, n_ccd)
        # Check for uniqueness.  If not unique, make it unique.
        CCDs = list(set(CCDs))
    else:
        CCDs = all_CCDs
    return CCDs


def _parse_WCS_inputs(world_pos, PA, date, PA_is_FPA, CCDs, SAA=None):
    """
    This routine parses the various input options to getWCS() and returns what the routine needs to
    do its job.  The reason to pull this out is so other helper routines can use it.
    """

    # Parse input position
    if not isinstance(world_pos, coord.CelestialCoord):
        raise TypeError("Position on the sky must be given as a galsim.CelestialCoord!")

    # Get the date. (Vernal equinox in 2025, taken from
    # http://www.astropixels.com/ephemeris/soleq2001.html, if none was supplied.)
    if date is None:
        date = datetime.datetime(2025, 3, 20, 9, 2, 0)

    # Are we allowed to look here?
    if not allowedPos(world_pos, date, SAA):
        raise GalSimError("Error, Euclid cannot look at this position on this date!")

    # If position angle was not given, then get the optimal one:
    if PA is None:
        PA_is_FPA = False
        PA = bestPA(world_pos, date, SAA)
    else:
        # Just enforce type
        if not isinstance(PA, coord.Angle):
            raise TypeError("Position angle must be a galsim.Angle!")

    # Check which CCDs are to be done using a helper routine in the euclidlike module.
    CCDs = _parse_CCDs(CCDs)

    # Compute position angle of FPA f2 axis, where positive corresponds to the angle east of North.
    if PA_is_FPA:
        pa_fpa = PA
        pa_obsy = PA - theta_fpa
    else:
        pa_obsy = PA
        pa_fpa = PA + theta_fpa

    return date, CCDs, pa_fpa, pa_obsy


def allowedPos(world_pos, date, SAA=None):
    """
    This routine can be used to check whether Euclid would be allowed to look at a particular
    position (``world_pos``) on a given ``date``.   This is determined by the angle of this position
    relative to the Sun.

    In general, Euclid can point at angles relative to the Sun in the range 90+20 & 90-3 degrees.
    Obviously, pointing too close to the Sun would result in overly high sky backgrounds.  It is
    less obvious why Euclid cannot look at a spot directly opposite from the Sun (180 degrees on the
    sky).  The reason is that the observatory is aligned such that if the observer is looking at
    some sky position, the solar panels are oriented at 90 degrees from that position.  So it's
    always optimal for the observatory to be pointing at an angle of 90 degrees relative to the
    Sun.  It is also permitted to look within [-3, + 20] degrees of that optimal position.

    Parameters:
        world_pos:      A galsim.CelestialCoord indicating the position at which the observer
                        wishes to look.
        date:           A python datetime object indicating the desired date of observation.
        SAA:            A `galsim.Angle` representing the Solar Aspect Angle
                        of the telescope for the observation. If not provided,
                        it will be computed internally.

    Returns:
        True or False, indicating whether it is permitted to look at this position on this date.
    """
    # Find the Sun's location on the sky on this date.
    lam = coord.util.sun_position_ecliptic(date)
    sun = coord.CelestialCoord.from_ecliptic(lam, 0*coord.radians, date.year)

    # Find the angle between that and the supplied position
    if SAA is None:
        angle_deg = abs(world_pos.distanceTo(sun)/coord.degrees)
    else:
        angle_deg = SAA/coord.degrees

    # Check if it's within tolerance.
    min_ang = 90. - min_sun_angle/coord.degrees
    max_ang = 90. + max_sun_angle/coord.degrees
    return min_ang <= angle_deg <= max_ang


def bestPA(world_pos, date, SAA=None):
    """
    This routine determines the best position angle for the observatory for a given observation date
    and position on the sky.

    The best/optimal position angle is determined by the fact that the solar panels are at 90
    degrees to the position being observed, and it is best to have those facing the Sun as directly
    as possible.  Note that if a given ``world_pos`` is not actually observable on the given
    ``date``, then this routine will return None.

    Parameters:
        world_pos:      A galsim.CelestialCoord indicating the position at which the observer
                        wishes to look.
        date:           A python datetime object indicating the desired date of observation.
        SAA:            A `galsim.Angle` representing the Solar Aspect Angle
                        of the telescope for the observation. If not provided,
                        it will be computed internally.

    Returns:
        the best position angle for the observatory, as a galsim.Angle, or None if the position
        is not observable.
    """
    # First check for observability.
    if not allowedPos(world_pos, date, SAA):
        return None

    # Find the location of the sun on this date.  +X_observatory points out into the sky, towards
    # world_pos, while +Z is in the plane of the sky pointing towards the sun as much as possible.
    lam = coord.util.sun_position_ecliptic(date)
    sun = coord.CelestialCoord.from_ecliptic(lam, 0*coord.radians, date.year)
    # Now we do a projection onto the sky centered at world_pos to find the (u, v) for the Sun.
    sun_tp_x, sun_tp_y = world_pos.project(sun, 'gnomonic')

    # We want to rotate around by 90 degrees to find the +Y obs direction.  Specifically, we want
    # (+X, +Y, +Z)_obs to form a right-handed coordinate system.
    y_obs_tp_x, y_obs_tp_y = -sun_tp_y, sun_tp_x
    y_obs = world_pos.deproject(y_obs_tp_x, y_obs_tp_y, 'gnomonic')

    # Finally the observatory position angle is defined by the angle between +Y_observatory and the
    # celestial north pole.  It is defined as position angle east of north.
    north = coord.CelestialCoord(y_obs.ra, 90.*coord.degrees)
    obs_pa = world_pos.angleBetween(y_obs, north)
    return obs_pa
