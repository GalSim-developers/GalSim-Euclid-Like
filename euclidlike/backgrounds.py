"""
This file includes any routines needed to define the background level, for which the main contribution (currently the only one implemented) is zodiacal light.
"""
import numpy as np
import galsim

def getSkyLevel(bandpass, world_pos=None, exptime=None, epoch=2025, date=None):
    """
    Get the expected sky level for a Euclid observation due to zodiacal light for this bandpass
    and position.

    This routine can take an arbitray galsim.Bandpass() and calculate the zodiacal background
    directly.  This is approximately 3x slower than using lookup tables.

    The numbers that are returned are in units of e-/arcsec^2.  The result can either be multiplied
    by the approximate pixel area to get e-/pix, or the result can be used with wcs.makeSkyImage()
    to make an image of the sky that properly includes the actual pixel area as a function of
    position on the detector.

    Note that in general results will depend on the adopted model for zodiacal light, and these are
    uncertain at the ~10% level.

    Positions should be specified with the ``world_pos`` keyword, which must be a CelestialCoord
    object.  If no ``world_pos`` is supplied, then the routine will use a default position that
    looks sensibly away from the sun.

    Parameters:
        bandpass:   A Bandpass object.
        world_pos:  A position, given as a CelestialCoord object.  If None, then the routine
                    will use an ecliptic longitude of 90 degrees with respect to the sun
                    position (as a fair compromise between 0 and 180), and an ecliptic latitude
                    of 30 degrees with respect to the sun position (decently out of the plane
                    of the Earth-sun orbit). [default: None]
        exptime:    Exposure time in seconds.  If None, use the longer Euclid exposure time.
                    [default: None]
        epoch:      The epoch to be used for estimating the obliquity of the ecliptic when
                    converting ``world_pos`` to ecliptic coordinates.  This keyword is only used
                    if ``date`` is None, otherwise ``date`` is used to determine the ``epoch``.
                    [default: 2025]
        date:       The date of the observation, provided as a python datetime object.  If None,
                    then the conversion to ecliptic coordinates assumes the sun is at ecliptic
                    coordinates of (0,0), as it is at the vernal equinox. [default: None]

    Returns:
        the expected sky level in e-/arcsec^2.

    """
    if exptime is None:
        from . import long_exptime as exptime

    # Check for proper type for position, and extract the ecliptic coordinates.
    if world_pos is None:
        # Use our defaults for the case of unspecified position.
        ecliptic_lat = 30.0 # degrees
        ecliptic_lon = 90.0 # degrees
    else:
        if not isinstance(world_pos, galsim.CelestialCoord):
            raise TypeError("world_pos must be supplied as a galsim.CelestialCoord.")
        if date is not None:
            epoch = date.year
        ecliptic_lon, ecliptic_lat = world_pos.ecliptic(epoch=epoch, date=date)
        ecliptic_lon = ecliptic_lon.deg
        ecliptic_lat = ecliptic_lat.deg

    if not isinstance(bandpass, galsim.Bandpass):
        raise TypeError("bandpass must be supplied as a galsim.Bandpass")
    
    # Then call our routine to calculate the sky background
    # For now, simply require that the blue/red limits be in nm, as will be the case in general for
    # how we've set up the bandpasses
    if bandpass.wave_type != "nm":
        raise ValueError("Error, bandpass wavelength units are not as assumed by sky level routine")
    sky_val = getZodiBackground(ecliptic_lat, ecliptic_lon, bandpass.blue_limit/1000.,
                                bandpass.red_limit/1000., bandpass.wave_list/1000.,
                                bandpass.func(bandpass.wave_list))
    
    # Now, convert to the right units, and return.  (See docstring for explanation.)
    # Multiply by exposure time.
    from . import gain
    from . import collecting_area
    sky_val *= exptime
    sky_val *= gain
    sky_val *= collecting_area/1e4

    # The result is now the sky level in e-/arcsec^2.
    return sky_val  

def getZodiBackground(ecl_lat, ecl_dlon, lambda_min, lambda_max, Tlambda, T):
    """
    This helper routine is used with permission from Chris Hirata's Exposure Time Calculator and enables the
    calculation of the zodiacal light in photons/m^2/arcsec^2/sec. The ETC may be found here:

    http://www.tapir.caltech.edu/~chirata/web/software/space-etc/

    The code was ported from C to python by Michael Troxel (matroxel on GitHub).  The units are
    exactly as in the original C code, and we convert to any other needed units outside of this
    routine, in user-facing code. 

    Parameters:
        ecl_lat:      Ecliptic latitude (degrees)
        ecl_dlon:     Ecliptic longitude (degrees)
        lambda_min:   Minimum wavelength for the bandpass (microns)
        lambda_max:   Maximum wavelength for the bandpass (microns)
        Tlambda:      Numpy array containing a grid of wavelength values for the bandpass (microns)
        T:            Numpy array containing the throughput (normalized to be between 0-1) for the
                      bandpass
    
    Returns:
        A floating point value for the zodiacal light in photons/m^2/arcsec^2/sec
    """
    deg = np.pi/180. # degrees 
    Nlambda = 100 # number of integration points in wavelength 
    # Sky brightness table: rows are varying ecliptic latitude, cols are varying longitude,
    # at the values shown.
    # 
    # This is at 0.5um in units of 1e-8 W/m^2/sr/um.
    # Electronic version of Table 17 of Leinert (1997), except for placeholders (1's) at
    # elongation <15 degrees (where we should not use this routine anyway!).
    nlat = 11
    nlon = 19
    betaTable = [0,5,10,15,20,25,30,45,60,75,90]
    dlonTable = [0,5,10,15,20,25,30,35,40,45,60,75,90,105,120,135,150,165,180]
    skyTable = [
        1,    1,    1, 3140, 1610, 985, 640, 275, 150, 100, 77,
        1,    1,    1, 2940, 1540, 945, 625, 271, 150, 100, 77,
        1,    1, 4740, 2470, 1370, 865, 590, 264, 148, 100, 77,
        11500, 6780, 3440, 1860, 1110, 755, 525, 251, 146, 100, 77,
        6400, 4480, 2410, 1410,  910, 635, 454, 237, 141,  99, 77,
        3840, 2830, 1730, 1100,  749, 545, 410, 223, 136,  97, 77,
        2480, 1870, 1220,  845,  615, 467, 365, 207, 131,  95, 77,
        1650, 1270,  910,  680,  510, 397, 320, 193, 125,  93, 77,
        1180,  940,  700,  530,  416, 338, 282, 179, 120,  92, 77,
        910,  730,  555,  442,  356, 292, 250, 166, 116,  90, 77,
        505,  442,  352,  292,  243, 209, 183, 134, 104,  86, 77,
        338,  317,  269,  227,  196, 172, 151, 116,  93,  82, 77,
        259,  251,  225,  193,  166, 147, 132, 104,  86,  79, 77,
        212,  210,  197,  170,  150, 133, 119,  96,  82,  77, 77,
        188,  186,  177,  154,  138, 125, 113,  90,  77,  74, 77,
        179,  178,  166,  147,  134, 122, 110,  90,  77,  73, 77,
        179,  178,  165,  148,  137, 127, 116,  96,  79,  72, 77,
        196,  192,  179,  165,  151, 141, 131, 104,  82,  72, 77,
        230,  212,  195,  178,  163, 148, 134, 105,  83,  72, 77
    ]
    # Solar spectrum: in units of W/m^2/sr/um at log10(lambda/um) = -0.80(0.01)+0.40
    # Ref: Colina, Bohlin, Castelli 1996 AJ 112, 307
    # 
    # V band (550 nm) is SolarSpec[54]
    SolarSpec = [
        1.87138e-01, 2.61360e-01, 4.08020e-01, 6.22197e-01, 9.02552e-01, 1.51036e+00, 2.25890e+00, 2.75901e+00, 4.03384e+00, 5.42817e+00, 
        7.26182e+00, 1.01910e+01, 2.01114e+01, 3.62121e+01, 4.31893e+01, 5.43904e+01, 4.91581e+01, 4.95091e+01, 4.95980e+01, 5.93722e+01, 
        5.27380e+01, 1.02502e+02, 1.62682e+02, 2.53618e+02, 2.01084e+02, 2.08273e+02, 4.05163e+02, 5.39830e+02, 5.31917e+02, 6.31200e+02, 
        7.06134e+02, 8.13653e+02, 1.00508e+03, 9.56536e+02, 9.50568e+02, 9.82400e+02, 1.06093e+03, 1.12669e+03, 1.09922e+03, 1.10224e+03, 
        1.36831e+03, 1.72189e+03, 1.74884e+03, 1.59871e+03, 1.74414e+03, 1.98823e+03, 2.02743e+03, 2.00367e+03, 2.03584e+03, 1.90296e+03, 
        1.93097e+03, 1.86594e+03, 1.86655e+03, 1.87957e+03, 1.87978e+03, 1.83915e+03, 1.84447e+03, 1.80371e+03, 1.76779e+03, 1.70796e+03, 
        1.66589e+03, 1.61456e+03, 1.53581e+03, 1.51269e+03, 1.44957e+03, 1.39215e+03, 1.34031e+03, 1.28981e+03, 1.24501e+03, 1.19548e+03, 
        1.15483e+03, 1.10546e+03, 1.06171e+03, 9.94579e+02, 9.54006e+02, 9.15287e+02, 8.63891e+02, 8.31183e+02, 7.95761e+02, 7.62568e+02, 
        7.27589e+02, 6.94643e+02, 6.60883e+02, 6.21830e+02, 5.83846e+02, 5.59624e+02, 5.34124e+02, 5.06171e+02, 4.80985e+02, 4.63139e+02, 
        4.39482e+02, 4.13122e+02, 3.94543e+02, 3.75591e+02, 3.56069e+02, 3.35294e+02, 3.16374e+02, 2.98712e+02, 2.82737e+02, 2.69581e+02, 
        2.49433e+02, 2.36936e+02, 2.21403e+02, 2.04770e+02, 1.87379e+02, 1.75880e+02, 1.60408e+02, 1.46210e+02, 1.36438e+02, 1.24412e+02, 
        1.16500e+02, 1.07324e+02, 9.89669e+01, 9.12134e+01, 8.28880e+01, 7.71064e+01, 7.06245e+01, 6.42367e+01, 5.87697e+01, 5.39387e+01, 
        4.98208e+01
    ]
    # Put longitude between 0 and 180
    ecl_dlon = np.abs(ecl_dlon)
    ecl_dlon -= 360*np.floor(ecl_dlon/360)
    if (ecl_dlon>180): ecl_dlon = 360-ecl_dlon
    if (ecl_dlon>180): ecl_dlon = 180
    # Set latitude to be positive
    ecl_lat = np.abs(ecl_lat)
    if (ecl_lat>90): ecl_lat = 90
    # Compute elongation (Sun angle). Complain if <15 degrees.
    z = np.cos(ecl_lat*deg)*np.cos(ecl_dlon*deg)
    if z>=1: elon = 0
    elif z<=-1: elon = 180
    else: elon = np.arccos(z)/deg
    if (elon<15):
        print("Error: get_zodi_bkgnd: elongation = "+str(elon)+" degrees out of valid range.\n");
        return -10
    # Compute sky brightness at 0.5um in units of 1e-8 W/m^2/sr/um.
    # Fit to Table 17 of Leinert (1997).
    ilat=0
    while betaTable[ilat+1]<ecl_lat and ilat<nlat-2: ilat+=1
    ilon=0
    while dlonTable[ilon+1]<ecl_dlon and ilon<nlon-2: ilon+=1
    frlat = (ecl_lat-betaTable[ilat])/(betaTable[ilat+1]-betaTable[ilat])
    frlon = (ecl_dlon-dlonTable[ilon])/(dlonTable[ilon+1]-dlonTable[ilon])
    sky05 = np.exp(
        np.log(skyTable[ilat  +(ilon  )*nlat]) * (1.-frlat) * (1.-frlon)
        +np.log(skyTable[ilat  +(ilon+1)*nlat]) * (1.-frlat) * (   frlon)
        +np.log(skyTable[ilat+1+(ilon  )*nlat]) * (   frlat) * (1.-frlon)
        +np.log(skyTable[ilat+1+(ilon+1)*nlat]) * (   frlat) * (   frlon)
    )
    # Integrate over wavelengths
    zodi_tot = 0.
    dlambda = (lambda_max-lambda_min)/float(Nlambda)
    for ilambda in range(Nlambda):
        lambda_ = lambda_min + (ilambda+0.5)/Nlambda * (lambda_max-lambda_min)
        # Solar spectrum at this wavelength: F_lambda/F_{0.5um}
        index_lambda = 100*np.log(lambda_)/np.log(10.) + 80
        ilam = int(np.floor(index_lambda))
        frlam = index_lambda - ilam
        sun_spec = (SolarSpec[ilam] + frlam*(SolarSpec[ilam+1]-SolarSpec[ilam]))/SolarSpec[50]

        # Color correction relative to solar
        if (lambda_>0.5):
            if elon>90: fco = 1. + (0.6)*np.log(lambda_/0.5)/np.log(10.)
            elif elon<30: fco = 1. + (0.8)*np.log(lambda_/0.5)/np.log(10.)
            else: fco = 1. + (0.8-0.2*(elon-30)/60.)*np.log(lambda_/0.5)/np.log(10.)
        else:
            if elon>90: fco = 1. + (0.9)*np.log(lambda_/0.5)/np.log(10.)
            elif elon<30: fco = 1. + (1.2)*np.log(lambda_/0.5)/np.log(10.)
            else: fco = 1. + (1.2-0.3*(elon-30)/60.)*np.log(lambda_/0.5)/np.log(10.)

        # The integral for the zodiacal foreground.
        # Here sky05*fco*sun_spec*dlambda is the power per unit area per unit solid angle in this
        # wavelength range (Units: 1e-8 W/m^2/sr).
        # 
        # The conversion from 1e-8 W --> photons/sec is 5.03411747e10*lambda(um).
        zodi_tot += sky05 * fco * sun_spec * dlambda * \
            5.03411747e10*lambda_* np.interp(lambda_,Tlambda,T)

    # We now have the zodi level in photons/m^2/sr/sec. Convert to photons/m^2/arcsec^2/sec.
    return (zodi_tot / 4.2545170296152206e10)

