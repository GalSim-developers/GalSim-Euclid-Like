## Reference

* `Euclid_VIS.vis.dat` file contains the Euclid VIS bandpass found [here](http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?mode=browse&gname=Euclid&gname2=VIS&asttype=)
* `NISP-PHOTO-PASSBANDS-V1-*_throughput.dat` (3 files) contain the NISP bandpasses from [here](https://euclid.esac.esa.int/msp/refdata/nisp/NISP-PHOTO-PASSBANDS-V1)
* `monopsfs_6_6.fits.gz` is the oversampled Euclid PSF at 17 different wavelengths given in `psf_wavelengths.dat`
    - The minimum wavelength is 540 nm and the maximum is 910 nm
    - The pixel size is 0.1 arcsec
    - The image size is 480 x 480
* `pv_coeffs.dat` contains the CD and PV coefficient to build Euclid-like WCS. At the moment they are set to 0 (no distortions).
* `ccd_data.dat` contains the CCD shifts with respect to the center of the focal plane. Derived from [Scaramella et al.](https://arxiv.org/abs/2108.01201)
* `config.zip` contains the configuration file for the Euclid PSF tool used to produce the stored PSF images.
