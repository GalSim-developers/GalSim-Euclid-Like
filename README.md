# GalSim-Euclid-Like

Helper functions to generate simulations of Euclid-like images using GalSim.

This repository contains information about the Euclid space telescope and survey that is needed to
produce simulations using [GalSim](https://github.com/GalSim-developers/GalSim).  Some of the
information provided is approximate, aimed towards fast simulations rather than full accuracy in
representation of Euclid images.  Places where the information is only approximate are flagged and
described in the docstring, and we particularly highlight that the PSF is only approximate;
for details, see the docstring of the `getPSF()` method.  This library should enable generation of
Euclid-like images of sufficient fidelity for preliminary exploration of object detection,
photometry, deblending, and joint analysis with ground-based observatories.  For
applications requiring high precision such as weak lensing, the higher fidelity simulations
available within the Euclid Consortium should be used.

The repository includes two distinct packages:
* `euclidlike` has basic observatory, intrumentation, and survey information for Euclid.
This package can be used on its own along with GalSim to produce Euclid-like simulations.
* `euclidlike_imsim` has configuration scripts to produce large-scale Euclid-like simulation runs
based on the information in `euclidlike`.  It is based heavily on [`roman_imsim`](https://github.com/matroxel/roman_imsim).

## References

For more information about [GalSim](https://github.com/GalSim-developers/GalSim), please see its README and documentation.

For more information about Euclid, please see the [Euclid Consortium website](https://www.euclid-ec.org/) and papers linked from there.

Attribution for software and data used by particular routines in this library is given in the docstring for the relevant routine.


## Installation

GalSim-Euclid-Like requires `python>3.10` and the following dependencies (see NOTE below):
```
numpy>=1.17,
galsim>=2.5,
astropy>=2.0,
```

NOTE:
At the moment `GalSim` has to be installed manually to have access to chromatic PSF interpolation that we need (see [here](https://github.com/GalSim-developers/GalSim/pull/1296)).  
To install GalSim:
```
git clone git@github.com:GalSim-developers/GalSim.git
cd GalSim
git checkout 1294_from_images
pip install .
```

The source code for GalSim-Euclid-Like has not been published to pypi. To install from source code:
```
git clone git@github.com:GalSim-developers/GalSim-Euclid-Like.git
```
and install by running
```
conda create -n euclidlike python=3.10
cd GalSim-Euclid-Like
conda activate euclidlike
pip install -e .
```

To make sure the installation is successful, do the following:
```
$ python
>>> import euclidlike
>>> euclidlike.getBandpasses()
```
## Downloading relevant data
The Euclid-like PSF is constructed from precomputed oversampled images on a grid in focal plane position and wavelength. To use the PSF functions within `GalSim-Euclidlike`, the images must be downloaded by running:
```
euclidlike_download_psf
```
in the terminal after installation of `GalSim-Euclidlike`. To install in an alternative directory to the default, use the `--dir` argument. 

## Getting started

Please see the examples/ directory for demos illustrating the use of this code.

## Communicating with the developers

Feel free to [open a GitHub issue](https://github.com/GalSim-developers/GalSim-Euclid-Like/issues) to reach the developers with questions, comments, and bug reports.  New contributors are also welcome and can indicate their interest in developing this code base through the Issues.

## Attribution

This software is open source and may be used according to the terms of its [license](LICENSE).

When using this software, please provide the URL to the repository in the resulting paper or note.  Once there is a Zenodo DOI or journal article, this README will be updated and we will ask those using the code in their research to cite the relevant journal article.

