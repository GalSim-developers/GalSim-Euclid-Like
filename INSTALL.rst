Installation Instructions
=========================
The GalSim-Euclid-Like repository contains two python libraries and currently supports Python versions 3.10 and above.

System requirements: Given the heavy dependance on GalSim, GalSim-Euclid-Like currently only supports Linux and Mac OSX. For 
further details on system requirements for GalSim see `GalSim Installation <https://github.com/GalSim-developers/GalSim/blob/main/INSTALL.rst>`_.

Dependencies
------------

GalSim-Euclid-Like requires ``python>=3.10`` and the following dependencies::

    numpy>=1.17,
    galsim>=2.6,
    astropy>=2.0
                                                                              
Install the code
---------------

The source code for GalSim-Euclid-Like has not been published to pypi. To install from source code::

    git clone git@github.com:GalSim-developers/GalSim-Euclid-Like.git

and install by running::

    conda create -n euclidlike python=3.10
    cd GalSim-Euclid-Like
    conda activate euclidlike
    pip install .

To make sure the installation is successful, do the following::

    $ python
    >>> import euclidlike
    >>> euclidlike.getBandpasses()

GalSim-Euclid-Like ImSim dependencies
-------------------------------------

To use the euclidlike_imsim package part of GalSim-Euclid-Like, other dependencies are required.
First, the `LSST Science pipelines <https://pipelines.lsst.io/index.html#>`_ need to be installed. The easiest way to install it (on Linux distribution) is using Stackvana as follow. For more options, we refer to the `installation instructions <https://pipelines.lsst.io/index.html#installation>`_::

    conda install stackvana

In addition to the LSST Science pipelines, the `SkyCatalogs <https://lsstdesc.org/skyCatalogs/>`_ needs to be installed from an independent fork to make use of the Euclid photometry. This can be achieve by running the command::

    pip install git+https://github.com/aguinot/skyCatalogs.git@euclid_band

SkyCatalogs requires other dependencies::

    conda install -c conda-forge dust_extinction

or::

    pip install dust_extinction

Extra files are also required, they can be downloaded following the instructions `here <https://lsstdesc.org/imSim/install.html#install-needed-data-files>`_. The path to those files needs to be set in a environment variable::

    export SIMS_SED_LIBRARY_DIR=path/to/rubin_sim_data/sims_sed_library

In a conda environment it is possible to automatically set environment variable upon environment activation by doing::

    conda env config vars set SIMS_SED_LIBRARY_DIR="path/to/rubin_sim_data/sims_sed_library"
