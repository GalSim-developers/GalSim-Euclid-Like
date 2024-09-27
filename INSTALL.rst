Installation Instructions
==================
GalSim-Euclid-Like is a python module implemented mainly in Python. It currently supports Python versions 3.10 and above.

System requirements: Given the heavy dependance on GalSim, GalSim-Euclid-Like currently only supports Linux and Mac OSX. For 
further details on system requirements for GalSim see `GalSim Installation <https://github.com/GalSim-developers/GalSim/blob/main/INSTALL.rst>`_.

Dependancies
----------------

GalSim-Euclid-Like requires ``python>=3.10`` and the following dependencies::

    numpy>=1.17,
    galsim>=2.6,
    astropy>=2.0
                                                                              
Installation
----------------
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
