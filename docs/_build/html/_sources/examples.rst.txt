Examples
=========

The ``GalSim-Euclid-Like/examples`` directory contains example files for how use the euclidlike module.  

End-to-end demo 
---------------

:gh-link:`end_to_end_demo.py <examples/end_to_end_demo.py>`

This first demo is the euclidlike-equivalent of  `demo13 <https://github.com/GalSim-developers/GalSim/blob/main/examples/demo13.py>`_ in ``GalSim``. This demo uses the Euclid-like PSF, WCS, and background noise to produce a realistic scene of galaxies and stars as observed from a Euclid-like Telescope. 

**Features introduced in the Python file**:

- euclidlike.getBandpasses(AB_zeropoint)
- euclidlike.getWCS(world_pos, CCDs_CCD, date)
- euclidlike.getPSF(use_CCD, filter_name, wcs)
- euclidlike.getSkyLevel(bandpass, world_pos)

The output generated from this file can be visualized by running the script :gh-link:`plot_VIS.py <examples/plot_VIS.py>`.


Focal Plane Layout 
------------------
:gh-link:`focal_plane_layout.ipynb <examples/focal_plane_layout.ipynb>`

This Jupyter Notebook shows the display of the focal plane used in the euclidlike package, along with the CCD centers and ID convention.

