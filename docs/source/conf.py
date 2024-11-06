# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys
sys.path.insert(0,os.path.abspath('.'))

project = 'GalSim-Euclid-Like'
copyright = '2024, Rachel Mandelbaum, Axel Guinot, Federico Berlfein, Andy Park, Xiangchong Li, Michael Troxel, Tianqing Zhang'
author = 'Rachel Mandelbaum, Axel Guinot, Federico Berlfein, Andy Park, Xiangchong Li, Michael Troxel, Tianqing Zhang'
release = '0.0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = []

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']

extensions = [
    'sphinx.ext.autodoc',  # Generate documentation from docstrings
    'sphinx.ext.viewcode',  # Show source code in the documentation
    'sphinx.ext.napoleon',  # Support for Google-style and NumPy-style docstrings
]
