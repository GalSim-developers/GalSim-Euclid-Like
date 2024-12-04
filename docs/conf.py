# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys
sys.path.insert(0,os.path.abspath('.'))

import euclidlike
import euclidlike_imsim

project = 'GalSim-Euclid-Like'
copyright = '2024, Rachel Mandelbaum, Axel Guinot, Federico Berlfein, Andy Park, Xiangchong Li, Michael Troxel, Tianqing Zhang'
author = 'Rachel Mandelbaum, Axel Guinot, Federico Berlfein, Andy Park, Xiangchong Li, Michael Troxel, Tianqing Zhang'

# red version and release
version = '.'.join(map(str,euclidlike.__version_info__[:2]))
release = euclidlike.__version__

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = []

templates_path = ['_templates']
exclude_patterns = ['**/.ipynb_checkpoints']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output


extensions = [
    'sphinx.ext.autodoc',  # Generate documentation from docstrings
    'sphinx.ext.viewcode',  # Show source code in the documentation
    'sphinx.ext.napoleon',  # Support for Google-style and NumPy-style docstrings
    'sphinx.ext.autosectionlabel',
    'gh-link',
]





# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# Custom sidebar templates, must be a dictionary that maps document names
# to template names.
#
# The default sidebars (for documents that don't match any pattern) are
# defined by theme itself.  Builtin themes are using these templates by
# default: ``['localtoc.html', 'relations.html', 'sourcelink.html',
# 'searchbox.html']``.
#
# html_sidebars = {}

# https://michaelgoerz.net/notes/extending-sphinx-napoleon-docstring-sections.html
# -- Extensions to the  Napoleon GoogleDocstring class ---------------------

# The master toctree document.
master_doc = 'index'

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# The reST default role (used for this markup: `text`) to use for all
# documents.
default_role = 'any'

from sphinx.ext.napoleon.docstring import GoogleDocstring

# first, we define new methods for any new sections and add them to the class
def parse_keys_section(self, section):
    return self._format_fields('Keys', self._consume_fields())
GoogleDocstring._parse_keys_section = parse_keys_section

def parse_attributes_section(self, section):
    return self._format_fields('Attributes', self._consume_fields())
GoogleDocstring._parse_attributes_section = parse_attributes_section

def parse_class_attributes_section(self, section):
    return self._format_fields('Class Attributes', self._consume_fields())
GoogleDocstring._parse_class_attributes_section = parse_class_attributes_section

# we now patch the parse method to guarantee that the the above methods are
# assigned to the _section dict
def patched_parse(self):
    self._sections['keys'] = self._parse_keys_section
    self._sections['class attributes'] = self._parse_class_attributes_section
    self._unpatched_parse()
GoogleDocstring._unpatched_parse = GoogleDocstring._parse
GoogleDocstring._parse = patched_parse