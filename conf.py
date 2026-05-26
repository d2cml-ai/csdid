# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

# Make the package importable so autodoc can pick it up if you switch
# the Reference page to use ``.. autoclass::`` / ``.. autofunction::``.
sys.path.insert(0, os.path.abspath(".."))


# -- Project information -----------------------------------------------------

project = "CSDID"
copyright = "2024, D2CML Team"
author = "Alexander Quispe, Carlos Guevara, Brantly Callaway, Pedro H.C. Sant'Anna"


# -- General configuration ---------------------------------------------------

extensions = [
    "nbsphinx",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.mathjax",
]

autosummary_generate = True
autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
}
napoleon_google_docstring = False
napoleon_numpy_docstring = True

nbsphinx_execute = "never"
nbsphinx_allow_errors = True

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "**.ipynb_checkpoints"]


# -- Options for HTML output -------------------------------------------------

html_theme = "alabaster"
html_static_path = ["_static"]
