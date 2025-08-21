import os
import sys
sys.path.insert(0, os.path.abspath('exts'))

doctest_test_doctest_blocks = ""
doctest_global_setup = '''
from rdkit import Chem
from rdkit.Chem import AllChem
'''

# -- Project information -----------------------------------------------------

project = 'RDKit'
copyright = '2025, Greg Landrum and others'
author = 'Greg Landrum'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx_rtd_theme',
    'sphinx_copybutton',
]

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']