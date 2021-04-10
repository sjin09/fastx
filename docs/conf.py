"""Sphinx configuration."""
from datetime import datetime


project = "fastx"
author = "Sangjin Lee"
copyright = f"{datetime.now().year}, {author}"
extensions = ["sphinx.ext.autodoc", "sphinx.ext.napoleon"]
autodoc_typehints = "description"
