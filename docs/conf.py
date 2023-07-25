import os
import re

project = "Sycomore"
copyright = "2019-2023, Université de Strasbourg-CNRS"
author = "Julien Lamy"

here = os.path.abspath(os.path.dirname(__file__))

html_theme = "furo"
html_static_path = ["_static"]
html_css_files = ["css/style.css"]
html_title = project

extensions = [
    "breathe", "sphinx.ext.autodoc", "sphinx.ext.doctest", "sphinx.ext.mathjax"]

autodoc_default_options = {
    "members": True,
    "special-members": True,
    "undoc-members": True,
    "exclude-members": (
        "__annotations__, __getstate__, "
        "__module__, __setstate__, __weakref__")
}

breathe_projects = { project: "_build/doxygen/xml" }
breathe_default_project = project
breathe_domain_by_extension = { "h" : "cpp" }
breathe_default_members = ("members", "undoc-members")

mathjax3_config = {
    "HTML-CSS": {"scale": 95},
}
