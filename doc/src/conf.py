#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import toml
import os
import sphinx_bootstrap_theme

# -- General configuration ------------------------------------------------

ROOT = os.path.join(os.path.dirname(__file__), "..")

# If your documentation needs a minimal Sphinx version, state it here.
needs_sphinx = '1.8'

extensions = [
    'sphinx.ext.mathjax',
    'sphinx.ext.ifconfig',
]

templates_path = [os.path.join(ROOT, "templates")]
source_suffix = '.rst'
master_doc = 'index'

copyright = '2019, the lumol developers'
author = 'The lumol developers'
project = 'Lumol'


def version():
    parsed = toml.loads(open(os.path.join("..", "..", "Cargo.toml")).read())
    release = parsed["package"]["version"]
    version = ".".join(release.split(".")[:2])
    return version, release


version, release = version()

language = None

exclude_patterns = []

pygments_style = 'default'
highlight_language = 'toml'


def setup(app):
    app.add_stylesheet("lumol.css")


# -- Options for HTML output ----------------------------------------------

html_theme = 'bootstrap'
html_theme_path = sphinx_bootstrap_theme.get_html_theme_path()

html_theme_options = {
    'navbar_site_name': "Navigation",
    'navbar_pagenav': False,
    'source_link_position': None,
    'bootswatch_theme': "sandstone",
    'bootstrap_version': "3",
}

html_sidebars = {
    '**': ['sidebar-toc.html', 'searchbox.html']
}

html_static_path = [os.path.join(ROOT, "static", "lumol.css")]

# -- Options for LaTeX output ---------------------------------------------

latex_elements = {
    'papersize': 'a4paper',
}

latex_documents = [
    (master_doc, 'Lumol.tex', 'Lumol user manual', author, 'howto'),
]
