#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
sys.path.insert(0, os.path.abspath('..'))

from globalconf import *

project = 'Lumol'

latex_documents = [
    (master_doc, 'Lumol.tex', 'Lumol user manual', author, 'howto'),
]
