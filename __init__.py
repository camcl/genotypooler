# -*- coding: utf-8 -*-
# flake8: noqa
from __future__ import absolute_import, print_function, division

import os, sys, importlib

home_dir = os.path.expanduser("~")
proj_dir = os.path.join(home_dir, '1000Genomes')
sys.path.insert(0, proj_dir)