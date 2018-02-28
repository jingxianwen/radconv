#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function, absolute_import, \
    unicode_literals
from radconv import RadConv

# Initialize
rc = RadConv()
rc.p_broad = True
rc.wv_tau = 0.001
rc.load()
rc.animate(rc.t, rc.q, rc.flux_sw, rc.flux_rad)


# linear tau such that CO2 amount is constant
# (tau decreases linearly with surf pressure)


# reduce linear tau and boost water tau so that they are comparable
# make surf temp 285
