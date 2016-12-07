#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function, absolute_import, unicode_literals
from radconv import RadConv

# Initialize
rc = RadConv()
rc.p_broad = True
rc.q_init = 0.1
rc.q_eps = 0.
rc.t_eps = 0.
rc.moistconv = True
rc.run()
rc.animate(rc.t, rc.q, rc.lw_tau)


# Add flux from ocean to atmosphere
# proportional to difference between vapor pressure and saturation vapor pressure
# Set q in stratosphere equal to q of tropopause? Or, to be simple, zero it out