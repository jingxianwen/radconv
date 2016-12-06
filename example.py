#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function, absolute_import, unicode_literals
from radconv import RadConv

# Initialize
rc = RadConv()
rc.q_init = 0.01
rc.moistconv = True
rc.run()
rc.animate(rc.t, rc.q)