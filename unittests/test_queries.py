# -*- coding: utf-8 -*-

import sys
sys.path.append('../')

from gdr3binaryorbits.orbits import NSS
import numpy as np

star=NSS()
star.query_source('5853193426917488128')
star.query_nss('SB1')
star.plot_gaia_sb1()
star.predict_next_rvminmax()
star.draw_from_sb1_model(draws=500)
star.plot_gaia_sb1_draws()
star.get_sb1_fm_dist()
star.get_asini_dist()
