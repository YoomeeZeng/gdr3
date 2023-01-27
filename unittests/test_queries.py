# -*- coding: utf-8 -*-

import sys
sys.path.append('../')

from gdr3binaryorbits.orbits import NSS


star=NSS()
star.query_source('4714104568778128256')
print(star.gdr3_source,star.gdr3_solution,star.gdr3_plx)
star.query_nss('SB1')
print(star.period,star.ecc,star.K1,star.gamma,star.t_peri)