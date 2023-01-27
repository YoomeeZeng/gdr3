# -*- coding: utf-8 -*-

import sys
sys.path.append('../')

from gdr3binaryorbits.orbits import NSS
import numpy as np

star=NSS()
star.query_source('4714104568778128256')
print(star.gdr3_source,star.gdr3_solution,star.gdr3_plx)
star.query_nss('SB1')
print(star.period,star.ecc,star.K1,star.gamma,star.t_peri)
#star.plot_gaia_sb1()
star.draw_from_sb1_model(draws=200)
star.plot_gaia_sb1_draws()

rvs=np.array([31.372,32.831,33.913,18.348,27.087,20.135])
times=np.array([2459873.72357,2459874.72558,2459894.63744,2455171.92410880,2459914.6283400003,2459935.57851])
rv_errs=np.ones(len(rvs))*0.05
data_src=['CHIRON','CHIRON','CHIRON','RAVE','CHIRON','CHIRON']

star.load_rv_observations(times,rvs,rv_errs,data_src)
star.plot_nss_sol_vs_data()
'''
star=NSS()
star.query_coords(267.65542,-21.50428)
print(star.gdr3_source,star.gdr3_solution,star.gdr3_plx)
star.query_nss('SB1')
print(star.period,star.ecc,star.K1,star.gamma,star.t_peri)
#star.plot_gaia_sb1()
star.draw_from_sb1_model(draws=200)
star.plot_gaia_sb1_draws()'''