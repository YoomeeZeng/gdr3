# -*- coding: utf-8 -*-

import sys
sys.path.append('../')

from gdr3binaryorbits.orbits import NSS
import numpy as np

star=NSS()
star.query_source('4714104568778128256')
star.query_nss('SB1')
#star.plot_gaia_sb1()
star.predict_next_rvminmax()
star.draw_from_sb1_model(draws=500)
star.plot_gaia_sb1_draws()
star.get_sb1_fm_dist()
star.get_asini_dist()

rvs=np.array([31.372,32.831,33.913,18.348,27.087,20.135])
times=np.array([2459873.72357,2459874.72558,2459894.63744,2455171.92410880,2459914.6283400003,2459935.57851])
rv_errs=np.ones(len(rvs))*0.05
data_src=['CHIRON','CHIRON','CHIRON','RAVE','CHIRON','CHIRON']

star.load_rv_observations(times,rvs,rv_errs,data_src)
star.plot_nss_sol_vs_data()

times=np.linspace(2459972.50476,2459972.50476+10,10)
a=star.get_predicted_sb1_rvs(times)
print(a)


star=NSS()
star.query_coords(267.65542,-21.50428)
star.query_nss('SB1')
star.predict_next_rvminmax()
#star.plot_gaia_sb1()
star.draw_from_sb1_model(draws=500)
star.plot_gaia_sb1_draws()
star.get_sb1_fm_dist()
star.get_asini_dist()

star=NSS()
star.query_source('5853193426917488128')
star.query_nss('SB1')