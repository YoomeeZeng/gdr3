# -*- coding: utf-8 -*-

import unittest
import sys
sys.path.append('../')

class QuickTest(unittest.TestCase):
    
    def test_dependencies(self):
        import pandas
        import numpy
        import astropy
        import matplotlib

    def test_create_NSSObject(self):
        
        from gdr3binaryorbits.orbits import NSS
        star=NSS()

    def test_query_source(self):
        
        from gdr3binaryorbits.orbits import NSS
        star=NSS()
        star.query_source('5853193426917488128')

    def test_query_coords(self):
        
        from gdr3binaryorbits.orbits import NSS
        star=NSS()
        star.query_coords(219.41107637177004,-63.36791102568197)
        
    def test_query_nss_sb1(self):
        
        from gdr3binaryorbits.orbits import NSS
        star=NSS()
        star.query_coords(219.41107637177004,-63.36791102568197)
        star.query_nss('SB1')
        
    def test_draw_from_orbit(self):
        
        from gdr3binaryorbits.orbits import NSS
        star=NSS()
        star.query_coords(219.41107637177004,-63.36791102568197)
        star.query_nss('SB1')        
        star.draw_from_sb1_model(draws=100)
        
if __name__ == '__main__':
    unittest.main()# -*- coding: utf-8 -*-

