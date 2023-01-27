# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from astroquery.vizier import Vizier
from astroquery.gaia import Gaia
import astropy.units as u
import astropy.coordinates as coord
import warnings
warnings.filterwarnings("ignore")

Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"

def get_mag_err(flux,flux_err):  
    
    '''
    Returns the error in magnitude given flux,flux_err
    
    Parameters
    ----------
    flux: float
        Flux
    flux_err: float
        Error in flux  
        
    Returns
    ----------
    mag_err: float
        Error in magnitude
    '''
    
    mag_err=(2.5/np.log(10))*(flux_err/flux)
    
    return mag_err

class NSS:
    
    def __init__(self):
        
        pass

    def query_coords(self,ra,dec,query_radius=2):
 
        self.ra = ra
        self.dec = dec
        
        #Create a SkyCoord
        try:
            self.coord=coord.SkyCoord(ra=ra, dec=dec,unit=(u.deg, u.deg),frame='icrs')
        except:
            print("Error creating SkyCoord. Are your coordinates in degrees?")
            
        try:
            
            gdr3_star_query = Gaia.cone_search_async(self.coord, query_radius*u.arcsec)
            gdr3_result=gdr3_star_query.get_results().to_pandas()
           
            if len(gdr3_result)!=0:
                
                print(f"Gaia DR3 match found {round(gdr3_result['dist'][0],2)} arcsec away.")
                
                self.gdr3_query_result=gdr3_result
                #Update Coordinates
                self.ra=gdr3_result['ra'][0]
                self.dec=gdr3_result['dec'][0]

                #Save GDR3 attributes
                self.gdr3_source=gdr3_result['source_id'][0]
                self.gdr3_plx=gdr3_result['parallax'][0]
                self.gdr3_plx_err=gdr3_result['parallax_error'][0]                      
                self.gdr3_rv=gdr3_result['radial_velocity'][0]       
                self.gdr3_rv=gdr3_result['radial_velocity_error'][0]             
                self.gdr3_solution=gdr3_result['solution_id'][0]  
                                            
                #Save magnitudes
                self.Gmag=gdr3_result['phot_g_mean_mag'][0]  
                self.Gmag_err=get_mag_err(gdr3_result['phot_g_mean_flux'][0] ,gdr3_result['phot_g_mean_flux_error'][0] )      
                self.BPmag=gdr3_result['phot_bp_mean_mag'][0]  
                self.BPmag_err=get_mag_err(gdr3_result['phot_bp_mean_flux'][0] ,gdr3_result['phot_bp_mean_flux_error'][0] )            
                self.RPmag=gdr3_result['phot_rp_mean_mag'][0]  
                self.RPmag_err=get_mag_err(gdr3_result['phot_rp_mean_flux'][0] ,gdr3_result['phot_rp_mean_flux_error'][0] ) 
                                 
                
            else:
                self.gdr3_query_success=False
                print("Gaia DR3 query returned 0 matches.")
                
            
        
        except:
            self.gdr3_query_success=False
            print('Querying Gaia DR3 failed.')     
            

    def query_source(self,source_id):
        
        self.gdr3_source=source_id
        try:
            source_query = Gaia.launch_job(f"select * from gaiadr3.gaia_source where source_id={self.gdr3_source}")
            gdr3_result=source_query.get_results().to_pandas()

            if len(gdr3_result)>0:
                print(f"Gaia DR3 match found for Source={source_id}.")
                                
                #Update Coordinates
                self.ra=gdr3_result['ra'][0]
                self.dec=gdr3_result['dec'][0]
 
                #Save GDR3 attributes
                self.gdr3_plx=gdr3_result['parallax'][0]
                self.gdr3_plx_err=gdr3_result['parallax_error'][0]                             
                self.gdr3_rv=gdr3_result['radial_velocity'][0]       
                self.gdr3_rv=gdr3_result['radial_velocity_error'][0]  
                self.gdr3_solution=gdr3_result['solution_id'][0]  
                
                #Save magnitudes
                self.Gmag=gdr3_result['phot_g_mean_mag'][0]  
                self.Gmag_err=get_mag_err(gdr3_result['phot_g_mean_flux'][0] ,gdr3_result['phot_g_mean_flux_error'][0] )      
                self.BPmag=gdr3_result['phot_bp_mean_mag'][0]  
                self.BPmag_err=get_mag_err(gdr3_result['phot_bp_mean_flux'][0] ,gdr3_result['phot_bp_mean_flux_error'][0] )            
                self.RPmag=gdr3_result['phot_rp_mean_mag'][0]  
                self.RPmag_err=get_mag_err(gdr3_result['phot_rp_mean_flux'][0] ,gdr3_result['phot_rp_mean_flux_error'][0] ) 
                
            else:
                print('GSPC returned 0 matches.')                
                
        except:
            print('Querying GSPC Failed.')
            
            
    def query_nss(self,solution_type):
        
        self.solution_type=solution_type
        
        try:
            nss_query = Gaia.launch_job(f"select * from gaiadr3.nss_two_body_orbit as t1 "
                                              f"join gaiadr3.gaia_source as t2 using (source_id) where source_id={self.gdr3_source}")
            nss_table=nss_query.get_results().to_pandas()
            
            if len(nss_table)>0:
                nss_table=nss_table.loc[nss_table['nss_solution_type']==solution_type.encode()]
                nss_table.reset_index(drop=True,inplace=True)
                
                if len(nss_table)>0:
                    if solution_type=='SB1' or solution_type=='SB1C':
                        print(f'Found orbit solution of type {solution_type} for Source={self.gdr3_source}')                        
                        #Orbit params
                        self.period=nss_table['period'][0]
                        self.period_err=nss_table['period_error'][0]                        
                        self.ecc=nss_table['eccentricity'][0]
                        self.ecc_err=nss_table['eccentricity_error'][0]
                        self.K1=nss_table['semi_amplitude_primary'][0]       
                        self.K1_err=nss_table['semi_amplitude_primary_error'][0]
                        self.gamma=nss_table['center_of_mass_velocity'][0]
                        self.gamma_err=nss_table['center_of_mass_velocity_error'][0] 
                        self.t_peri=nss_table['t_periastron'][0]
                        self.t_peri_err=nss_table['t_periastron_error'][0]
                        
                        
                        #Solution params
                        self.total_rvs_primary=int(nss_table['rv_n_obs_primary'][0]       )    
                        self.total_good_rvs_primary=int(nss_table['rv_n_good_obs_primary'][0])
                        self.solution_gof=nss_table['goodness_of_fit'][0]      
                        self.solution_efficiency=nss_table['efficiency'][0]        
                        self.solution_significance=nss_table['significance'][0]     
                        self.solution_flags=nss_table['flags'][0]      
                          
                #Add more solution types: SB2                                                   
                else:
                    print(f'No orbit solutions of type {solution_type} for Source={self.gdr3_source}')
            
            else:
                print('GDR3 NSS query returned 0 matches.')                
                
        except:
            print('Querying GDR3 NSS Failed.')
            
                
            
    