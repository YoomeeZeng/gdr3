# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from astroquery.vizier import Vizier
from astroquery.gaia import Gaia
import astropy.units as u
import astropy.coordinates as coord
from astropy.time import Time
from datetime import datetime

import warnings
warnings.filterwarnings("ignore")

from .misc_utils import *
from .plot_utils import *

Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"

class NSS:
    
    def __init__(self):
        self.ra=np.nan
        self.dec=np.nan
        self.gdr3_source=np.nan
        self.rv_orbit_loaded=False
        self.orbits_sampled=False
        self.observations_loaded=False
        self.current_jd= Time(datetime.now().isoformat(), format='isot').jd

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
                print('Gaia DR3 Query returned 0 matches.')                
                
        except:
            print('Gaia DR3 Query Failed.')
            
            
    def query_nss(self,solution_type):
        
        self.solution_type=solution_type
        
        try:
            nss_query = Gaia.launch_job(f"select * from gaiadr3.nss_two_body_orbit as t1 "
                                              f"join gaiadr3.gaia_source as t2 using (source_id) where source_id={self.gdr3_source}")
            nss_table=nss_query.get_results().to_pandas()

            if len(nss_table)>0:
                nss_table=nss_table.loc[nss_table['nss_solution_type']==solution_type]
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
                        self.arg_per=nss_table['arg_periastron'][0]
                        self.arg_per_err=nss_table['arg_periastron_error'][0]                         
                        self.t_peri=nss_table['t_periastron'][0]
                        self.t_peri_err=nss_table['t_periastron_error'][0]
                        self.t_peri_jd=epoch=2457388.5+self.t_peri
                        
                        #Lucy-Sweeney test for near circular binaries
                        self.ecc_over_err=self.ecc/self.ecc_err
                        if not self.ecc_over_err>5 and self.ecc<0.05:
                            print(f'Orbital eccentricity might be negligible: e/de={round(self.ecc_over_err,2)}')
                            print('Setting arg_periastron=0, ecc=0')
                            self.arg_per=0
                            self.arg_per_err=0
                            self.ecc=0
                            self.ecc_err=0                            
                        #Definition of t_peri is time of RV max for SB1C
                        
                        #Solution params
                        self.total_rvs_primary=int(nss_table['rv_n_obs_primary'][0]       )    
                        self.total_good_rvs_primary=int(nss_table['rv_n_good_obs_primary'][0])
                        self.solution_gof=nss_table['goodness_of_fit'][0]      
                        self.solution_efficiency=nss_table['efficiency'][0]        
                        self.solution_significance=nss_table['significance'][0]     
                        self.solution_flags=nss_table['flags'][0]    
                        
                        self.params_dict={'gdr3_source':self.gdr3_source,'period':self.period,'ecc':self.ecc,'K1':self.K1,
                                          'gamma':self.gamma,'arg_per':self.arg_per,'t_peri':self.t_peri_jd}
                        
                        self.rv_orbit_loaded=True
                        
                        print('-----------------------------------')
                        print(f'P_orb={self.period}+/-{self.period_err} days')
                        print(f'e_orb={self.ecc}+/-{self.ecc_err}' )    
                        print(f'K1={self.K1}+/-{self.K1_err} km/s')  
                        print(f'v_gamma={self.gamma}+/-{self.gamma_err} km/s')        
                        print(f'omega_per={self.arg_per}+/-{self.arg_per_err} deg')  
                        print(f't_peri={self.t_peri_jd}+/-{self.t_peri_err}')
                        print('-----------------------------------')                         
                #Add more solution types: SB2                                                   
                else:
                    print(f'No orbit solutions of type {solution_type} for Source={self.gdr3_source}')
            
            else:
                print('GDR3 NSS query returned 0 matches.')                
                
        except:
            print('Querying GDR3 NSS Failed.')

    def _get_rvdf(self):
    
        phases=np.linspace(0,1,100)
        rvs=get_phased_sb1_rvs(phases,self.K1,self.ecc,self.arg_per,self.gamma)  
        
        self.rv_df=pd.DataFrame({'phase':phases,'RV':rvs}) 

    def plot_gaia_sb1(self):
    
        self._get_rvdf()
        plot_sb1_rv(self.rv_df,self.params_dict)         

    def plot_gaia_sb1_draws(self):
    
        if self.rv_orbit_loaded:
            
            if self.orbits_sampled:
                
                self._get_rvdf()
                plot_sb1_rv(self.rv_df,self.params_dict,self.sb1_draws)  
                
            else:
                
                print('This object does not have its orbit sampled. Sample with draw_from_sb1_model()')
            
        else:
            
            print('This object does not have a Gaia RV orbit. Try querying NSS.')
            
    def draw_from_sb1_model(self,draws=100):
        
        if self.rv_orbit_loaded:
            period_dist=np.random.normal(self.period,self.period_err,draws)
            ecc_dist=np.random.normal(self.ecc,self.ecc_err,draws) 
            K1_dist=np.random.normal(self.K1,self.K1_err,draws) 
            arg_per_dist=np.random.normal(self.arg_per,self.arg_per_err,draws) 
            gamma_dist=np.random.normal(self.gamma,self.gamma_err,draws)    
            t_peri_dist=np.random.normal(self.t_peri_jd,self.t_peri_err,draws) 
            
            
            self.sb1_draws=pd.DataFrame({'period':period_dist,'ecc':ecc_dist,'K1':K1_dist,'arg_per':arg_per_dist,
                                         'gamma':gamma_dist,'t_peri':t_peri_dist})
            self.sb1_draws['fm']=np.vectorize(get_fm)(self.sb1_draws.period,self.sb1_draws.ecc,self.sb1_draws.K1)
            self.sb1_draws['asini']=np.vectorize(get_asini)(self.sb1_draws.period,self.sb1_draws.ecc,self.sb1_draws.K1)            
            self.orbits_sampled=True
        else:
            print('This object does not have a Gaia RV orbit. Try querying NSS.')
            
    def predict_next_rvminmax(self):
        

        self.current_orbital_phase=((self.current_jd-self.t_peri_jd)/self.period)%1
        
        print(f'At current JD={self.current_jd}, orbital phase is {round(self.current_orbital_phase,2)}')
                
        self.rv_max,self.phase_rv_max,self.rv_min,self.phase_rv_min=get_rv_extrema(self.K1,self.ecc,self.arg_per,self.gamma)
        
        phase_diff_rvmin=self.phase_rv_min-self.current_orbital_phase
        phase_diff_rvmax=self.phase_rv_max-self.current_orbital_phase
        
        #if phase_diffs are negative, extrema are in next cycle, add 1
        
        if phase_diff_rvmin<0:
            phase_diff_rvmin+=1
        if phase_diff_rvmax<0:
            phase_diff_rvmax+=1            
            
            
        self.jd_next_rv_min=(phase_diff_rvmin*self.period)+self.current_jd
        self.jd_next_rv_max=(phase_diff_rvmax*self.period)+self.current_jd            
            
          
        print(f'Next RV maximum of {round(self.rv_max,1)} km/s will occur on JD={round(self.jd_next_rv_max,5)} in {round(self.jd_next_rv_max-self.current_jd,1)} days.')     
        print(f'Next RV minimum of {round(self.rv_min,1)} km/s will occur on JD={round(self.jd_next_rv_min,5)} in {round(self.jd_next_rv_min-self.current_jd,1)} days.')       
    
    def get_predicted_sb1_rvs(self,times):
        
        if self.rv_orbit_loaded:
            
            if self.orbits_sampled:
                
                expected_rv=[]
                expected_rv_std=[]
                
                for time in times:
                    expected_rv_i,expected_rv_std_i=get_predicted_rv(time,self.sb1_draws)
                    expected_rv.append(expected_rv_i)
                    expected_rv_std.append(expected_rv_std_i)
                    
                
                predictions_df=pd.DataFrame({'time':times,'predicted_rv':expected_rv,'predicted_rv_std':expected_rv_std})
                
                return predictions_df

            else:
                
                print('This object does not have its orbit sampled. Sample with draw_from_sb1_model()')
            
        else:
            
            print('This object does not have a Gaia RV orbit. Try querying NSS.')


    def get_sb1_fm_dist(self):
        
        if self.rv_orbit_loaded:
            
            if self.orbits_sampled:
                plot_fm_dist(self.sb1_draws,self.params_dict)
                
                self.fm_50=self.sb1_draws.fm.quantile(0.5)
                self.fm_84=self.sb1_draws.fm.quantile(0.84)
                self.fm_16=self.sb1_draws.fm.quantile(0.16)       
                self.fm_err=(np.abs(self.fm_84-self.fm_50)+np.abs(self.fm_16-self.fm_50))/2
                
                print(f'f(M)={self.fm_50} +/- {self.fm_err} Msun')
            else:
                
                print('This object does not have its orbit sampled. Sample with draw_from_sb1_model()')
            
        else:
            
            print('This object does not have a Gaia RV orbit. Try querying NSS.')
            
    def get_asini_dist(self):
        
        if self.rv_orbit_loaded:
            
            if self.orbits_sampled:
                
                if self.solution_type=="SB1" or self.solution_type=='SB1C':
                    plot_sb1_asini_dist(self.sb1_draws,self.params_dict)
                    
                    self.asini1_50=self.sb1_draws.asini.quantile(0.5)
                    self.asini1_84=self.sb1_draws.asini.quantile(0.84)
                    self.asini1_16=self.sb1_draws.asini.quantile(0.16)       
                    self.asini1_err=(np.abs(self.asini1_84-self.asini1_50)+np.abs(self.asini1_16-self.asini1_50))/2
                
                    print(f'asini_1={self.asini1_50} +/- {self.asini1_err} Rsun')
            else:
                
                print('This object does not have its orbit sampled. Sample with draw_from_sb1_model()')
            
        else:
            
            print('This object does not have a Gaia RV orbit. Try querying NSS.')
            
    def load_rv_observations(self,times,rvs,rv_errs,data_source):
        
        if self.rv_orbit_loaded:
            
            phases=((times-self.t_peri_jd)/self.period)%1 
            self.data_df=pd.DataFrame({'JD':times,'phase':phases,'RV':rvs,'RV_err':rv_errs,'obs_src':data_source})
            self.observations_loaded=True
            self.calculate_residuals_from_observations()
            
        else:
            
            print('No orbit to calculate phases.')
            self.data_df=pd.DataFrame({'JD':times,'phase':np.nan*np.zeros(len(times)),'RV':rvs,'RV_err':rv_errs,'obs_src':data_source})
            self.observations_loaded=True
        
        print(f'RV data loaded: {len(self.data_df)} epochs')
        
    def calculate_residuals_from_observations(self):
        
        if self.observations_loaded:
            
            if self.rv_orbit_loaded:
                    
                self.data_df['RV_model']=get_phased_sb1_rvs(self.data_df.phase,self.K1,self.ecc,self.arg_per,self.gamma)
                self.data_df['RV_resid']=self.data_df['RV']-self.data_df['RV_model']
                    
            else:
                
                print('This object does not have a Gaia RV orbit. Try querying NSS.')
     
        else:
            print('Observations not loaded!')
            
    def plot_nss_sol_vs_data(self):
    
        if self.observations_loaded:
            
            if self.rv_orbit_loaded:
                
                if self.orbits_sampled:
                    
                    self._get_rvdf()
                    plot_sb1_orbit_data_comparison(self.rv_df,self.data_df,self.params_dict,self.sb1_draws)
                    
                else:
                    
                    print('This object does not have its orbit sampled. Sample with draw_from_sb1_model()')
                
            else:
                
                print('This object does not have a Gaia RV orbit. Try querying NSS.')
     
        else:
            print('Observations not loaded!')
