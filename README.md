<!-- Improved compatibility of back to top link: See: https://github.com/othneildrew/Best-README-Template/pull/73 -->
<a name="readme-top"></a>
<!--
*** Thanks for checking out the Best-README-Template. If you have a suggestion
*** that would make this better, please fork the repo and create a pull request
*** or simply open an issue with the tag "enhancement".
*** Don't forget to give the project a star!
*** Thanks again! Now go create something AMAZING! :D
-->

#gdr3binaryorbits

<!-- ABOUT THE PROJECT -->
## About The Project

Python code to visualize the radial velocity orbits published in Gaia DR3.

Note: Currently only supports queries for SB1 orbits. SB2 orbits will be added shortly!

Details on the NSS orbits: [Gaia DR3 NSS_TWO_BODY_ORBIT Documentation](https://gea.esac.esa.int/archive/documentation/GDR3/Gaia_archive/chap_datamodel/sec_dm_non--single_stars_tables/ssec_dm_nss_two_body_orbit.html)

<!-- GETTING STARTED -->

### Installation

1. Clone the repo
   ```sh
   git clone https://github.com/tjayasinghe/gdr3binaryorbits.git
   ```
2. Install with python
   ```sh
   python setup.py install
   ```


<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- USAGE EXAMPLES -->
## Usage

Here is an example on how to use this project to retrieve and visualize an SB1 RV orbit published in Gaia DR3.

1. Create a NSS object
   ```sh
   star=NSS()
   ```
2. A star can be loaded either through a cone search with RAJ2000 and DEJ2000 coordinates or through a direct query based on the Gaia DR3 source_id.
   ```sh
   star.query_source('5853193426917488128') #Query by source_id
   star.query_coords(219.41107637177004,-63.36791102568197) #Cone search
   ```
   
4. Once the star is loaded, query the Gaia DR3 NSS database and search for a SB1 orbit.
   ```sh
   star.query_nss('SB1')  
   ```
   
   We can visualize the orbit
     ```sh
   star.plot_gaia_sb1()
   ``` 
   ![image](https://user-images.githubusercontent.com/20095290/215894156-09a84b90-76c6-493b-b9a1-d782b43997ef.png)

5. We can better understand the quality of the published orbit by looking at the published uncertainties. Draw from the published Gaia DR3 posteriors for the RV orbital parameters

   ```sh
   star.draw_from_sb1_model(draws=500)  
   ```
   Visualize the orbit with the 500 orbits drawn from the posteriors.
![image](https://user-images.githubusercontent.com/20095290/215894345-306d8b02-b8fb-4570-a0ad-4f8260d0f4b3.png)

6. We can also look at the distribution of asini and f(M) for SB1 orbits.

![image](https://user-images.githubusercontent.com/20095290/215894577-3262defe-4634-45bb-af84-edbf162581c2.png)


![image](https://user-images.githubusercontent.com/20095290/215894514-ad047bef-334b-46d7-8b73-e0f6fb816e75.png)



<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE.txt` for more information.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

