# "Vegetation-based climate mitigation in a warmer and greener World"  

To estimate the plant biophysics in response to the future climate we first estimate the temperature sensitivity to LAI from observations over 2003-2014 under different combination of key environmental drivers such as snow cover, solar radiation and evaporation. Such datasets are used at 0.05 degree spatial resolution. We then used this sensitivity with future LAI values and climate conditions as simulated by an ensemble of CMIP6 experiments under different SSPs at their common 2x2 degree spatial resolution, in order to estimate the future evolution of plant biophysical impacts on climate. To account for the relationship between water use efficiency and atmospheric CO2 concentration, we used evaporation rates rather than soil moisture as a driver of dT/dLAI. In fact, plant evaporation simulated in climate models already includes the effect of CO2 fertilization on stomatal conductance.
This study was made using the in 3 folowimg steps. 

# 1  Estimte the sensitivity dT/dLAI as fuction of evaporation snow cover and solar radiation from observation 

The goal of this part of the work is to quantify the local climate impacts of observed changes in vegetation density over the period 2003-2014, for which combined observation of leaf area index (GLASS) and air temperature (inferred from satellite observations by Hooker et al. 2018) are available. For this purpose, the variation of air surface temperatures induced by the change in LAI has been factored out from the natural climate variability using the temperature signal from neighbouring areas with stable LAI. 
Our methodology is similar to the one presented in Alkama and Cescatti 2016, based on Equation 1, which assumes that, for a given grid cell, the difference in temperature between two years is equal to the sum of the temperature variation induced by LAI change (dT(lai)) plus the residual signal (dT(res)) due to the natural inter-annual climate variability.

==>  dT=dT(res)+dT(lai)         (1)

From Equation 1 it follows that the lai  can be quantified as the difference between the observed temperature variations (dT) and the natural climate variability (dT(res)). This latter term is estimated from nearby reference grid cells located within 50km distance (di) and with stable LAI (i.e. less than 0.1 m2/m2 variation in LAI during the observation period). For these grid cells, we can assume dT(lai)=0 and consequently dT(res)=dT. The residual temperature signal is therefore the temporal variation in temperature observed in the surroundings of the target grid cell in areas with stable land cover. We used an inverse distance weighting to estimate  from the n reference grid cells (50 km radius), according to Equation 2.
 dT(res)=sum(dT/di)/sum(1/di)     (2)
From equations 1 and 2, we estimate  for all grid cells where we observe LAI change larger than 0.1 m2/m2. This calculation of these differences is repeated between two years for each individual month, for all the 66 pairs of years available for the period 2003-2014. 

To do this we use "job.py" python program that import a function "nearby" which is based on a fortran program "nearby.f90" that should be compiled using f2py. The result is computaion of dT/dLAI.

We then used Obs.py python program to get the average dT/dLai from the the all 66 pairs of years available. The result of this program is stored in "Obs.nc" that is shown in Fig1c.

Once the monthly local sensitivity is estimated, we split the data in two regions (with and without snow using 1% monthly snow cover as threshold). The snow-covered areas are known to be dominated by the radiative effect (i.e. due to the contrast in albedo between vegetation and snow), while the snow-free areas are dominated by the partitioning of the available energy in turbulent fluxes (i.e. evaporation versus sensible heat). The sensitivity dT/dLAI in the first areas are then expressed as a function of monthly snow cover (SC in %) estimates from MODIS (MYD10CM.006) and solar radiation (SWdown in W/m2) from the ERA5 reanalysis using the bivariate quadratic least square regression and the resulted function is as follow.  

This is done using "evap_solarRadiatio.py" and "snow_solarRadiation.py" which allow us to have dt/dLAI as fuction of evaporatio-SolarRadiatio and SnowCover-SolarRadiation respectevely.
 
# 2  Use this sensitivity to estimate the impact of future greening on temperature (this is based on CMIP6 experiments)
We use the sensitivity with simulated LAI, evaporation, solar radiation and snow cover for each of the four SSPs (SSP585,SSP370, SSP245 and SSP126) over 2015-2100. 
Program and figures are besed on python notebook programm "main_figures.ipynb".

# 3  Biochemical effect of vegetation greening on air temperature
A previous study (Leduc et al. 2016 Nat Clim Change) shows a strong linear relationship between atmospheric carbon concentration and regional surface air temperature. Here, we combined the strong linearity of the regional climate response over most land regions presented from Leduc et al.59 and the simulated variation in vegetation carbon stock by CMIP6 climate models, to drive the global-scale biochemical climate impacts. Basically, Leduc et al. 2016 find an increase of land temperature of 2.2±0.5o per 1 Terra ton of carbon (Tt C) in the atmosphere. In our case we used the total increase of carbon in plants (ΔB in Tt C) between 2015 and 2100 coming from the CMIP6 archive to estimate the biochemical effect.
Program and figures are besed on python notebook programm "main_figures.ipynb".
