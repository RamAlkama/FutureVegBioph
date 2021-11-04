# "Vegetation-based climate mitigation in a warmer and greener World" published in Nature Comunications 2021

This study was made in 3 steps. 

# 1  Estimte the sensitivity dT/dLAI as fuction of evaporation snow cover and solar radiation from observation 
First we use "Job.py" python program that import a function "nearby" with is a fortran program "nearby.f90" that is compiled using f2py in order to estimate dT/dLAI.
Second we use "evap_solarRadiatio.py" and "snow_solarRadiation.py" to get dt/dLAI as fuction of evaporatio-SolarRadiatio and SnowCover-SolarRadiation respectevely.
 
# 2  Use this sensitivity to estimate the impact of future greening on tenperature (this is based on CMIP6 experiments)
We use the observed sensitivity 

# 3  Estimate the biochemical effect of simulated future greening

