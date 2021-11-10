import numpy as np, glob
from itertools import combinations
import scipy.linalg
import cdms2 as cdms, numpy as N
import MV2 as MV
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from mpl_toolkits.basemap import Basemap, shiftgrid, cm

cdms.setNetcdfShuffleFlag(0)
cdms.setNetcdfDeflateFlag(0)
cdms.setNetcdfDeflateLevelFlag(0)

def readDTperLAI(ys,ye,mth,lati,loni):
  f=cdms.open('/ESS_EarthObs/SATELLITE_RS/GLASS/LAI/0.05deg/monthly_GLASS_LAI_0.05_%i.nc'%(ys),'r')
  lai_1    = f('LAI',cdms.timeslice(mth-1,mth),lat=(lati-1,lati+1),lon=(loni-1,loni+1),squeeze=1)
  f.close()
  f=cdms.open('/ESS_EarthObs/SATELLITE_RS/GLASS/LAI/0.05deg/monthly_GLASS_LAI_0.05_%i.nc'%(ye),'r')
  lai_2    = f('LAI',cdms.timeslice(mth-1,mth),lat=(lati-1,lati+1),lon=(loni-1,loni+1),squeeze=1)
  f.close()
  if lati<0:
    f=cdms.open('/ESS_EarthObs/SATELLITE_RS/GLASS/LAI/AST/OUTPUT/%i-%i/Tair_scater_%i%i%02i_aqua_56S_0N.nc'%(ys,ye,ys,ye,mth))
  else:
    f=cdms.open('/ESS_EarthObs/SATELLITE_RS/GLASS/LAI/AST/OUTPUT/%i-%i/Tair_scater_%i%i%02i_aqua_0N_80N.nc'%(ys,ye,ys,ye,mth))
  dt1=f('Dtas',lat=(lati-1,lati+1),lon=(loni-1,loni+1),squeeze=1)
  f.close()
  dlai=MV.masked_equal(lai_2-lai_1,0)
  return dt1/dlai, dlai

nblat=len(range(-55,80,2))
nblon=len(range(-179,180,2))
#sensitivity=N.zeros((12,66,nblat,nblon))
sensitivity=N.zeros((66,nblat,nblon))
STD=N.zeros((66,nblat,nblon))
for mth in [11,12]:#range(1,13,1):
  print mth
  latk=-1
  for lati in range(-55,80,2):
    latk=latk+1
    lonk=-1
    for loni in range(-179,180,2): 
      lonk=lonk+1
      yk=-1
      f=cdms.open('/ESS_Datasets/EXT_ESM/CMIP6/Historical/SE/seamask.nc')
      land_mask=int(f('topo',lat=lati,lon=loni))
      f.close()
      if land_mask==1:
        for ys in range(2003,2014,1):
          for ye in range(ys+1,2015,1):
            yk=yk+1
            dtdlai,dlai=readDTperLAI(ys,ye,mth,lati,loni)
            dtdlai=MV.masked_equal(dtdlai,0)
            msk=(dtdlai/dtdlai)*(dlai/dlai)
            dtdlai=(dtdlai*msk).compressed()
            dlai=(dlai*msk).compressed()
            # mask 70% of the smallest lai change, we assume beter accuracy over gridcells with larger lai change 
            if len(dtdlai)>29:
              dlai70=N.percentile(dlai,70, axis=0)
              msk=MV.masked_less(dlai,dlai70)
              msk=msk/msk
              dtdlai=(dtdlai*msk).compressed()
              dlai=(dlai*msk).compressed()            
              #if len(dtdlai)>29:sensitivity[mth-1,yk,latk,lonk]=N.percentile(dtdlai,50, axis=0)
              if len(dtdlai)>3:
                sensitivity[yk,latk,lonk]=N.percentile(dtdlai,50, axis=0)
                STD[yk,latk,lonk]=(N.percentile(dtdlai,75, axis=0)-N.percentile(dtdlai,25, axis=0))/2.

  tabax=N.array(range(1,12+1,1),dtype='f')
  mtax = cdms.createAxis(tabax,id='months')
  mtax.long_name = 'months' 
  mtax.units = '-'

  tabax=N.array(range(1,66+1,1),dtype='f')
  ytax = cdms.createAxis(tabax,id='years')
  ytax.long_name = 'years' # Yi
  ytax.units = '-'
       
  f=cdms.open('/ESS_Datasets/EXT_ESM/CMIP6/Historical/SE/seamask.nc')
  land_mask=f('topo',lat=(-56,80),lon=(-180,180))
  f.close()
  latax=land_mask.getLatitude()
  lonax=land_mask.getLongitude()

  sensitivity=MV.masked_equal(sensitivity,0)
  STD=MV.masked_equal(STD,0)
  #sensitivity=cdms.createVariable(sensitivity,fill_value = 1.e+20,dtype='f',axes =[mtax,ytax,latax,lonax])
  sensitivity=cdms.createVariable(sensitivity,fill_value = 1.e+20,dtype='f',axes =[ytax,latax,lonax])
  STD=cdms.createVariable(STD,fill_value = 1.e+20,dtype='f',axes =[ytax,latax,lonax])
  f=cdms.open('sensitivity_%02i.nc'%mth,'w')
  f.write(sensitivity,id='dT_dLAI')
  f.write(STD,id='STD')
  f.close()

os.system("cdo ensmean sensitivity_??.nc sensitivity_yearly.nc")
