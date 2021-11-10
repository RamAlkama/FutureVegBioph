"""

 dT/dlai as function of evaporation and solar radiation 

"""

# ----------------------------------- #
 # --------------------------------- #
  #         import les moduls       #
 # --------------------------------- #
# ----------------------------------- #

from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import numpy as np
import matplotlib.pyplot as plt
import cdms2 as cdms, numpy as N
import MV2 as MV, cdutil
#from mpl_toolkits.basemap import Basemap, shiftgrid, cm
import numpy.polynomial as poly
import gdal,glob
import scatt

cdms.setNetcdfShuffleFlag(0)
cdms.setNetcdfDeflateFlag(0)
cdms.setNetcdfDeflateLevelFlag(0)

latrange=(0,80,'cc') 
lonrange=(-180,180,'cc')

def read_snowcoverLAI(yr_,mth_):
  f=cdms.open('/ESS_EarthObs/REANALYSIS/ERA5/radiative/SWdown_005deg_%i.nc'%(yr_))
  temp=f('ssrd',cdms.timeslice(mth_-1,mth_),lat=latrange,lon=lonrange,squeeze=1)/86400. # J/m2 to W/m2
  f.close()
  msk=(temp+999)/(temp+999)
  temp=MV.where(temp<0,0,temp)
  temp=MV.where(temp>460.,460.,temp) # 460 W/m2 is the maximum obseerved SWdown
  temp=100.*((temp-0.)/(460.-0.))*msk
  #read hdf_snowcover
  filename=glob.glob('/ESS_EarthObs/SATELLITE_RS/MODIS_ARC/MODIS/MYD10CM.006/%i.%02i.01/MYD10CM.*.hdf' %(yr_,mth_))[0]
  hdf_file = gdal.Open(filename) 
  subDatasets = hdf_file.GetSubDatasets()
  subDatasets
  snowcover = gdal.Open(subDatasets[0][0])
  snowcover=snowcover.ReadAsArray()
  snowcover=snowcover[int(10/0.05):-int(90/0.05),:] 
  snowcover=MV.masked_greater(snowcover[::-1,:],1)#mask snowcover >1%
  snowcover=(snowcover+999)/(snowcover+999)
  # read latent heat
  f = cdms.open('/ESS_EarthObs/DATA_PRODUCTS/GLEAM/GLEAM3a/TMP/E_%i_GLEAM_v3.1a_005deg_monthly.nc'%yr_,'r')
  sm  = f('E',cdms.timeslice(mth_-1,mth_),lat=latrange,lon=lonrange,squeeze=1)*snowcover
  f.close()
  del snowcover
  # ERA5
  msk=msk*(sm+999)/(sm+999)
  sm=MV.where(sm<0,0,sm)
  sm=MV.where(sm>7.,7.,sm)    # 7 mm/day  is the maximum obseerved evaporation 
  sm=100.*((sm-0.)/(7.-0.))*msk
  # read LAI
  f = cdms.open('/ESS_EarthObs/SATELLITE_RS/GLASS/LAI/0.05deg/monthly_GLASS_LAI_0.05_%i.nc'%(yr_),'r')
  lai_  = f('LAI',cdms.timeslice(mth_-1,mth_),lat=latrange,lon=lonrange,squeeze=1)
  f.close()
  lai_=MV.masked_equal(lai_,0)
  return sm,lai_,temp,msk


def andiff(yd,yf):
  for mm in range(1,13,1):
    print 'month=',mm
    f=cdms.open('/ESS_EarthObs/SATELLITE_RS/GLASS/LAI/AST/OUTPUT/%i-%i/Tair_2D_%i%i_56S_80N.nc' %(yd,yf,yd,yf))
    dt=f('Dtas',cdms.timeslice(mm-1,mm),lat=latrange,lon=lonrange,squeeze=1)
    f.close()
    latent1,lai1,temp1,msk1=read_snowcoverLAI(yd,mm)
    latent2,lai2,temp2,msk2=read_snowcoverLAI(yf,mm)
    msk=msk1*msk2
    del msk1,msk2
    abstemp=(temp1+temp2)/2.
    latent=(latent1+latent2)/2.
    dlai=N.log(MV.masked_equal(lai2/lai1,1)) 
    del lai1,lai2,temp1,temp2,latent1,latent2
    dt_dlai=MV.masked_less(dt/dlai,-20)   # normaly never happend but in case of strange numbers due to the logarithmic equation
    dt_dlai=MV.masked_greater(dt_dlai,20) # normaly never happend but in case of strange numbers due to the logarithmic equation
    del dlai, dt
    ari=area[:,:]
    msk=msk*dt_dlai*abstemp*latent*ari
    msk=(msk+999)/(msk+999)
    ari=ari*msk
    abstemp=abstemp*msk
    dt_dlai=dt_dlai*msk
    latent=latent*msk
    sesitivity_tmp,surface_tmp=scatt.ytoy(latent.compressed(),abstemp.compressed(),dt_dlai.compressed(),ari.compressed())
    if mm ==1:
      Ysesitivity,Ysurface=sesitivity_tmp*surface_tmp,surface_tmp
    else:
      Ysesitivity=Ysesitivity+sesitivity_tmp*surface_tmp
      Ysurface=Ysurface+surface_tmp

  Ysesitivity=MV.where(Ysurface>0,Ysesitivity/Ysurface,0)
  return Ysesitivity,Ysurface


f=cdms.open('/workdir/DATA/FOREST/3Min/areacella.nc')
area = f('areacella',lat=latrange,lon=lonrange)/1.e6 # m2 to km2
SENS,AREA=N.zeros((66,50,50)),N.zeros((66,50,50))
jj=-1
for yd in range(2003,2014,1):  
  for yf in range(yd+1,2015,1): 
     jj=jj+1
     print yd,yf
     sesitivity_tmp,surface_tmp=andiff(yd,yf)
     AREA[jj,:,:]=surface_tmp[:,:]
     SENS[jj,:,:]=sesitivity_tmp[:,:]   

sesitivity=N.percentile(SENS,50,axis=0)

tabax=N.array(range(1,50+1,1),dtype='f')
zonax = cdms.createAxis(tabax,id='to')
zonax.long_name = 'SWdown' # Yi+1
zonax.units = '-'

zonax2 = cdms.createAxis(tabax,id='from')
zonax2.long_name = 'Evap' # Yi
zonax2.units = '-'

tabax3=N.array(range(1,66+1,1),dtype='f')
zonax3 = cdms.createAxis(tabax3,id='years')
zonax3.long_name = 'couple of years'
zonax3.units = '-'

tabax4=N.array(range(1,21,1),dtype='f')
zonax4 = cdms.createAxis(tabax4,id='temperature')
zonax4.long_name = 'temperature range from 0 to 40 with 2 deg step'
zonax4.units = 'deg C'

SENS=MV.masked_greater(SENS,1.e19)
AREA=MV.masked_greater(AREA,1.e19)
sesitivity=MV.masked_greater(sesitivity,1.e19)
SENS=cdms.createVariable(SENS,fill_value = 1.e+20,dtype='f',axes =[zonax3,zonax,zonax2])
AREA=cdms.createVariable(AREA,fill_value = 1.e+20,dtype='f',axes =[zonax3,zonax,zonax2])
sesitivity=cdms.createVariable(sesitivity,fill_value = 1.e+20,dtype='f',axes =[zonax,zonax2])
f=cdms.open('Evap_SWdown_sensitivity_logLAI.nc','w')
f.write(AREA,id='area')
f.write(sesitivity,id='meandT_dLai')
f.close()


