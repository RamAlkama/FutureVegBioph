# -*- coding: latin-1 -*- 

import cdms2 as cdms, MV2 as MV, os.path
import numpy as N, cdutil, os, sys
from nearby import *
import os.path


cdms.setNetcdfShuffleFlag(0) 
cdms.setNetcdfDeflateFlag(0) 
cdms.setNetcdfDeflateLevelFlag(0)

# ----------------------------- #
 # --------------------------- #
  #       Parametres          #
 # --------------------------- #
# ----------------------------- #

modis='aqua' #  'terra'    'aqua'

nbzone = 5 # 5  16

latrange=(-56,80,'cc')# (-56,0,'cc') (0,80,'cc')
lonrange=(-180,180,'cc')

# ----------------------------- #
 # --------------------------- #
  #            def            #
 # --------------------------- #
# ----------------------------- #

if modis=='aqua' :
  tasfile='/ESS_EarthObs/DATA_PRODUCTS/AirTemp_fromLST/Tair_SG_%i.nc'

resolution=0.05 # 0.05 deg

if nbzone == 5 : zones=[1,2,3,4] # koppen-Geiger main climate zones "Boreal, Arid, Temperate and Tropical"

if latrange[0]<0: 
  lims=str(abs(latrange[0]))+'S'
else:
  lims=str(abs(latrange[0]))+'N'

if latrange[1]<0: 
  limn=str(abs(latrange[1]))+'S'
else:
  limn=str(abs(latrange[1]))+'N'


# ----------------------------- #
 # --------------------------- #
  #      functions            #
 # --------------------------- #
# ----------------------------- #


def read(filein,varin):
      fin=cdms.open(filein,'r')
      tab=fin(varin,lat=latrange,lon=lonrange)
      fin.close()
      return tab

def read2(filein,varin,mm_):
      fin=cdms.open(filein,'r')
      tab=fin(varin,cdms.timeslice(mm_-1,mm_),lat=latrange,lon=lonrange,squeeze=1)
      fin.close()
      return tab

def scatter(temp,treef2,treef1,surf,msk):

      tab    =N.zeros((nbzone,101,101))
      tabsurf=N.zeros((nbzone,101,101))
      TAS2D=N.zeros((temp.shape[0],temp.shape[1]))
      lati=temp.getLatitude()[:]
      loni=temp.getLongitude()[:]
      latt=N.multiply.outer(lati,loni/loni)
      lonn=N.multiply.outer(lati/lati,loni)
      lonn=cdms.createVariable(lonn,fill_value = 1.e+20,dtype='f',axes =surf.getAxisList())
      latt=cdms.createVariable(latt,fill_value = 1.e+20,dtype='f',axes =surf.getAxisList())

      for jj in zones :

        mms = MV.masked_not_equal(msk,jj)
        mms = mms/mms
        Ytmp= temp  * mms
        Xtmp2= treef2 * mms
        Xtmp1= treef1 * mms
        mms = (Ytmp*Xtmp1*Xtmp2+9999.)/(Ytmp*Xtmp1*Xtmp2+9999.)

        Ytmp= Ytmp*mms
        Xtmp1= Xtmp1*mms
        Xtmp2= Xtmp2*mms
        Stmp= surf* mms
        lon=lonn*mms
        lat=latt*mms

        # --------------------
        #  start fortran loop
        # --------------------

        Ytmp=Ytmp.compressed()
        Xtmp1=Xtmp1.compressed()
        Xtmp2=Xtmp2.compressed()
        Stmp=Stmp.compressed()
        lat=lat.compressed()
        lon=lon.compressed()

        if len(lon)>10:

          totnorm,Surf_step,tas2d=nearby.paris(Ytmp,Xtmp1,Xtmp2,Stmp,lat,lon,2,temp.shape[0],temp.shape[1],resolution,latrange[0],lonrange[0])

          tab[jj-1,:,:]=MV.where((Surf_step[:,:]+tabsurf[jj-1,:,:])!=0,(tab[jj-1,:,:]*tabsurf[jj-1,:,:]+totnorm[:,:]*Surf_step[:,:])/(Surf_step[:,:]+tabsurf[jj-1,:,:]),0)

          tabsurf[jj-1,:,:]=tabsurf[jj-1,:,:]+Surf_step[:,:]

        else:

          tas2d =N.zeros((temp.shape[0],temp.shape[1]))

        TAS2D=TAS2D[:,:]+tas2d[:,:]

        # ------------------
        #  end fortran loop
        # ------------------

      tab     = MV.masked_equal(tab    ,0)
      tabsurf = MV.masked_equal(tabsurf,0)

      return tab, tabsurf,TAS2D


tabax=N.array(range(1,101+1,1),dtype='f')
zonax = cdms.createAxis(tabax,id='to')
zonax.long_name = 'percent of tree'
zonax.units = '-'

zonax2 = cdms.createAxis(tabax,id='from')
zonax2.long_name = 'percent of tree'
zonax2.units = '-'

tabax3=N.array(range(1,nbzone+1,1),dtype='f')
zonax3 = cdms.createAxis(tabax3,id='climzonaxe')
zonax3.long_name = 'climate zone axe'
zonax3.units = '-'

zclim=read('/workdir/DATA/FOREST/koppen-Geiger_7200x3600_%s_masked.nc' %(str(nbzone)),'climzone')
zclim=MV.masked_equal(zclim,0)

surface = read('/workdir/DATA/FOREST/3Min/areacella.nc','areacella')

  # ----------------------------- #
   # --------------------------- #
    #          Calculs          #
   # --------------------------- #
  # ----------------------------- #

mask_lakes=read('/workdir/DATA/MODIS_Land_Cover/NETCDF/mask_terre_mer.nc','mask_water')
for yrdeb in range(2003,2014,1):#(2013,2002,-1):  
 for yf in range(yrdeb+1,2015,1):  
  print yrdeb,yf
 
  for mois in [1,2,3,4,5,6,7,8,9,10,11,12]:
    print mois

    if mois < 1 or mois > 12 :
      print 'mois %i must be between 1 and 12' %(mois)
      sys.exit()

    else:

      #ffracobs = cdms.open('/ESS_EarthObs/SATELLITE_RS/GLASS/LAI/0.05deg/yearly_GLASS_LAI_0.05_%i.nc'%(yrdeb),'r')
      #ftree_1    = ffracobs('LAI',lat=latrange,lon=lonrange,squeeze=1)
      ffracobs = cdms.open('/ESS_EarthObs/SATELLITE_RS/GLASS/LAI/0.05deg/monthly_GLASS_LAI_0.05_%i.nc'%(yrdeb),'r')
      ftree_1    = ffracobs('LAI',cdms.timeslice(mois-1,mois),lat=latrange,lon=lonrange,squeeze=1)
      ffracobs.close()
      ftree_1=(MV.where(ftree_1>6,6,ftree_1))*100/6.    # translate LAImax to 100 (LAImax=6m2/m2)

      #ffracobs = cdms.open('/ESS_EarthObs/SATELLITE_RS/GLASS/LAI/0.05deg/yearly_GLASS_LAI_0.05_%i.nc'%(yf),'r')
      #ftree_2    = ffracobs('LAI',lat=latrange,lon=lonrange,squeeze=1)
      ffracobs = cdms.open('/ESS_EarthObs/SATELLITE_RS/GLASS/LAI/0.05deg/monthly_GLASS_LAI_0.05_%i.nc'%(yf),'r')
      ftree_2    = ffracobs('LAI',cdms.timeslice(mois-1,mois),lat=latrange,lon=lonrange,squeeze=1)
      ffracobs.close()
      ftree_2=(MV.where(ftree_2>6,6,ftree_2))*100/6. 

      ftree_1=ftree_1*mask_lakes
      ftree_2=ftree_2*mask_lakes

      rep='/ESS_EarthObs/SATELLITE_RS/GLASS/LAI/AST/OUTPUT/%i-%i/'%(yrdeb,yf)
      os.system("mkdir -p %s"%(rep))
      
      for var in ['Tair']:
        if not (os.path.exists(rep+'%s_scater_%i%i%02i_%s_%s_%s.nc' %(var,yrdeb,yf,mois,modis,lims,limn))):
          fw=cdms.open(rep+'%s_scater_%i%i%02i_%s_%s_%s.nc' %(var,yrdeb,yf,mois,modis,lims,limn),'w')
          print var

          LSTd_2  = read2(tasfile %(yf),var,mois)
          LSTd_2  = MV.masked_greater(LSTd_2,400)

          LSTd_1  = read2(tasfile %(yrdeb),var,mois)
          LSTd_1  = MV.masked_greater(LSTd_1,400)

          dASTmax = cdms.createVariable(LSTd_2-LSTd_1,fill_value = 1.e+20,dtype='f',axes =ftree_1.getAxisList())

          del LSTd_2,LSTd_1

          tot,totsurf,TAS2D=scatter(dASTmax,ftree_1,ftree_2,surface,zclim)

          tot=MV.masked_greater(tot,50)
          tot=MV.masked_less(tot,-50)
          totsurf=MV.masked_greater(totsurf,1.e15)
          totsurf=MV.masked_less(totsurf,-1.e15)

          tot=cdms.createVariable(tot,fill_value = 1.e+20,dtype='f',axes =[zonax3,zonax,zonax2])
          totsurf=cdms.createVariable(totsurf,fill_value = 1.e+20,dtype='f',axes =[zonax3,zonax,zonax2])

          fw.write(tot,    id='dlta_tmp')
          fw.write(totsurf,id='area')
          TAS2D=N.ma.masked_invalid(TAS2D)
          TAS2D[TAS2D.mask]=0
          TAS2D=cdms.createVariable(TAS2D,fill_value = 1.e+20,dtype='f',axes =ftree_1.getAxisList())
          fw.write(TAS2D,id='Dtas')

          fw.close()

