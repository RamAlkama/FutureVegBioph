!
!f2py -c --f90flags="-fopenmp -lm " --opt="-O3 -ffast-math" -lgomp --fcompiler=gnu95 nearby.f90 -m nearby
!
      MODULE nearby

      implicit none
      include "omp_lib.h"

      CONTAINS

      SUBROUTINE set_num_threads(n)
      !f2py threadsafe
      !f2py intent(in) n
      INTEGER :: n
      CALL OMP_SET_NUM_THREADS(n)
      END SUBROUTINE set_num_threads

!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine paris(tas,frac1,frac2,surf,lat,lon,n, &
                       ny,nx,reso,latmin,lonmin,nxy,   &
                       tas_out,srf_out,tas2d)

      implicit none
      
      integer,               intent(in ):: n, nxy,nx,ny
      real,                  intent(in ):: reso,latmin,lonmin
      real, dimension(nxy),  intent(in ):: tas
      real, dimension(nxy),  intent(in ):: frac1,frac2
      real, dimension(nxy),  intent(in ):: lat,lon,surf
      real, dimension(101,101),intent(out):: tas_out,srf_out
      real, dimension(ny,nx),intent(out):: tas2d
      real, dimension(:),   allocatable :: msk,srf_tmp,frac1lu,frac2lu
      real, dimension(:),   allocatable :: lu_dtas,lu_lat,lu_lon,msk2
      real, dimension(:),   allocatable :: nolu_dtas,nolu_lat,nolu_lon
      real    :: dist,wdi,wti, zz
      integer :: ii, kk, ll, jj, nlu,nnolu
!f2py intent(in) tas,frac1,frac2,surf,lat,lon,n
!f2py intent(in)  ny,nx,reso,latmin,lonmin
!f2py intent(out) tas_out,srf_out,tas2d
!f2py intent(hide),depend(tas) :: nxy=shape(tas,0)

      allocate(msk(nxy))
      allocate(msk2(nxy))

      msk       = 0.
      tas_out   = 0.
      srf_out   = 0.
      tas2d     = 0.

      ! delta : tas, frac

      where(abs(frac1(:)-frac2(:))<2)
        msk(:)=1.
      endwhere

      zz    =sum(msk(:))
      nnolu =int(zz)
      nlu   =nxy-nnolu

      ! landuse

      allocate(srf_tmp(nlu))
      srf_tmp=PACK(surf, mask=msk<0.5)

      allocate(lu_lat(nlu))
      lu_lat=PACK(lat, mask=msk<0.5)

      allocate(lu_lon(nlu))
      lu_lon=PACK(lon, mask=msk<0.5)

      allocate(lu_dtas(nlu))
      lu_dtas=PACK(tas(:), mask=msk<0.5)

      do ii=1,nlu
         wdi=0.
         wti=0.
         ! no landuse
         msk2= 0.
         where(abs(lat(:)-lu_lat(ii))<0.5 .and. &
               abs(lon(:)-lu_lon(ii))<1. .and. msk>0.5)
               msk2=1.
         endwhere
         nnolu=int(sum(msk2(:)))
         allocate(nolu_lat(nnolu))
         nolu_lat=PACK(lat, mask=msk2>0.5)

         allocate(nolu_lon(nnolu))
         nolu_lon=PACK(lon, mask=msk2>0.5)

         allocate(nolu_dtas(nnolu))
         nolu_dtas=PACK(tas(:), mask=msk2>0.5)

         do kk=1,nnolu
            call distance(lu_lat(ii),nolu_lat(kk),lu_lon(ii), &
                          nolu_lon(kk),dist)
            if (dist<50)then
               wdi=wdi+1/dist
               wti=wti+nolu_dtas(kk)/dist
            endif
         enddo

         deallocate(nolu_lat)
         deallocate(nolu_lon)
         deallocate(nolu_dtas)
         if (wdi /=0. )then
           lu_dtas(ii)=lu_dtas(ii)-wti/wdi
         else
           lu_dtas(ii)=0.
           srf_tmp(ii)=0.
         endif
         jj=nint((lu_lat(ii)-latmin-reso/2.)*1./reso)
         ll=nint((lu_lon(ii)-lonmin-reso/2.)*1./reso)
         tas2d(jj,ll)=lu_dtas(ii)
      enddo

      deallocate(msk2)
      deallocate(lu_lat)
      deallocate(lu_lon)

      allocate(frac1lu(nlu))
      frac1lu=PACK(frac1, mask=msk<0.5)

      allocate(frac2lu(nlu))
      frac2lu=PACK(frac2, mask=msk<0.5)
      deallocate(msk)
      !
      do ii=1,nlu
         ll=nint(frac1lu(ii))+1
         kk=nint(frac2lu(ii))+1
         tas_out(ll,kk)=tas_out(ll,kk)+lu_dtas(ii)*srf_tmp(ii)
         srf_out(ll,kk)=srf_out(ll,kk)+srf_tmp(ii)
      enddo
      !
      deallocate(srf_tmp)
      deallocate(frac1lu)
      deallocate(frac2lu)
      deallocate(lu_dtas)

      where(srf_out(:,:) /=0.)
        tas_out(:,:)=tas_out(:,:)/srf_out(:,:)
      endwhere
      
      end subroutine paris

!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine tab1d_to_tab2d(nx,ny,reso,latmin,lonmin,lu_lat,lu_lon,&
                                lu_dtas,nlu,tas2d)

      integer,               intent(in ):: nlu,nx,ny
      real,                  intent(in ):: reso,latmin,lonmin
      real, dimension(nlu),  intent(in ):: lu_lat,lu_lon,lu_dtas
      real, dimension(ny,nx),intent(out):: tas2d
      integer :: ii,jj,ll
!f2py intent(in)  nx,ny,reso,latmin,lonmin,lu_lat,lu_lon,lu_dtas
!f2py intent(out) tas2d
!f2py intent(hide),depend(lu_dtas) :: nlu=shape(lu_dtas,0)

      tas2d=0.
      do ii=1,nlu
         jj=nint((lu_lat(ii)-latmin-reso/2.)*1./reso)
         ll=nint((lu_lon(ii)-lonmin-reso/2.)*1./reso)
         tas2d(jj,ll)=lu_dtas(ii)
      enddo

      end subroutine tab1d_to_tab2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine distance(lata,latb,lona,lonb,dist)
      REAL, INTENT(in)   :: lata,latb,lona,lonb
      REAL, INTENT(out)  :: dist
      REAL :: R, pi, a, b, c, d
      pi=3.14159265359
      R=6367.445
      a=(pi*lata)/180.
      b=(pi*latb)/180.
      c=(pi*lona)/180.
      d=(pi*lonb)/180.   
      dist=R*acos( sin(a)*sin(b)+cos(a)*cos(b)*cos(c-d) )
      end subroutine distance

!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module

!!!!!!!!!!!!!!!!!!!!!!!!!!!

