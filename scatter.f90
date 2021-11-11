!
!f2py -c --f90flags="-fopenmp -lm " --opt="-O3 -ffast-math" -lgomp --fcompiler=gnu95 scatter.f90 -m scatt
!

      subroutine doit(dtas,dbiom,biom,nxy,dout,mnbr)

      implicit none
      
      integer,              intent(in):: nxy
      real, dimension(nxy), intent(in):: dtas,dbiom,biom
      real, dimension(40,151),intent(out):: dout
      real, dimension(40,151),intent(out):: mnbr
      integer :: ii, jj, ij

!f2py intent(in)  dtas,dbiom,biom
!f2py intent(out) dout
!f2py intent(hide),depend(tas) :: nxy=shape(dtas,0)
  
      dout   = 0.
      mnbr   = 0.

      do ij=1,nxy

        ii=nint(dbiom(ij)+75)+1
        if (ii<1)ii=1
        if (ii>151)ii=151
        jj=nint(biom(ij)/4.)
        if (jj<1)jj=1
        if (jj>40)jj=40
        dout(jj,ii)=dout(jj,ii)+dtas(ij)
        mnbr(jj,ii)=mnbr(jj,ii)+1.

      enddo

      where(mnbr(:,:) /=0.)
        dout(:,:)=dout(:,:)/mnbr(:,:)
      endwhere

      end subroutine doit
!--------------------------------------------------------
!  save data by climate (temp,prec) class 
!--------------------------------------------------------
      subroutine climate(tas,prec,ari,biom,nxy,dout,mnbr)

      implicit none
      
      integer,              intent(in):: nxy
      real, dimension(nxy), intent(in):: tas,prec,biom,ari
      real, dimension(50,50),intent(out):: dout,mnbr
      integer :: ii, jj, ij

!f2py intent(in)  tas,prec,biom,ari
!f2py intent(out) dout,mnbr
!f2py intent(hide),depend(tas) :: nxy=shape(tas,0)
  
      dout   = 0.
      mnbr   = 0.

      do ij=1,nxy

        ii=nint(prec(ij)/50)+1
        jj=nint(tas(ij)-253)
        if (jj<1)jj=1
        if (jj>50)jj=50
        if (ii<1)ii=1
        if (ii>50)ii=50
        dout(jj,ii)=dout(jj,ii)+biom(ij)*ari(ij)
        mnbr(jj,ii)=mnbr(jj,ii)+ari(ij)

      enddo

      where(mnbr(:,:) /=0.)
        dout(:,:)=dout(:,:)/mnbr(:,:)
      endwhere

      end subroutine climate
!--------------------------------------------------------
!  save data by climate (temp,prec) class 
!--------------------------------------------------------
      subroutine yeartoyear(yr1,yr2,delta,ari,nxy,dout,mnbr)

      implicit none
      
      integer,              intent(in):: nxy
      real, dimension(nxy), intent(in):: yr1,yr2,delta,ari
      real, dimension(50,50),intent(out):: dout,mnbr
      integer :: ii, jj, ij

!f2py intent(in)  yr1,yr2,delta,ari
!f2py intent(out) dout,mnbr
!f2py intent(hide),depend(yr1) :: nxy=shape(yr1,0)
  
      dout   = 0.
      mnbr   = 0.

      do ij=1,nxy

        ii=nint((yr1(ij)-18.)/2)+1
        jj=nint((yr2(ij)-18.)/2)+1
        if (jj<1)jj=1
        if (jj>50)jj=50
        if (ii<1)ii=1
        if (ii>50)ii=50
        dout(jj,ii)=dout(jj,ii)+delta(ij)*ari(ij)
        mnbr(jj,ii)=mnbr(jj,ii)+ari(ij)

      enddo

      where(mnbr(:,:) /=0.)
        dout(:,:)=dout(:,:)/mnbr(:,:)
      endwhere

      end subroutine yeartoyear
!--------------------------------------------------------
!  save data by climate (temp,prec) class 
!--------------------------------------------------------
      subroutine ytoy(yr1,yr2,delta,ari,nxy,dout,mnbr)

      implicit none
      
      integer,              intent(in):: nxy
      real, dimension(nxy), intent(in):: yr1,yr2,delta,ari
      real, dimension(50,50),intent(out):: dout,mnbr
      integer :: ii, jj, ij

!f2py intent(in)  yr1,yr2,delta,ari
!f2py intent(out) dout,mnbr
!f2py intent(hide),depend(yr1) :: nxy=shape(yr1,0)
  
      dout   = 0.
      mnbr   = 0.

      do ij=1,nxy

        ii=nint((yr1(ij))/2)+1
        jj=nint((yr2(ij))/2)+1
        if (jj<1)jj=1
        if (jj>50)jj=50
        if (ii<1)ii=1
        if (ii>50)ii=50
        dout(jj,ii)=dout(jj,ii)+delta(ij)*ari(ij)
        mnbr(jj,ii)=mnbr(jj,ii)+ari(ij)

      enddo

      where(mnbr(:,:) /=0.)
        dout(:,:)=dout(:,:)/mnbr(:,:)
      endwhere

      end subroutine ytoy

!--------------------------------------------------------
!  dT/dLAI as function of snow cover and soil moisture 
!--------------------------------------------------------
      subroutine dtdlai(mrso,snc,sncfunc,mrsofunc,nt,ny,nx,n,dtdl)

      implicit none
      
      integer,              intent(in):: n,nx,ny,nt
      real, dimension(nt,ny,nx), intent(in):: mrso,snc
      real, dimension(n), intent(in):: sncfunc,mrsofunc
      real, dimension(nt,ny,nx),intent(out):: dtdl
      integer :: t, ii, jj

!f2py intent(in)  mrso,snc,sncfunc,mrsofunc
!f2py intent(out) dtdl
!f2py intent(hide),depend(snc) :: nt=shape(snc,0),ny=shape(snc,1)
!f2py intent(hide),depend(snc) :: nx=shape(snc,2)
!f2py intent(hide),depend(snc_func) :: n=shape(sncfunc,0)
      dtdl   = 1.e20
      do t=1,nt
       do jj=1,ny
        do ii=1,nx
         if (snc(t,jj,ii)>=2 .and. snc(t,jj,ii)<=100)then
           dtdl(t,jj,ii)=sncfunc(int(snc(t,jj,ii)))+&
           (sncfunc(int(snc(t,jj,ii))+1)-sncfunc(int(snc(t,jj,ii))))*&
           (snc(t,jj,ii)-int(snc(t,jj,ii)))
         elseif (mrso(t,jj,ii)>=0 .and. mrso(t,jj,ii)<=100) then
           dtdl(t,jj,ii)=mrsofunc(int(mrso(t,jj,ii)))+&
           (mrsofunc(int(mrso(t,jj,ii))+1)-mrsofunc(int(mrso(t,jj,ii)))&
           )*(mrso(t,jj,ii)-int(mrso(t,jj,ii)))
         endif
        enddo
       enddo
      enddo
      !
      end subroutine dtdlai
!-------------------------------------
!  dT/dLAI from high to low resolution using median  
!-------------------------------------
      subroutine dtdlailow(dtdlai,lath,lonh,land,latl,lonl,ny,nx,n,&
                           dtdlmedian)
      implicit none
      
      integer,              intent(in):: n,nx,ny
      real, dimension(ny,nx), intent(in):: land
      real, dimension(ny), intent(in)::latl
      real, dimension(nx), intent(in)::lonl
      real, dimension(n), intent(in):: dtdlai,lath,lonh
      real, dimension(ny,nx),intent(out)::dtdlmedian
      real, dimension(:),   allocatable :: msk,msk2,a
      integer :: ii, jj,nlu,i,passo
      real :: temp
!f2py intent(in)  dtdlai,lath,lonh,land,latl,lonl
!f2py intent(out) dtdlmedian
!f2py intent(hide),depend(land) :: ny=shape(land,0),nx=shape(land,1)
!f2py intent(hide),depend(dtdlai) :: n=shape(dtdlai,0)
      dtdlmedian   = 1.e20
      allocate(msk(n))
      allocate(msk2(n))
      do jj=1,ny
        do ii=1,nx
          if (land(jj,ii)>0)then
           msk       = 1.
           msk2      = 1.
           where((lath(:)<latl(jj)-1.5) .or. (lath(:)>latl(jj)+1.5))
             msk(:)=0.
           endwhere
           where((lonh(:)<lonl(ii)-4.) .or. (lonh(:)>lonl(ii)+4.))
             msk2(:)=0.
           endwhere
           msk=msk*msk2
           nlu=int(sum(msk(:)))
           if (nlu>19)then
             allocate(a(nlu))
             a=PACK(dtdlai, mask=msk<0.5)
             ! first sort the data
             do passo=1,nlu 
               do i = 1,nlu
                 if(a(i) .gt. a(i+1)) then
                   temp = a(i)
                   a(i) = a(i+1)
                   a(i+1) = temp
                 endif
               enddo
             enddo
             ! second median 
             if ((nlu/2)*2 .EQ. nlu)then
               dtdlmedian(jj,ii)=(a(nlu/2)+a(nlu/2+1))/2.
             else
               dtdlmedian(jj,ii)=a(nlu/2+1)
             endif
             deallocate(a)
           endif
          endif
        enddo
      enddo
      deallocate(msk)
      deallocate(msk2)
      !
      end subroutine dtdlailow
!-------------------------------------
!  dT/dLAI from high to low resolution using linear lest squares
!-------------------------------------
      subroutine dtdlaiLLS(dt,dlai,lath,lonh,land,latl,lonl,ny,nx,n,&
                           dtdlLLS,dtdlLLS_inter,dtdlLLS_std)
      implicit none
      
      integer,              intent(in):: n,nx,ny
      real, dimension(ny,nx), intent(in):: land
      real, dimension(ny), intent(in)::latl
      real, dimension(nx), intent(in)::lonl
      real, dimension(n), intent(in):: dt,dlai,lath,lonh
      real, dimension(ny,nx),intent(out)::dtdlLLS,dtdlLLS_std
      real, dimension(ny,nx),intent(out):: dtdlLLS_inter
      real, dimension(:),   allocatable :: msk,msk2,X,Y
      integer :: ii, jj,nlu
      real :: a,b,d
!f2py intent(in)  dt,dlai,lath,lonh,land,latl,lonl
!f2py intent(out) dtdlLLS,dtdlLLS_std,dtdlLLS_inter
!f2py intent(hide),depend(land) :: ny=shape(land,0),nx=shape(land,1)
!f2py intent(hide),depend(dt) :: n=shape(dt,0)
      dtdlLLS    = 1.e20
      dtdlLLS_std= 1.e20
      dtdlLLS_inter=1.e20
      allocate(msk(n))
      allocate(msk2(n))
      do jj=1,ny
        do ii=1,nx
          if (land(jj,ii)>0)then
           msk(:)       = 1.
           msk2(:)      = 1.
           where((lath(:)<latl(jj)-3.) .or. (lath(:)>latl(jj)+3.))
             msk(:)=0.
           endwhere
           where((lonh(:)<lonl(ii)-4.) .or. (lonh(:)>lonl(ii)+4.))
             msk2(:)=0.
           endwhere
           msk=msk(:)*msk2(:)
           nlu=int(sum(msk(:)))
           if (nlu>29)then
             allocate(X(nlu))
             allocate(Y(nlu))
             Y=PACK(dt, mask=msk<0.5)
             X=PACK(dlai, mask=msk<0.5)
             call Least_Square(nlu,X,Y,a,b,d)
             ! first sort the data
             ! second median 
             dtdlLLS(jj,ii)=b
             dtdlLLS_std(jj,ii)=d
             dtdlLLS_inter(jj,ii)=a
             deallocate(X)
             deallocate(Y)
             !
           endif
          endif
        enddo
      enddo
      deallocate(msk)
      deallocate(msk2)
      !
      end subroutine dtdlaiLLS

!!!!!!!!!!!!!

!!!!!!!!!!!!!
      subroutine bubble_sort(datato,counto)
      integer,              intent(in):: counto
      real, dimension(counto), intent(inout):: datato    
      integer :: i, passo, sorted 
      real :: temp
      passo = 1
      sorted = 0
      do while(sorted .eq. 0) 
        sorted = 1
        do 2 i = 1,counto-passo
          if(datato(i) .gt. datato(i+1)) then
            temp = datato(i)
            datato(i) = datato(i+1)
            datato(i+1) = temp
            sorted = 0
          endif
 2      continue
        passo = passo +1
      end do
      return
      end subroutine bubble_sort
!!!!!!!!!!!!!
      subroutine ram_sort(datato,counto)
      integer,              intent(in):: counto
      real, dimension(counto), intent(inout):: datato    
      integer :: i, passo
      real :: temp
      do passo=1,counto 
        do i = 1,counto
          if(datato(i) .gt. datato(i+1)) then
            temp = datato(i)
            datato(i) = datato(i+1)
            datato(i+1) = temp
          endif
        enddo
      enddo
      end subroutine ram_sort
!!!!!!!!!!!!!
!***************************************************
!*        Linear least squares subroutine          *
!* ----------------------------------------------- *
!* The input data set is X(m), Y(m). The number of *
!* data points is n (n must be > 2). The returned  *
!* parameters are: a,b, coefficients of equation   *
!* Y = a + b X, and d, standard deviation of fit.  *
!***************************************************
      SUBROUTINE Least_Square(n,X,Y,a,b,d)
      integer,              intent(in):: n
      real, dimension(n), intent(in):: X,Y   
      real, intent(out):: a,b,d
      integer :: m
      real :: a1,a2,b0,b1,d1
      a1 = 0
      a2 = 0
      b0 = 0
      b1 = 0
      DO m = 1, n
        a1 = a1 + X(m)
        a2 = a2 + X(m) * X(m)
        b0 = b0 + Y(m)
        b1 = b1 + Y(m) * X(m)
      END DO
      a1 = a1 / n
      a2 = a2 / n
      b0 = b0 / n
      b1 = b1 / n
      d = a1 * a1 - a2
      a = a1 * b1 - a2 * b0
      a = a / d
      b = a1 * b0 - b1
      b = b / d
      !  Evaluation of standard deviation d (unbiased estimate) 
      d = 0
      DO m = 1, n
        d1 = Y(m) - a - b * X(m)
        d  = d + d1 * d1
      END DO
      d = SQRT(d / (n - 2))
      RETURN
      END subroutine Least_Square
