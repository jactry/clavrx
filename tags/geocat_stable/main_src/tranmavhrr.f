!$Id$
c--------------------------------------------------------------------------
c Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
c
c NAME: tranmavhrr.f (src)
c
c PURPOSE: Calculate AVHRR transmittances
c
c DESCRIPTION: (see below)
c
c AUTHORS:
c  Andrew Heidinger, Andrew.Heidinger@noaa.gov
c  Andi Walther, CIMSS, andi.walther@ssec.wisc.edu
c  Denis Botambekov, CIMSS, denis.botambekov@ssec.wisc.edu
c  William Straka, CIMSS, wstraka@ssec.wisc.edu
c
c COPYRIGHT
c THIS SOFTWARE AND ITS DOCUMENTATION ARE CONSIDERED TO BE IN THE PUBLIC
c DOMAIN AND THUS ARE AVAILABLE FOR UNRESTRICTED PUBLIC USE. THEY ARE
c FURNISHED "AS IS." THE AUTHORS, THE UNITED STATES GOVERNMENT, ITS
c INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND AGENTS MAKE NO WARRANTY,
c EXPRESS OR IMPLIED, AS TO THE USEFULNESS OF THE SOFTWARE AND
c DOCUMENTATION FOR ANY PURPOSE. THEY ASSUME NO RESPONSIBILITY (1) FOR
c THE USE OF THE SOFTWARE AND DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL
c SUPPORT TO USERS.
c
c--------------------------------------------------------------------------
       subroutine tranmavhrr(ancil_data_path,temp,wvmr,ozmr,zena, 
     &                      craft,kban,taut,*)

! * Calculate AVHRR transmittances
! .... version of 09.02.09

! .... Andrew Heidinger 11-2012 removed hardwire of clast so should skip
!                               duplicate coefficient read when possible
!
! .... Chia Moeller 3-2012 removed flip

! ++++ for TIROS-N, NOAA-6 ... NOAA-19, METOP-A, ff.

! * LarrabeeStrow/HalWoolf/PaulVanDelst regression model based on
! *       LBLRTM line-by-line transmittances.
! * Input temperatures, and water-vapor and ozone mixing ratios, must
! *       be defined at the pressure levels in array 'pstd'
! *       (see block data 'reference_atmosphere').
! * Units: temperature, deg-K; water vapor, g/kg; ozone, ppmv.
! * Logical units 71-75 are used for coefficient files. (NO, GET_LUN IS USED - AKH)
! * Component tau's are returned through common, product in 'taut'.
 
! * Inputs:
!        temp = profile of temperature ............ degK
!        wvmr = profile of H2O mixing ratio ....... g/kg
!        ozmr = profile of  O3 mixing ratio ....... ppmv
!        zena = local zenith angle ................ deg
!        rco2 = CO2 mixing ratio .................. ppmv
!        craft = spacecraft: tirosn,noaa06...noaa19,metopa... upper or lower case
!        kban = band number ....................... 3,4,5

! * Outputs:
!              taut = total transmittance
!                 * = alternate return if any coefficient-file I/O problems
!          in common/taudwo/
!              taud = transmittance due to uniformly mixed gases
!              tauw = transmittance due to water vapor
!              tauo = transmittance due to ozone

      parameter (lfac=4,nk=5,nl=101,nm=nl-1,nr=3)
      parameter (nxc= 4,ncc=nxc+1,lencc=ncc*nm,lenccb=lencc*lfac)
      parameter (nxd= 8,ncd=nxd+1,lencd=ncd*nm,lencdb=lencd*lfac)
      parameter (nxo= 9,nco=nxo+1,lenco=nco*nm,lencob=lenco*lfac)
      parameter (nxl= 2,ncl=nxl+1,lencl=ncl*nm,lenclb=lencl*lfac)
      parameter (nxs=11,ncs=nxs+1,lencs=ncs*nm,lencsb=lencs*lfac)
      parameter (nxw=nxl+nxs)
      parameter (rco2=380.0)
      common/stdatm/pstd(nl),tstd(nl),wstd(nl),ostd(nl)
      common/taudwo/taud(nl),tauw(nl),tauo(nl)
      dimension temp(*),wvmr(*),ozmr(*),taut(*)
      dimension coefd(ncd,nm,nr),coefo(nco,nm,nr),coefl(ncl,nm,nr)
      dimension coefs(ncs,nm,nr),coefc(ncc,nm,nr), cbufs(lencs)
      dimension pavg(nm),tref(nm),wref(nm),oref(nm)
      dimension tavg(nm),wamt(nm),oamt(nm),secz(nm)
      dimension tauc(nl),tlas(nl),wlas(nl),olas(nl)
      dimension xdry(nxd,nm),xozo(nxo,nm),xwet(nxw,nm),xcon(nxc,nm)
      character*12 cfile/'avhrncom.dat'/
      character*72 ancil_data_path,pfaast_path
      character*6 craft
      character*6 clast/'999999'/
      character*3 comp(nk)/'dry','ozo','wco','wtl','wts'/
      integer*4 iuc(nk)
      integer*4 lencf(nk)/lencdb,lencob,lenccb,lenclb,lencsb/
      integer*4 get_lun
      logical flip(nk), newang,newatm
      data init/0/,tlas/nl*0./,wlas/nl*0./,olas/nl*0./,zlas/-999./
      secant(z)=1./cos(0.01745329*z)

      !--- this path is specific to CLAVR-x
      pfaast_path = trim(ancil_data_path)//"static/pfaast/"

       if(craft.ne.clast) then     !begin check clast if test
          do l=1,nk
             cfile(6:8)=comp(l)
             iuc(l)=GET_LUN()
             open(iuc(l),file=trim(pfaast_path)//cfile,
     +            recl=lencf(l),access='direct',
     +             status='old',err=200)
          enddo
          call getnumsc(craft,noff)
          koff=noff*(nr+1)+1
          krec=koff
! * determine if coefficients need to be byte-flipped
       do k=1,nk
            lenk=lencf(k)/lfac
            read(iuc(k),rec=krec) (cbufs(j),j=1,lenk)
            flip(k)=cbufs(1).lt.1.0
       enddo


! * read in coefficients
          do k=1,nr
             krec=k+koff
             read(iuc(1),rec=krec) ((coefd(i,j,k),i=1,ncd),j=1,nm)
             read(iuc(2),rec=krec) ((coefo(i,j,k),i=1,nco),j=1,nm)
             read(iuc(3),rec=krec) ((coefc(i,j,k),i=1,ncc),j=1,nm)
             read(iuc(4),rec=krec) ((coefl(i,j,k),i=1,ncl),j=1,nm)
             read(iuc(5),rec=krec) ((coefs(i,j,k),i=1,ncs),j=1,nm)
          enddo
          do l=1,nk
             close(iuc(l))
          enddo
! * flip if necessary
         if( flip(1) ) call flip_rtc(coefd,ncd,nm,nr)
         if( flip(2) ) call flip_rtc(coefo,nco,nm,nr)
         if( flip(3) ) call flip_rtc(coefc,ncc,nm,nr)
         if( flip(4) ) call flip_rtc(coefl,ncl,nm,nr)
         if( flip(5) ) call flip_rtc(coefs,ncs,nm,nr)
         call conpir(pstd,tstd,wstd,ostd,nl,1,pavg,tref,wref,oref)
         clast=craft
       endif !end check clast if test

       dt=0.
       dw=0.
       do=0.
       do j=1,nl
          dt=dt+abs(temp(j)-tlas(j))
          tlas(j)=temp(j)
          dw=dw+abs(wvmr(j)-wlas(j))
          wlas(j)=wvmr(j)
          do=do+abs(ozmr(j)-olas(j))
          olas(j)=ozmr(j)
          taud(j)=1.0
          tauw(j)=1.0
          tauc(j)=1.0
          tauo(j)=1.0
          taut(j)=1.0
       enddo
       datm=dt+dw+do
       newatm=datm.ne.0.
! * if atmosphere has changed, convert atmospheric parameters
       if(newatm) then
          call conpir(pstd,temp,wvmr,ozmr,nl,1,pavg,tavg,wamt,oamt)
       endif

! * if angle has changed, recalculate secant
       newang=zena.ne.zlas
       if(newang) then
          zsec=secant(zena)
          do l=1,nm
             secz(l)=zsec
          enddo
          zlas=zena
       endif

! * if atmosphere OR angle has changed, recalculate predictors
       if(newang.or.newatm) then
          call calpir(tref,wref,oref,tavg,wamt,oamt,pavg,secz,
     +                      nm,nxd,nxw,nxo,nxc,xdry,xwet,xozo,xcon)
       endif

       k=kban-2
! * dry
       call taudoc(ncd,nxd,nm,coefd(1,1,k),xdry,taud)
       if(rco2.ne.380.) then
! ... adjust for CO2 variation from model basis
          ratio=rco2/380.
          do j=1,nl
             if(taud(j).gt.0.0 .and. taud(j).lt.1.0) then
                taud(j)=taud(j)**ratio
             endif
          enddo
       endif
! * ozo
       call taudoc(nco,nxo,nm,coefo(1,1,k),xozo,tauo)
! * wet
       call taudoc(ncc,nxc,nm,coefc(1,1,k),xcon,tauc)
       call tauwtr(ncs,ncl,nxs,nxl,nxw,nm,coefs(1,1,k),
     +                     coefl(1,1,k),xwet,tauw)
       do j=1,nl
          tauw(j)=tauw(j)*tauc(j)
       enddo
! * total
       do j=1,nl
          taut(j)=taud(j)*tauo(j)*tauw(j)
       enddo
       return
200       return1
       end
