c--------------------------------------------------------------------------
c Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
c
c NAME: tran_viirsm.f (src)
c
c PURPOSE: 
c * VIIRS -- calculate transmittances
c
c DESCRIPTION:  (see below)
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
      subroutine tran_viirsm(ancil_data_path,temp,wvmr,ozmr,theta, 
     &                      rco2,kban,taut,*)
c * VIIRS -- calculate transmittances
c .... version of 12.09.01
c
c ... modified by Chia Moeller 11-2011

! * LarrabeeStrow/HalWoolf/PaulVanDelst regression model based on
! *     LBLRTM line-by-line transmittances.
! * Input temperatures, and water-vapor and ozone mixing ratios, must
! *     be defined at the 101 pressure levels in array 'pstd' (see block
! *   data 'reference_atmosphere' in file 'irtsubn101.f').
! * Units: temperature, deg-K; water vapor, g/kg; ozone, ppmv.
! * Logical unit numbers 71-75 are used for coefficient files.

! * Input
!        ancil_data_path = location of PFAAST viirscom.dat files
!        temp = profile of temperature ............ degK
!        wvmr = profile of H2O mixing ratio ....... g/kg
!        ozmr = profile of  O3 mixing ratio ....... ppmv
!       theta = local zenith angle ................ deg
!        rco2 = CO2 mixing ratio .................. ppmv
!       craft = spacecraft, upper or lower case ... tirosn,noaa06,...,noaa18,metopa
!         kan = channel number .................... 1 - 16

! * Output
!        taut = profile of total transmittance (components are returned through common)
!           * = error return in case of coefficient-file I/O trouble

! WCS edits - Added pfaast_path, GET_LUN calls.

c ... parameters ..
c     nk = number of comp
c     nr = number of bands

      parameter (nk=5,nl=101,nm=nl-1,nr=5)
      parameter (nxd= 8,ncd=nxd+1,lencd=ncd*nm,lencdb=lencd*4)
      parameter (nxo= 9,nco=nxo+1,lenco=nco*nm,lencob=lenco*4)
      parameter (nxc= 4,ncc=nxc+1,lencc=ncc*nm,lenccb=lencc*4)
      parameter (nxl= 2,ncl=nxl+1,lencl=ncl*nm,lenclb=lencl*4)
      parameter (nxs=11,ncs=nxs+1,lencs=ncs*nm,lencsb=lencs*4)
      parameter (nxw=nxl+nxs)
      common /stdatm/pstd(nl),tstd(nl),wstd(nl),ostd(nl)
      common /taudwo/taud(nl),tauw(nl),tauo(nl)
      dimension temp(*),wvmr(*),ozmr(*),taut(*) 
      dimension coefd(ncd,nm,nr),coefo(nco,nm,nr),coefc(ncc,nm,nr)
      dimension coefl(ncl,nm,nr),coefs(ncs,nm,nr),lencf(nk),iuc(nk)
      dimension cbufs(lencs)
      dimension pavg(nm),tref(nm),wref(nm),oref(nm)
      dimension tavg(nm),wamt(nm),oamt(nm),secz(nm)
      dimension tauc(nl),tlas(nl),wlas(nl),olas(nl)
      dimension xdry(nxd,nm),xozo(nxo,nm),xwet(nxw,nm),xcon(nxc,nm)

      character*12 cfile
      character*72 ancil_data_path,pfaast_path !WCS3
      character*3 comp(nk)
      logical newang,newatm, flip(nk)
      data olas/nl*0./,tlas/nl*0./,wlas/nl*0./,zlas/-999./
      data lencf/lencdb,lencob,lenccb,lenclb,lencsb/
      data cfile/'viirscom.dat'/
      data comp/'dry','ozo','wco','wtl','wts'/
      integer*4 get_lun   !needed by Get_lun() call by AKH     

      secant(z)=1./cos(0.01745329*z)

      pfaast_path = trim(ancil_data_path)//"static/pfaast/" !WCS3
 
        koff=0
! * define and open the coefficient files
        iux=70
        do l=1,nk
           cfile(6:8)=comp(l)
           iux=iux+1
           iux=get_lun() !WCS3
           open(iux,file=trim(pfaast_path)//cfile,recl=lencf(l),
     +             access='direct', status='old',err=200)
           iuc(l)=iux
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

        call conpir(pstd,tstd,wstd,ostd,nl,1,pavg,tref,wref,oref)
       
	dt=0.
	dw=0.
	do=0.
	do j=1,nl
	   taud(j)=1.0
	   tauw(j)=1.0
	   tauc(j)=1.0
	   tauo(j)=1.0
	   taut(j)=1.0
	   dt=dt+abs(temp(j)-tlas(j))
	   tlas(j)=temp(j)
	   dw=dw+abs(wvmr(j)-wlas(j))
	   wlas(j)=wvmr(j)
	   do=do+abs(ozmr(j)-olas(j))
	   olas(j)=ozmr(j)
	enddo

	datm=dt+dw+do
	newatm=datm.ne.0.

	if(newatm) then
	   call conpir(pstd,temp,wvmr,ozmr,nl,1,pavg,tavg,wamt,oamt)
	endif

	newang=theta.ne.zlas
	if(newang) then
	   zsec=secant(theta)
	   do l=1,nm
	      secz(l)=zsec
	   enddo
	   zlas=theta
	endif

	if(newang.or.newatm) then
	   call calpir(tref,wref,oref,tavg,wamt,oamt,pavg,secz,
     +			 nm,nxd,nxw,nxo,nxc,xdry,xwet,xozo,xcon)
	endif

	k=kban
c * dry
	call taudoc(ncd,nxd,nm,coefd(1,1,k),xdry,taud)

        if(rco2.ne.380.) then
! .... adjust for CO2 variation from model basis
           ratio=rco2/380.
           do j=1,nl
              if(taud(j).gt.0.0 .and. taud(j).lt.1.0) then
                 taud(j)=taud(j)**ratio
              endif
           enddo
        endif

c * ozo
	call taudoc(nco,nxo,nm,coefo(1,1,k),xozo,tauo)
c * wet
c ...continuum
	call taudoc(ncc,nxc,nm,coefc(1,1,k),xcon,tauc)
c ...without continuum
	call tauwtr(ncs,ncl,nxs,nxl,nxw,nm,coefs(1,1,k),
     +		     coefl(1,1,k),xwet,tauw)
c ...total water vapor
	do j=1,nl
	   tauw(j)=tauc(j)*tauw(j)
	enddo
c * total transmittance
	do j=1,nl
	   taut(j)=taud(j)*tauo(j)*tauw(j)
	enddo
	return
200	return1
	end
