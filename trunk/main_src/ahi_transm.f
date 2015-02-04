	subroutine ahi_transm(ancil_data_path,temp,wvmr,ozmr,
     &             theta,clbl,kban, taut,*)
! * GOES-AHI: dry/wet/ozo transmittance at 101 AIRS levels
! .... version of 19.11.08
!
! .... modified 04-2014 Chia Moeller
!      adapted from ABI to AHI

! * LarrabeeStrow/HalWoolf/PaulVanDelst regression model based on
! *	LBLRTM line-by-line transmittances.
! * Input temperatures, and water-vapor and ozone mixing ratios, must
! *	be defined at the 101 levels in array 'pstd' (see block
! *   data 'reference_atmosphere' in file 'irtsubn101.f').
! * Logical unit numbers 11-15 are used for coefficient files.

! * Input
!	 temp = profile of temperature ........ degK
!	 wvmr = profile of H2O mixing ratio ... g/kg
!	 ozmr = profile of  O3 mixing ratio ... ppmv
!	 clbl = lblrtm version ................ 6[.12], 7[.04], 8[.30], or 11[.3]
!	theta = local zenith angle ............ deg
!	 kban = band or channel number ........ 7 - 16 

! * Output
!	 taut = profile of total transmittance (components are returned through common)
!	    * = error return in case of coefficient-file I/O trouble

	parameter (lfac=4,nk=5,nl=101,nm=nl-1,nr=10)
	parameter (nxd= 8,ncd=nxd+1,lencd=ncd*nm,lencdb=lencd*lfac)
	parameter (nxo= 9,nco=nxo+1,lenco=nco*nm,lencob=lenco*lfac)
	parameter (nxc= 4,ncc=nxc+1,lencc=ncc*nm,lenccb=lencc*lfac)
	parameter (nxl= 2,ncl=nxl+1,lencl=ncl*nm,lenclb=lencl*lfac)
	parameter (nxs=11,ncs=nxs+1,lencs=ncs*nm,lencsb=lencs*lfac)
	parameter (nxw=nxl+nxs)
	common/stdatm/pstd(nl),tstd(nl),wstd(nl),ostd(nl)
	common/taudwo/taud(nl),tauw(nl),tauo(nl)
	dimension temp(*),wvmr(*),ozmr(*),taut(*)
	dimension coefd(ncd,nm,nr),coefo(nco,nm,nr),coefc(ncc,nm,nr)
	dimension coefl(ncl,nm,nr),coefs(ncs,nm,nr)
	dimension pavg(nm),tref(nm),wref(nm),oref(nm)
	dimension tavg(nm),wamt(nm),oamt(nm),secz(nm)
	dimension tauc(nl),tlas(nl),wlas(nl),olas(nl)
	dimension xdry(nxd,nm),xozo(nxo,nm),xcon(nxc,nm),xwet(nxw,nm)
	character*32 cpath/'/home/chial/AHI/101_levels/coef/'/
	character*13 cfile/'ahixxx101.dat'/
        character*72 ancil_data_path,pfaast_path !akh
	character*3 comp(nk)/'dry','ozo','wco','wtl','wts'/
!	character*3 clblrtm(6:11)/'612','704','830','9xx','10x','113'/
	character*2 clbl
!	integer*4 iuc(nk)/11,12,13,14,15/
 	integer*4 iuc(nk)
	integer*4 lengcf(nk)/lencdb,lencob,lenccb,lenclb,lencsb/
	logical big_endian,newang,newatm
	data init/1/,tlas/nl*0./,wlas/nl*0./,olas/nl*0./,zlas/-999./
	secant(z)=1./cos(0.01745329*z)

        pfaast_path = trim(ancil_data_path)//"static/pfaast/" !akh

	if(init.ne.0) then
	   read(clbl,*) lbl
	   do l=1,nk
	      cfile(4:6)=comp(l)
	      lencf=lengcf(l)

           iuc(l)=get_lun() !akh
           open(iuc(l),file=trim(pfaast_path)//cfile,recl=lencf,
     +             access='direct', status='old',err=200)

	   enddo
	   do k=1,nr
  	      read(iuc(1),rec=k) ((coefd(i,j,k),i=1,ncd),j=1,nm)
  	      read(iuc(2),rec=k) ((coefo(i,j,k),i=1,nco),j=1,nm)
  	      read(iuc(3),rec=k) ((coefc(i,j,k),i=1,ncc),j=1,nm)
  	      read(iuc(4),rec=k) ((coefl(i,j,k),i=1,ncl),j=1,nm)
  	      read(iuc(5),rec=k) ((coefs(i,j,k),i=1,ncs),j=1,nm)
	   enddo
	   do l=1,nk
	      close(iuc(l))
	   enddo
	   if(big_endian()) then
	      call flip_rtc(coefd,ncd,nm,nr)
	      call flip_rtc(coefo,nco,nm,nr)
	      call flip_rtc(coefc,ncc,nm,nr)
	      call flip_rtc(coefl,ncl,nm,nr)
	      call flip_rtc(coefs,ncs,nm,nr)
	   endif
	   call conpir(pstd,tstd,wstd,ostd,nl,1,pavg,tref,wref,oref)
	   init=0
	endif

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

	k=kban-6

c * dry
	call taudoc(ncd,nxd,nm,coefd(1,1,k),xdry,taud)

c * ozo
	call taudoc(nco,nxo,nm,coefo(1,1,k),xozo,tauo)

c * wet
	call taudoc(ncc,nxc,nm,coefc(1,1,k),xcon,tauc)
	call tauwtr(ncs,ncl,nxs,nxl,nxw,nm,coefs(1,1,k),
     +		     coefl(1,1,k),xwet,tauw)
	do j=1,nl
	   tauw(j)=tauw(j)*tauc(j)
	enddo

c * total
	do j=1,nl
	   taut(j)=taud(j)*tauo(j)*tauw(j)
	enddo
	return
200	return1
	end
