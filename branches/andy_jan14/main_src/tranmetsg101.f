c--------------------------------------------------------------------------
c Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
c
c NAME: tranmetsg101.f (src)
c
c PURPOSE: 
c * METEOSAT SECOND GENERATION transmittance calculation
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
        subroutine tranmetsg101(ancil_data_path,temp,wvmr,ozmr, zena, 
     &                        mets,kban,taut,*)

! RCS Version: $Id:

c * METEOSAT SECOND GENERATION transmittance calculation
c .... version of 04.05.06

c * LarrabeeStrow/HalWoolf/PaulVanDelst regression model based on
c *     LBLRTM line-by-line transmittances.
c * Input temperatures, and water-vapor and ozone mixing ratios, must
c *     must be defined at the 101 pressure levels in array 'pstd'
c *     (see block data 'reference_atmosphere').
c * Units: temperature, deg-K; water vapor, g/kg; ozone, ppmv.
c * Logical units 71-75 are used for coefficient files.
c * Component tau's are returned through common, product in 'taut'.

c * Inputs:
c               temp = temperature profile (degK)
c               wvmr = water-vapor mixing-ratio profile (g/kg)
c               ozmr = ozone mixing-ratio profile (ppmv)
c               zena = local zenith angle in degrees
c               mets = Meteosat number (8...)
c               kban = band number: 4...11 corresponding to
c                          3.9,6.2,7.3,8.7,9.7,10.8,12.0,13.4 microns
c * Outputs:
c               taut = total transmittance
c                  * = alternate return if any coefficient-file I/O problems
c          in common/taudwo/
c               taud = transmittance due to uniformly mixed gases
c               tauw = transmittance due to water vapor
c               tauo = transmittance due to ozone

        parameter (lfac=4,nk=5,nl=101,nm=nl-1,nr=8)
        parameter (nxc= 4,ncc=nxc+1,lencc=ncc*nm,lenccb=lencc*lfac)
        parameter (nxd= 8,ncd=nxd+1,lencd=ncd*nm,lencdb=lencd*lfac)
        parameter (nxo= 9,nco=nxo+1,lenco=nco*nm,lencob=lenco*lfac)
        parameter (nxl= 2,ncl=nxl+1,lencl=ncl*nm,lenclb=lencl*lfac)
        parameter (nxs=11,ncs=nxs+1,lencs=ncs*nm,lencsb=lencs*lfac)
        parameter (nxw=nxl+nxs)
        common/stdatm/pstd(nl),tstd(nl),wstd(nl),ostd(nl)
        common/taudwo/taud(nl),tauw(nl),tauo(nl)
        dimension temp(*),wvmr(*),ozmr(*),taut(*)
        dimension coefd(ncd,nm,nr),coefo(nco,nm,nr),coefl(ncl,nm,nr)
        dimension coefs(ncs,nm,nr),coefc(ncc,nm,nr),cbufs(lencs)
        dimension pavg(nm),tref(nm),wref(nm),oref(nm)
        dimension tavg(nm),wamt(nm),oamt(nm),secz(nm)
        dimension tauc(nl),tlas(nl),wlas(nl),olas(nl)
        dimension xdry(nxd,nm),xozo(nxo,nm),xwet(nxw,nm),xcon(nxc,nm)
        character*16 cfile
        character*72 ancil_data_path,pfaast_path
        character*3 comp(nk)
        data comp /'dry','ozo','wco','wtl','wts'/
        integer iuc(nk),lencf(nk)
        data lencf /lencdb,lencob,lenccb,lenclb,lencsb/
        logical flip(nk),newang,newatm
        data init/0/,tlas/nl*0./,wlas/nl*0./,olas/nl*0./,zlas/-999./
        integer*4 get_lun    !needed by Get_lun() call by AKH    
        
        secant(z)=1./cos(0.01745329*z)

        cfile = 'metsecgencom.dat'
        pfaast_path = trim(ancil_data_path)//"static/pfaast/"

        if(mets.ne.init) then
           !iux=70
           do l=1,nk
              cfile(10:12)=comp(l)
!              iux=iux+1
!              iuc(l)=iux
              iuc(l)=get_lun()    !akh added
              open(iuc(l),file=trim(pfaast_path)//cfile,
     +             recl=lencf(l),access='direct',
     +             status='old',err=200)
           enddo
c * determine if coefficients need to be byte-flipped
           koff=(mets-8)*(nr+1)+1
           krec=koff
           do k=1,nk
              lenk=lencf(k)/lfac
              read(iuc(k),rec=krec) (cbufs(j),j=1,lenk)
              flip(k)=cbufs(1).lt.1.0
           enddo
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
           if( flip(1) ) call flip_rtc(coefd,ncd,nm,nr)
           if( flip(2) ) call flip_rtc(coefo,nco,nm,nr)
           if( flip(3) ) call flip_rtc(coefc,ncc,nm,nr)
           if( flip(4) ) call flip_rtc(coefl,ncl,nm,nr)
           if( flip(5) ) call flip_rtc(coefs,ncs,nm,nr)
           call conpir(pstd,tstd,wstd,ostd,nl,1,pavg,tref,wref,oref)
           init=mets
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

        newang=zena.ne.zlas
        if(newang) then
           zsec=secant(zena)
           do l=1,nm
              secz(l)=zsec
           enddo
           zlas=zena
        endif

        if(newang.or.newatm) then
           call calpir(tref,wref,oref,tavg,wamt,oamt,pavg,secz,
     +                   nm,nxd,nxw,nxo,nxc,xdry,xwet,xozo,xcon)
        endif

        k=kban-3
c * dry
        call taudoc(ncd,nxd,nm,coefd(1,1,k),xdry,taud)
c * ozo
        call taudoc(nco,nxo,nm,coefo(1,1,k),xozo,tauo)
c * wet
        call tauwtr(ncs,ncl,nxs,nxl,nxw,nm,coefs(1,1,k),
     +                  coefl(1,1,k),xwet,tauw)
        call taudoc(ncc,nxc,nm,coefc(1,1,k),xcon,tauc)
        do j=1,nl
           tauw(j)=tauw(j)*tauc(j)
        enddo
c * total
        do j=1,nl
           taut(j)=taud(j)*tauo(j)*tauw(j)
        enddo
        return

200     return1
        end
