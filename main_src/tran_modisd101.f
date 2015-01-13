c--------------------------------------------------------------------------
c Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
c
c NAME: tran_modisd101.f (src)
c
c PURPOSE: 
c * MODIS band/detector 101-level fast transmittance routine
c
c DESCRIPTION: 
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
      subroutine tran_modisd101(ancil_data_path, 
     +                          temp,wvmr,ozmr,theta,
     +                          craft,kban,jdet,
     +                          taut,iok)
c * MODIS band/detector 101-level fast transmittance routine
c .... version of 03.04.07

c      temp = temperature (Kelvin) profile
c      wvmr = water-vapor mixing-ratio (g/kg) profile
c      ozmr = ozone mixing-ratio (ppmv) profile
c     theta = local zenith angle
c     craft = TERRA, AQUA (either upper or lower case)
c      kban = band number     (20...36)
c      jdet = detector number (0...10) [ Product Order ]
c            detector 0 is based on band-average response functions

c            taut = total transmittance (see note below)
c             iok = 0 if successful, 1 if I/O problem

c * NOTE: for kban = 26, return-arrays are filled with 1.0

c * PLOD/PFAAST regression model based on LBLRTM line-by-line transmittances.
c * Input temperatures, and water-vapor and ozone mixing ratios, must
c *   be defined at the pressure levels in array 'pstd'
c *    (see block data 'reference_atmosphere').
c * Units: temperature, deg-K; water vapor, g/kg; ozone, ppmv.
c * Logical units 31-35 are used for coefficient files.
c * Logical units 81-85 are used for coefficient files. (AKH change)
c * Component tau's are returned through common, product in 'taut'.
c
c        incorporated get_lun routine from file_utility.f90 to remove
c        conflicts with fixed lun values.  Note, I had to extract this
c        function from the module because I could not successfully
c        access the function within the module
c        2010/11/23      A.Heidinger
 
      parameter (nd=10,nk=5,nl=101,nm=nl-1,koff=19,nr=17)
      parameter (nxc= 4,ncc=nxc+1,lencc=ncc*nm,lenccb=lencc*4)
      parameter (nxd= 8,ncd=nxd+1,lencd=ncd*nm,lencdb=lencd*4)
      parameter (nxo= 9,nco=nxo+1,lenco=nco*nm,lencob=lenco*4)
      parameter (nxl= 2,ncl=nxl+1,lencl=ncl*nm,lenclb=lencl*4)
      parameter (nxs=11,ncs=nxs+1,lencs=ncs*nm,lencsb=lencs*4)
      parameter (ndt=nd+1,nrps=nr*ndt,nxw=nxl+nxs)
      common/stdatm/pstd(nl),tstd(nl),wstd(nl),ostd(nl)
      common/taudwo/taud(nl),tauw(nl),tauo(nl)
      dimension temp(*),wvmr(*),ozmr(*),taut(*)
      dimension coefd(ncd,nm,0:nd,nr),coefo(nco,nm,0:nd,nr),
     +           coefl(ncl,nm,0:nd,nr),coefs(ncs,nm,0:nd,nr),
     +           coefc(ncc,nm,0:nd,nr)
      dimension bufc(lencc),bufd(lencd),bufo(lenco),
     +           bufl(lencl),bufs(lencs)
      dimension pavg(nm),tref(nm),wref(nm),oref(nm)
      dimension tavg(nm),wamt(nm),oamt(nm),secz(nm)
      dimension tauc(nl),tlas(nl),wlas(nl),olas(nl)
      dimension xdry(nxd,nm),xozo(nxo,nm),xwet(nxw,nm),xcon(nxc,nm)
      character*24 cfile(nk),xfile/'modisdet.com.101.xxx_end'/
      character*72 ancil_data_path,path
c     character*14 path/'./pfaast_coef/'/
      character*6 craft,cinit/'zzzzzz'/
      character*6 cbt/'TERRA'/,cba/'AQUA'/
      character*6 cst/'terra'/,csa/'aqua'/
      character*3 comp(nk)/'dry','ozo','wts','wtl','wco'/
      character*3 cbe/'big'/,cle/'lit'/
      integer*4 lengcf(nk)/lencdb,lencob,lencsb,lenclb,lenccb/
      integer*4 lengcx(nk)/lencd,lenco,lencs,lencl,lencc/
      integer*4 iuc(nk)
      integer*4 get_lun
      logical big_endian,newang,newatm
      data tlas/nl*0./,wlas/nl*0./,olas/nl*0./,zlas/-999./
      secant(z)=1./cos(0.01745329*z)

      save

      !--- initialize this (AKH)
      iok = 0

      path = trim(ancil_data_path) // "static/pfaast/"

      if(craft .ne. cinit) then
         if(craft == cbt .or. craft == cst) then
            ksat=1
         elseif(craft == cba .or. craft == csa) then
            ksat=2
         else
            write(0,'(''In tran_modisd101 -- unknown spacecraft '',a6)') craft
            stop
            go to 200
         endif

c * determine which coefficient files to use
         if(big_endian()) then
            xfile(18:20)=cbe
         else
            xfile(18:20)=cle
         endif

c * define and open the files
       ! iux=30
       ! iux=80
         do m=1,nk
       !    iux=iux+1
            iux = get_lun()
            xfile(10:12)=comp(m)
            lencf=lengcf(m)
c           open(iux,file=xfile,recl=lencf,access='direct',
c             print *, 'opening ', path//xfile
            open(iux,file=trim(path)//xfile,recl=lencf,access='direct',
     +                status='old',err=200)
            iuc(m)=iux
            cfile(m)=xfile
         enddo

c * first read each file's fill-record for band 26/det 0
c           and verify satellite number stored in word 1
c * note: number of levels is in word 2, creation date (yyyyddd) is in word 3
         ikrec=nrps*(ksat-1)
         krecx=ikrec+7
         do k=1,nk
            lencx=lengcx(k)
            read(iuc(k),rec=krecx) (bufs(j),j=1,lencx)
            nsat=bufs(1)
            if(nsat.ne.ksat) then
               xfile=cfile(k)
               go to 100
            endif
         enddo

c * now read in the coefficients
         krec=ikrec
         do l=0,nd
            do k=1,nr
              krec=krec+1
              read(iuc(1),rec=krec) ((coefd(i,j,l,k),i=1,ncd),j=1,nm)
              read(iuc(2),rec=krec) ((coefo(i,j,l,k),i=1,nco),j=1,nm)
              read(iuc(3),rec=krec) ((coefs(i,j,l,k),i=1,ncs),j=1,nm)
              read(iuc(4),rec=krec) ((coefl(i,j,l,k),i=1,ncl),j=1,nm)
              read(iuc(5),rec=krec) ((coefc(i,j,l,k),i=1,ncc),j=1,nm)
            enddo
         enddo
         do k=1,nk
            close(iuc(k))
         enddo

         call conpir(pstd,tstd,wstd,ostd,nl,1,pavg,tref,wref,oref)
         cinit=trim(craft)
         iok=0
      endif

c * if ozone profile is null, put in std-atm
      if(ozmr(1) == 0.) then
         do j=1,nl
            ozmr(j)=ostd(j)
         enddo
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
     +                   nm,nxd,nxw,nxo,nxc,xdry,xwet,xozo,xcon)
      endif

      if(kban == 26) then
         do j=1,nl
            taud(j)=1.0
            tauo(j)=1.0
            tauw(j)=1.0
            taut(j)=1.0
         enddo
         return
      endif

      j=jdet
      k=kban-koff
c * dry
      call taudoc(ncd,nxd,nm,coefd(1,1,j,k),xdry,taud)
c * ozo
      call taudoc(nco,nxo,nm,coefo(1,1,j,k),xozo,tauo)
c * wet
      call tauwtr(ncs,ncl,nxs,nxl,nxw,nm,coefs(1,1,j,k),
     +                 coefl(1,1,j,k),xwet,tauw)
      call taudoc(ncc,nxc,nm,coefc(1,1,j,k),xcon,tauc)
      do j=1,nl
         tauw(j)=tauw(j)*tauc(j)
      enddo
c * total
      do j=1,nl
         taut(j)=taud(j)*tauo(j)*tauw(j)
      enddo

c       print *, "end of trans"

      return
100   write(0,'(''In tran_modisd101 ... requested data for '',
     +           ''satellite '',i1/'' but read data for '',
     +           ''satellite '',i1,'' from file '',a24)')
     +       ksat,nsat,xfile
200   iok=1
      return
      end
