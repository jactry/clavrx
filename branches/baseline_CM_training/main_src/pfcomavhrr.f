c--------------------------------------------------------------------------
c Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
c
c NAME: pfcomavhrr.f (src)
c
c PURPOSE: 
c * Input AVHRR Planck-function and band-correction coefficients.
c ++++ for TIROS-N, NOAAA-6 ... NOAA-18, METOP-A, ff.
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
      subroutine pfcomavhrr(ancil_data_path,craft,*)
c .... version of 15.06.06

c     craft = spacecraft: tirosn,noaa06...noaa18, metopa, ... upper or lower case

      parameter (iuc=70,lenc=100,nb=3,nt=2,lenp=nb*(nt+3))
      parameter (lencb=lenc*4)
c     parameter (kb=3,ke=5)       !commented in original
c     common/avhpfc/wnum(kb:ke),fk1(kb:ke),fk2(kb:ke),tc(nt,kb:ke)  !commented in original
      common/avhpfc/pbuf(lenp)
      dimension cbuf(lenc)
      character*72 ancil_data_path,pfaast_path
      character*12 cfile
      character*6 craft

      cfile="avhrnbnd.dat"
      pfaast_path = trim(ancil_data_path)//"static/pfaast/"

      open(iuc,file=trim(pfaast_path)//cfile, 
     +          recl=lencb,access='direct',status='old',
     +            err=100)
      call getnumsc(craft,noff)
      irec=noff+1
      read(iuc,rec=irec) cbuf
      close(iuc)

c * determine if coefficients need to be byte-flipped
      if(cbuf(lenc).lt.1.0) then
         call re4flip(cbuf,cbuf,lenc)
      endif

      do i=1,lenp
         pbuf(i)=cbuf(i)
      enddo
      return

100   return1
      end
