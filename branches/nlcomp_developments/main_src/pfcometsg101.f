c--------------------------------------------------------------------------
c Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
c
c NAME: pfcometsg101.f (src)
c
c PURPOSE: 
c * Input METEOSAT SECOND GENERATION Planck-function and band-correction 
c   coefficients.
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

        subroutine pfcometsg101(ancil_data_path,meteo,*)

! RCS Version: $Id:

c * Input METEOSAT SECOND GENERATION Planck-function and band-correction coefficients.
c .... version of 29.06.05

c.... ancilary path

c * meteo = Meteosat number (8...10)

        parameter (iuc=15,lenc=100,nk=8,nt=2,lenp=nk*(nt+3))
        parameter (lencb=lenc*4)
c       parameter (kb=4,ke=11)
c       common/plnmsg/cwn(kb:ke),fk1(kb:ke),fk2(kb:ke),tc(nt,kb:ke)
        common/plnmsg/pbuf(lenp)
        dimension cbuf(lenc)
        character*16 cfile
        character*72 ancil_data_path,pfast_path 
        
        cfile = 'metsecgenbnd.dat'
        pfast_path = trim(ancil_data_path)//"pfast/"

c       print *, pfast_path
        open(iuc,file=trim(pfast_path)//cfile,
     +          recl=lencb,access='direct',status='old',
     +          err=100)
        irec=meteo-7
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
100     return1
        end
