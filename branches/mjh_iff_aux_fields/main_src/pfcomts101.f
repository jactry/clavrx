c--------------------------------------------------------------------------
c Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
c
c NAME: pfcomts101.f (src)
c
c PURPOSE: 
c * Input MTSAT Planck-function and band-correction coefficients
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
	subroutine pfcomts101(ancil_data_path,isat,*)
c * Input MTSAT Planck-function and band-correction coefficients
c .... version of 20.08.07
c added ancil and pfaast path - WCS3

        parameter (iuc=70,irec=1,lenc=100,nc=5,nt=2)
        parameter (lenp=nc*(nt+3),lenb=lenc*4)
        common/mtsatpfc/cwn(nc),fk1(nc),fk2(nc),tc(nt,nc)
        dimension cbuf(lenc),pbuf(lenp)
        equivalence (cwn(1),pbuf(1))
        character*16 cfile
        character*72 ancil_data_path,pfaast_path 
        
    
        cfile = 'mtsatbnd_101.dat'
        pfaast_path = trim(ancil_data_path)//"static/pfaast/"

        open(iuc,file=trim(pfaast_path)//cfile,
     +          recl=lenb,access='direct',status='old',
     +          err=100)
        krec=isat
        read(iuc,rec=krec) cbuf
        close(iuc)

        do i=1,lenp
           pbuf(i)=cbuf(i)
        enddo
        return

100     return1
	end
