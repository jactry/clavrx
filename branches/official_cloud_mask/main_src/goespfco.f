c--------------------------------------------------------------------------
c Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
c
c NAME: goespfco.f (src)
c
c PURPOSE: 
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
        subroutine goespfco(ancil_data_path,ngoes,*)

! RCS Version: $Id:

c * Input GOES/I-M Planck-function, band-cor'n & TSKIN coeff's
c        plus gammas, deltas, epsilons, and 'use' flags
c        NGOES = GOES satellite number, e.g. 8 (GOES/I)
c .... version of 15.02.06

c   This routine MUST be called before using GOESBRIT, GOESPLAN, GOESDBDT, or GOESSKIN.

        parameter (iuc=15,lenc=400,nk=26,nt=2,lenp=nk*(nt+3),lent=30)
        parameter (leng=nk*3,lenu=nk*nt,lenr=lenc*4)
        parameter (mbg=200,mbu=300)
c       common/gimgdx/gamma(nk),delta(nk),epsln(nk)
        common/gimgdx/gbuf(leng)
c       common/plncgx/cwn(nk),fk1(nk),fk2(nk),tc(nt,nk)
        common/plncgx/pbuf(lenp)
        common/tskcfx/tbuf(lent)
        common/use/ibuf(lenu)
        dimension cbuf(lenc)
        character*12 cfile
        character*72 ancil_data_path,pfast_path 
        
        cfile = 'goescbnd.dat'
        pfast_path = trim(ancil_data_path)//"pfast/"

        open(iuc,file=trim(pfast_path)//cfile,recl=lenr,
     +          access='direct',status='old',
     +          err=100)
        irc=ngoes-7
        read(iuc,rec=irc) cbuf
        close(iuc)

        do l=1,lenp
           pbuf(l)=cbuf(l)
        enddo
        m=lenp
        do l=1,lent
           m=m+1
           tbuf(l)=cbuf(m)
        enddo
        m=mbg
        do l=1,leng
           m=m+1
           gbuf(l)=cbuf(m)
        enddo
        m=mbu
        do l=1,lenu
           m=m+1
           ibuf(l)=cbuf(m)
        enddo
        return
100     return1
        end
