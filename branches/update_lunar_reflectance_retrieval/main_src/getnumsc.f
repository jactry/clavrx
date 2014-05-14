c--------------------------------------------------------------------------
c Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
c
c NAME: getnumsc.f (src)
c
c PURPOSE: Get "Satellite number" for use in r/t code
c
c DESCRIPTION: 
c  * For all spacecraft from TIROSN through METOPC!
c 
c  ... Input name may contain any combination of upper and lower case.
c 
c  ... NOAA s/c prior to NOAA10 may be input as either noaa0n or noaa-n
c 
c  * If given s/c not found, routine returns number = -1
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
	subroutine getnumsc(craft,number)
! .... version of 27.07.06

	implicit none
	integer isat,ksat,nsat,number
	integer k,kc,len,nc
	character*6 craft,craftl
	character*3 csat

	parameter (nsat=18)
	character*6 crafts(nsat)
	character cc,char,cm,cz

	data cm/'-'/,cz/'0'/

	data crafts/'tirosn','noaa06','noaa07','noaa08','noaa09',
     +              'noaa10','noaa11','noaa12','noaa13','noaa14',
     +              'noaa15','noaa16','noaa17','noaa18',
     +              'metopa','noaa19','metopb','metopc'/

! Convert spacecraft name to lower case if necessary

	nc=len(craft)
	do k=1,nc
	   cc=craft(k:k)
	   kc=ichar(cc)
	   if(kc.gt.64.and.kc.lt.91) cc=char(kc+32)
	   craftl(k:k)=cc
	enddo
	if(craftl(5:5) == cm) craftl(5:5)=cz

! Search list for given spacecraft

	ksat=0
	do isat=1,nsat
	   if(craftl == crafts(isat)) then
	      ksat=isat
	   endif
	enddo

	number=ksat-1

	return
	end
