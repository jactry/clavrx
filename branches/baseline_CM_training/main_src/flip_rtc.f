c--------------------------------------------------------------------------
c Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
c
c NAME: flip_rtc.f (src)
c
c PURPOSE: Byte-flip radiative-transfer coefficients
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
	subroutine flip_rtc(cbuf,nco,nle,nch)
c .... version of 25.08.03

c * nco = number of coefficients
c * nle = number of levels
c * nch = number of channels

	dimension cbuf(nco,nle,nch)

	do k=1,nch
	   do j=1,nle
	      call re4flip(cbuf(1,j,k),cbuf(1,j,k),nco)
	   enddo
	enddo

	return
	end
