c--------------------------------------------------------------------------
c Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
c
c NAME: big_endian.f (src)
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
	logical function big_endian()
c .... version of 29.03.05
        implicit none
        integer*4 long
        integer*2 short(2)
        equivalence (short(1),long)

        long  = 0
        short(1) = 0
        short(2) = 1

        if(long == 1) then
           big_endian = .true.
        else
           big_endian = .false.
        endif

        end
