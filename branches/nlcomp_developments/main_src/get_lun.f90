!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: get_lun.f90 (src)
!       get_lun (program)
!
! PURPOSE: PUBLIC function to obtain a free logical unit number for file access
!
! DESCRIPTION: 
!       The search for a free logical unit number begins at 10. The logical
!       unit number if tested to see if it is connected to an open file. If
!       so, it is incremented by 1. This is repeated until a free logical
!       unit number is found.
!
! AUTHORS:
!  Andrew Heidinger, Andrew.Heidinger@noaa.gov
!  Andi Walther, CIMSS, andi.walther@ssec.wisc.edu
!  Denis Botambekov, CIMSS, denis.botambekov@ssec.wisc.edu
!  William Straka, CIMSS, wstraka@ssec.wisc.edu
!
! COPYRIGHT
! THIS SOFTWARE AND ITS DOCUMENTATION ARE CONSIDERED TO BE IN THE PUBLIC
! DOMAIN AND THUS ARE AVAILABLE FOR UNRESTRICTED PUBLIC USE. THEY ARE
! FURNISHED "AS IS." THE AUTHORS, THE UNITED STATES GOVERNMENT, ITS
! INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND AGENTS MAKE NO WARRANTY,
! EXPRESS OR IMPLIED, AS TO THE USEFULNESS OF THE SOFTWARE AND
! DOCUMENTATION FOR ANY PURPOSE. THEY ASSUME NO RESPONSIBILITY (1) FOR
! THE USE OF THE SOFTWARE AND DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL
! SUPPORT TO USERS.
!
! CALLING SEQUENCE:
!       result = get_lun()
!
! INPUT ARGUMENTS:
!       None.
!
! OUTPUT ARGUMENTS:
!       None.
!
! FUNCTION RESULT:
!       Function returns a default integer that can be used as a logical unit
!       number to open and access a file.
!
! CALLS:
!       None.
!
! EXTERNALS:
!       None
!
! COMMON BLOCKS:
!       None.
!
! SIDE EFFECTS:
!       None known.
!
! RESTRICTIONS:
!       None.
!--------------------------------------------------------------------------------------
  FUNCTION get_lun() RESULT( lun )


    ! -----------------
    ! Type declarations
    ! -----------------
 
    INTEGER :: lun
    LOGICAL :: file_open


    ! --------------------------------------------
    ! Initialise logical unit number and file_open
    ! --------------------------------------------

    lun = 9
    file_open = .TRUE.


    ! ------------------------------
    ! Start open loop for lun search
    ! ------------------------------

    lun_search: DO

      ! -- Increment logical unit number
      lun = lun + 1

      ! -- Check if file is open
      INQUIRE( lun, OPENED = file_open )

      ! -- Is this lun available?
      IF ( .NOT. file_open ) EXIT lun_search

    END DO lun_search

  END FUNCTION get_lun
