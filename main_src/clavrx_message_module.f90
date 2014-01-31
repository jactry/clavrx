! - $Header$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: clavrx_message_module.f90 (src)
!       clavrx_message_module (program)
!
! PURPOSE: 
!
! DESCRIPTION: 
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
!--------------------------------------------------------------------------------------

module clavrx_message_module

integer , parameter :: QUIET_LEV = 0
integer , parameter :: DEBUG_LEV =9

character(len = *) , parameter :: PROMPT = 'CLAVR-x>>'
integer :: VERBOSE_LEVEL = QUIET_LEV ! from quiet (0) to verbose(10)

public :: mesg

contains

subroutine mesg ( text, level , color , stop_evt )
  character (len = *) , intent(in) :: text
  integer, optional, intent(in) :: level
  integer, optional, intent(in) :: color
  logical, optional, intent(in) :: stop_evt
  character( len = 2) :: color_string
  
 
  write(color_string,'(I2)') color 
 
  if ( level .LE. VERBOSE_LEVEL) then
   print*,PROMPT//achar(27)//'['//color_string//'m '//text//achar(27)//'[0m'
  end if
end subroutine mesg


end module clavrx_message_module
