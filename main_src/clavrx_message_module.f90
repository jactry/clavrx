! - $Header: https://svn.ssec.wisc.edu/repos/cloud_team_clavrx/baseline_CM_training/main_src/clavrx_message_module.f90 1403 2015-10-15 21:40:38Z dbotambekov $
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

type verbose_type
   integer :: QUIET   =  0
   integer :: ERROR =  1
   integer :: MINIMAL =  2
   
   integer :: WARNING =  4
   integer :: DEFAULT =  5
   integer :: VERBOSE =  9 
end type verbose_type

type (verbose_type ) , save :: verb_lev 


character(len = *) , parameter :: PROMPT = 'CLAVR-x >'
integer :: VERBOSE_LEVEL = 5  ! from quiet (0) to verbose(10)

public :: mesg

interface mesg
   module procedure mesg_pure
   module procedure mesg_1r
   module procedure mesg_1d
end interface mesg


contains


   subroutine do_it ( text, message_level  )
      use file_tools, only: file_test
      character ( len = * ) , intent (in) :: text 
      integer , intent ( in ) :: message_level
      integer :: verbose_level
      
      verbose_level = verb_lev % DEFAULT
      if ( file_test ( 'verbose_level.txt') ) then
         open ( 37, file = 'verbose_level.txt' )
         read (37,'(i1)') verbose_level
         close (37)
      end if
      
   
      
      if ( message_level <= verbose_level ) then
         
         print*,PROMPT,' ', text 
        
      end if

   end subroutine do_it

   subroutine mesg_pure ( text, level , color )
      character (len = *) , intent(in) :: text
      integer, optional, intent(in) :: level
      integer, optional, intent(in) :: color
      character( len = 2) :: color_string
      integer :: lev
      
      lev = verb_lev % DEFAULT
      if ( present (level)) lev = level
      color_string=''  
      if (present(color)) write(color_string,'(I2)') color 
      
      call do_it ( trim(text), lev ) 
      
     
 
   end subroutine mesg_pure

   subroutine mesg_1r ( text,  param_r, level , color )
      character (len = *) , intent(in) :: text
      real, intent(in) :: param_r
      integer, optional, intent(in) :: level
      integer, optional, intent(in) :: color
      character( len = 2) :: color_string
      integer :: lev
      character ( len =100 ) :: string_100
      character ( len =200) :: text_1
 
  
      write ( string_100, '(f20.4)') param_r
  
      lev = verb_lev % DEFAULT
      if ( present (level)) lev = level
      text_1 = text//trim(string_100)
      color_string=''  
      if (present(color)) write(color_string,'(I2)') color 
      
      call do_it ( trim(text_1), lev ) 
      
    
 
   end subroutine mesg_1r

   subroutine mesg_1d ( text,  param_d, level , color )
      character (len = *) , intent(in) :: text
      real(8), intent(in) :: param_d
      integer, optional, intent(in) :: level
      integer, optional, intent(in) :: color
      character( len = 2) :: color_string
      integer :: lev
      character ( len =100 ) :: string_100
      character ( len =200) :: text_1


      write ( string_100, '(f25.4)') param_d

      lev = verb_lev % DEFAULT
      if ( present (level)) lev = level
      text_1 = text//trim(string_100)
      color_string=''
      if (present(color)) write(color_string,'(I2)') color

      call do_it ( trim(text_1), lev )

   end subroutine mesg_1d


end module clavrx_message_module
