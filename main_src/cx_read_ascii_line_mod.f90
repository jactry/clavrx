! $Id$
module cx_read_ascii_line_mod

   interface read_ascii_line
      module procedure  read_option_line_logical 
      module procedure read_option_line_int 
      module procedure  read_option_line_string 
      module procedure read_option_line_l6 
       module procedure read_option_line_real_9
   end interface read_ascii_line
   
contains   
   
   subroutine read_option_line_int ( lun , data_i )
      integer :: lun
      integer :: data_i
      
      read(unit=lun,fmt=*)  data_i
      
   end subroutine read_option_line_int
   
   subroutine read_option_line_logical( lun , data_l) 
      integer :: lun
      logical :: data_l
      
      integer :: int_dummy
      
      read(unit=lun,fmt=*)  int_dummy
      data_l = int_dummy == 1
   end subroutine read_option_line_logical 
   
    subroutine read_option_line_string ( lun , data_i )
      integer :: lun
      character ( len = * ) :: data_i
      
      read(unit=lun,fmt=*)  data_i
      
   end subroutine read_option_line_string
   
   subroutine read_option_line_l6 ( lun , data_i )
      implicit none
      integer :: lun
      logical :: data_i (6)
      integer :: int_dummy(6)
      
  
      read(unit=lun,fmt=*)  int_dummy
      
      
      data_i = int_dummy == 1
      
   end subroutine read_option_line_l6
   
      subroutine read_option_line_real_9 ( lun , data_i )
      implicit none
      integer :: lun
      real :: data_i (9)
      
      
  
      read(unit=lun,fmt=*)  data_i
      
      
     
      
   end subroutine read_option_line_real_9

end module cx_read_ascii_line_mod
