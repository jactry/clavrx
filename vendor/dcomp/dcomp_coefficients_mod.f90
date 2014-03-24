! $Id$
module dcomp_coefficients
   
   type sensor_constants_type 
      character ( len = 10) :: sat_name 
      real :: nu_20    
      real :: a1_20
      real :: a2_20
      real :: solar
      real :: ew
   end type sensor_constants_type
   
   
contains

   subroutine read_instr_constants(Instr_Const_file)
      implicit none
	  
      character(len=*), intent(in):: Instr_Const_file
      integer:: ios0, erstat
      integer:: Instr_Const_lun
      real:: dummy

      Instr_Const_lun = get_lun ( )
      open(unit=Instr_Const_lun,file=trim(Instr_Const_file) &
	    ,status="old",position="rewind",action="read",iostat=ios0)

      print *, "opening ", trim(Instr_Const_file)
      erstat = 0
      if (ios0 /= 0) then
         erstat = 19 
         print *, EXE_PROMPT &
		      , "Error opening VIIRS constants file, ios0 = ", ios0
         stop 19
      end if

      read(unit=Instr_Const_lun,fmt="(a3)") sat_name
      read(unit=Instr_Const_lun,fmt=*) Solar_Ch20
      read(unit=Instr_Const_lun,fmt=*) Ew_Ch20
      read(unit=Instr_Const_lun,fmt=*) a1_20, a2_20,nu_20
      read(unit=Instr_Const_lun,fmt=*) a1_22, a2_22,nu_22
      read(unit=Instr_Const_lun,fmt=*) a1_29, a2_29,nu_29
      read(unit=Instr_Const_lun,fmt=*) a1_31, a2_31,nu_31
      read(unit=Instr_Const_lun,fmt=*) a1_32, a2_32,nu_32
      read(unit=Instr_Const_lun,fmt=*) a1_40, a2_40,nu_40
      read(unit=Instr_Const_lun,fmt=*) a1_41, a2_41,nu_41
      read(unit=Instr_Const_lun,fmt=*) b1_day_mask,b2_day_mask,b3_day_mask,b4_day_mask
      close(unit=Instr_Const_lun)

      !-- convert solar flux in channel 20 to mean with units mW/m^2/cm^-1
      Solar_Ch20_Nu = 1000.0 * Solar_Ch20 / Ew_Ch20

   end subroutine READ_VIIRS_INSTR_CONSTANTS




end module dcomp_coefficients
