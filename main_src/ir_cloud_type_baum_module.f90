!$Id$
module IR_CLOUD_TYPE_BAUM_MODULE
!====================================================================
! Module Name: IR_CLOUD_TYPE_BAUM_MODULE
!
! Function:  House routines needed for Dr Bryan Baum's
!            IR Cloud Phase
!
! Description:
!
! Dependencies:
!        none
!
! Restrictions:  None
!
! Reference:
!
! Author: Bryan "Buzz" Baum, UW/SSEC
!         Andrew Heidinger, NOAA/NESDIS
!
!
!
!Cloud Mask Values
! CLEAR = 0
! PROB_CLEAR = 1
! PROB_CLOUDY = 2
! CLOUDY = 3
!
!Cloud Type Values
! CLEAR_TYPE = 0
! PROB_CLEAR_TYPE = 1
! FOG_TYPE = 2
! WATER_TYPE = 3
! SUPERCOOLED_TYPE = 4
! MIXED_TYPE = 5
! OPAQUE_ICE_TYPE = 6
! CIRRUS_TYPE = 7
! OVERLAP_TYPE = 8
! OVERSHOOTING_TYPE = 9
! UNKNOWN_TYPE = 10
! DUST_TYPE = 11
! SMOKE_TYPE = 12
! FIRE_TYPE = 13
!
!Cloud Phase Values
! CLEAR_PHASE = 0
! WATER_PHASE = 1
! SUPERCOOLED_PHASE = 2
! MIXED_PHASE = 3
! ICE_PHASE = 4
! UNKNOWN_PHASE = 5
!====================================================================
  use CONSTANTS
  use PIXEL_COMMON, only: Image,  &
                          Ch,     &
                          Geo,    &
                          Sensor, &
                          Sfc,    &
                          Cld_Mask, &
                          Bad_Pixel_Mask, &
                          I_Nwp, &
                          J_Nwp, &
                          Zen_Idx_Rtm, &
                          Beta_11um_12um_Tropo_Rtm, &
                          Beta_11um_85um_Tropo_Rtm, &
                          Beta_11um_67um_Tropo_Rtm, &
                          Beta_11um_133um_Tropo_Rtm, &
                          Beta_11um_133fusum_Tropo_Rtm, &
                          Cld_Test_Vector_Packed, &
                          Cld_Phase_IR, &
                          Cld_Type_IR, &
                          Diag_Pix_Array_1, &
                          Diag_Pix_Array_2, &
                          Diag_Pix_Array_3

  implicit none

  public::  IR_CLOUD_TYPE_BAUM, SET_IR_CLOUD_TYPE_VERSION

  CONTAINS
!====================================================================
!  record cvs version as a global variable for output to hdf
!====================================================================
subroutine SET_IR_CLOUD_TYPE_VERSION()
   Cloud_Type_IR_Version = "$Id$"
end subroutine SET_IR_CLOUD_TYPE_VERSION

!====================================================================
! Subroutine Name: INFRARED_CLOUD_PHASE
!
! Function: Compute cloud phase from IR measurements
!
! Description:
! This is a modification of the MODIS C6 IR phase module:
!     a. Uses 13.3-micron channel rather than MODIS 7.3-micron channel
!        (7.3-um was used because it was a much cleaner channel with almost
!         no artifacts from striping or radiometric noise)
!     b. Assumes availability of 13.3-um channel for AVHRR-HIRS and VIIRS-CrIS
!        based on data fusion approach
!
! Dependencies: RTM data structure needed
!
! Restrictions:  11 um channel needs to be valid
!
! Reference:
!
! Author: Bryan Baum, credentials and character questionable at best
!
!====================================================================
subroutine IR_CLOUD_TYPE_BAUM()

  implicit none

  integer:: Line_Idx
  integer:: Elem_Idx
  integer:: Num_Elem
  integer:: Num_Line
  integer:: Nwp_Lon_Idx
  integer:: Nwp_Lat_Idx
  integer:: Vza_Idx 
  logical:: Fire_Flag

  real(kind=real4) :: BTD8511
  real(kind=real4) :: Beta_11um_133um

  !--- begin executable code
  Num_Line = Image%Number_Of_Lines_Per_Segment
  Num_Elem = Image%Number_Of_Elements

  Cld_Phase_IR = sym%UNKNOWN_PHASE
  Cld_Phase_IR = sym%UNKNOWN_TYPE

  Fire_Flag = .false.


  !------------------------------------------------------------------
  ! Step #0: Check the needed data are present
  !------------------------------------------------------------------
  if (Sensor%Chan_On_Flag_Default(29) == sym%NO) return
  if (Sensor%Chan_On_Flag_Default(31) == sym%NO) return
  if (Sensor%Chan_On_Flag_Default(32) == sym%NO) return
  if (Sensor%Chan_On_Flag_Default(33) == sym%NO .and. &
      Sensor%Chan_On_Flag_Default(45) == sym%NO) return

  !------------------------------------------------------------------
  ! Step #1: Check for non-cloud conditions and
  !          compute parameters needed in future steps
  !------------------------------------------------------------------
  Line_Loop: do Line_Idx = 1, Num_Line
  Element_Loop: do Elem_Idx = 1, Num_Elem

     !--- skip bad pixels
     if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) cycle

     !--- save indices
     Nwp_Lon_Idx = I_Nwp(Elem_Idx,Line_Idx)
     Nwp_Lat_Idx = J_Nwp(Elem_Idx,Line_Idx)
     Vza_Idx = Zen_Idx_Rtm(Elem_Idx,Line_Idx)

     !-------------------------------------------------------------
     ! Determine if a non-cloud type has been determined in the
     ! cloud mask, if so, set the type flag and exit
     !-------------------------------------------------------------
     Fire_Flag = BTEST(Cld_Test_Vector_Packed(2,Elem_Idx,Line_Idx), 7)
     if (Fire_Flag) then
        Cld_Type_IR(Elem_Idx,Line_Idx) = sym%FIRE_TYPE
        Cld_Phase_IR(Elem_Idx,Line_Idx) = sym%UNKNOWN_PHASE
        cycle
     endif

     !--- set clear to clear phase
     if (Cld_Mask(Elem_Idx,Line_Idx) == sym%CLEAR) then
          Cld_Phase_IR(Elem_Idx,Line_Idx) = sym%CLEAR_PHASE
          Cld_Type_IR(Elem_Idx,Line_Idx) = sym%CLEAR_TYPE
          cycle
     endif

     !--- set probably clear to clear phase
     if (Cld_Mask(Elem_Idx,Line_Idx) == sym%PROB_CLEAR) then
          Cld_Phase_IR(Elem_Idx,Line_Idx) = sym%CLEAR_PHASE
          Cld_Type_IR(Elem_Idx,Line_Idx) = sym%PROB_CLEAR_TYPE
          cycle
     endif


     ! ------------------------------------------
     ! the following block of code added by BB
     ! IR phase computation for all pixels not bad/clear/probably clear/fire
     !

     !---- compute brightness temp difference for convenience
     BTD8511 = ch(29)%Bt_Toa(Elem_Idx,Line_Idx) - ch(31)%Bt_Toa(Elem_Idx,Line_Idx)

     !--- use 13.3 fusion beta when availble
     if (Sensor%Chan_On_Flag_Default(33) == sym%YES) then
          Beta_11um_133um = Beta_11um_133um_Tropo_Rtm(Elem_Idx,Line_Idx)
     endif
     if (Sensor%Chan_On_Flag_Default(45) == sym%YES) then
          Beta_11um_133um = Beta_11um_133fusum_Tropo_Rtm(Elem_Idx,Line_Idx)
     endif

    ! one path for pixels over land; the other over water
    !    tests use Sfc%Land for this determination where
    !    Sfc%Land = 1  means land surface; all other values are water
    !
    if (ch(31)%Bt_Toa(Elem_Idx,Line_Idx) <= 238.) then
    ! NOTE: test not currently dependent on land/water surface type

          Cld_Phase_IR(Elem_Idx,Line_Idx) = sym%ICE_PHASE

    else if (Sfc%Land(Elem_Idx,Line_Idx) /= 1 .and. &
            BTD8511 >= 0.5) then

          cld_Phase_IR(Elem_Idx,Line_Idx) = sym%ICE_PHASE

    else if (Sfc%Land(Elem_Idx,Line_Idx) /= 1 .and. &
            BTD8511  >= -0.75    .and. &
            Beta_11um_85um_Tropo_Rtm(Elem_Idx,Line_Idx) <=  1.1     .and. &
            Beta_11um_12um_Tropo_Rtm(Elem_Idx,Line_Idx) >= 0.95) then

         cld_Phase_IR(Elem_Idx,Line_Idx) = sym%ICE_PHASE

    else if (Sfc%Land_Mask(Elem_Idx,Line_Idx) == 1 .and. &
            BTD8511              >= 1.0    .and. &
            Beta_11um_85um_Tropo_Rtm(Elem_Idx,Line_Idx) <= 1.1    .and. &
            Beta_11um_12um_Tropo_Rtm(Elem_Idx,Line_Idx) >= 0.95   .and. &
            Beta_11um_133um >= 0.5) then

        cld_Phase_IR(Elem_Idx,Line_Idx) = sym%ICE_PHASE

    else if (Sfc%Land_Mask(Elem_Idx,Line_Idx) == 1 .and. &
            BTD8511          >= 0.5    .and. &
            Beta_11um_85um_Tropo_Rtm(Elem_Idx,Line_Idx) < 1.1     .and. &
            Beta_11um_85um_Tropo_Rtm(Elem_Idx,Line_Idx) > 0.75    .and. &
            Beta_11um_12um_Tropo_Rtm(Elem_Idx,Line_Idx) >= 0.95   .and. &
            Beta_11um_133um >= 0.5) then

        cld_Phase_IR(Elem_Idx,Line_Idx) = sym%ICE_PHASE

    ! now done with ice phase tests; moving on to liquid water tests
    ! NOTE: the following tests not currently dependent on land/water surface type
    else if (ch(31)%Bt_Toa(Elem_Idx,Line_Idx) > 238.0               .and. &
            BTD8511 <= -1.0) then

        cld_Phase_IR(Elem_Idx,Line_Idx) = sym%WATER_PHASE

    else if (ch(31)%Bt_Toa(Elem_Idx,Line_Idx) > 285.0               .and. &
            BTD8511 <= -0.5) then

        cld_Phase_IR(Elem_Idx,Line_Idx) = sym%WATER_PHASE

    else if (Beta_11um_133um < 0.5   .and. &
             Beta_11um_133um > 0.0)  then

       cld_Phase_IR(Elem_Idx,Line_Idx) = sym%WATER_PHASE

    else
    ! everything left will be labeled as unknown phase even though the array
    ! was initialized to unknown phase; this step is likely unnecessary

       cld_Phase_IR(Elem_Idx,Line_Idx) = sym%UNKNOWN_PHASE

    endif
! end of block of code added by BB
! ------------------------------------------

  enddo Element_Loop
  enddo Line_Loop

  !--------------------------------------------------------------------
  !--- construct cloud types from cloud phases
  !--------------------------------------------------------------------
  Line_Loop_2: do Line_Idx = 1, Num_Line
  Element_Loop_2: do Elem_Idx = 1, Num_Elem
 
  if (Cld_Type_IR(Elem_Idx,Line_Idx) == sym%FIRE_TYPE) cycle
  if (Cld_Type_IR(Elem_Idx,Line_Idx) == sym%CLEAR_TYPE) cycle
  if (Cld_Type_IR(Elem_Idx,Line_Idx) == sym%PROB_CLEAR_TYPE) cycle

  if (Cld_Phase_IR(Elem_Idx,Line_Idx) == sym%CLEAR_PHASE) then

           if (Cld_Mask(Elem_Idx,Line_Idx)== sym%CLEAR) then 
              Cld_Type_IR(Elem_Idx,Line_Idx) = sym%CLEAR_TYPE
           else
              Cld_Type_IR(Elem_Idx,Line_Idx) = sym%PROB_CLEAR_TYPE
           endif

   elseif (Cld_Phase_IR(Elem_Idx,Line_Idx) == sym%WATER_PHASE) then

           Cld_Type_IR(Elem_Idx,Line_Idx) = sym%WATER_TYPE

   elseif (Cld_Phase_IR(Elem_Idx,Line_Idx) == sym%SUPERCOOLED_PHASE) then

           Cld_Type_IR(Elem_Idx,Line_Idx) = sym%SUPERCOOLED_TYPE

   elseif (Cld_Phase_IR(Elem_Idx,Line_Idx) == sym%MIXED_PHASE) then

           Cld_Type_IR(Elem_Idx,Line_Idx) = sym%SUPERCOOLED_TYPE

   elseif (Cld_Phase_IR(Elem_Idx,Line_Idx) == sym%ICE_PHASE) then

           if (ch(31)%Emiss_Tropo(Elem_Idx,Line_Idx) > 0.95) then 
              Cld_Type_IR(Elem_Idx,Line_Idx) = sym%OVERSHOOTING_TYPE
           elseif (ch(31)%Emiss_Tropo(Elem_Idx,Line_Idx) > 0.80) then 
              Cld_Type_IR(Elem_Idx,Line_Idx) = sym%OPAQUE_ICE_TYPE
           else
              Cld_Type_IR(Elem_Idx,Line_Idx) = sym%CIRRUS_TYPE
           endif

   elseif (Cld_Phase_IR(Elem_Idx,Line_Idx) == sym%UNKNOWN_PHASE) then

            Cld_Type_IR(Elem_Idx,Line_Idx) = sym%UNKNOWN_TYPE

   else

            Cld_Type_IR(Elem_Idx,Line_Idx) = sym%UNKNOWN_TYPE

  endif

  enddo Element_Loop_2
  enddo Line_Loop_2


end subroutine IR_CLOUD_TYPE_BAUM

!----------------------------------------------------------------------
! End of Module
!----------------------------------------------------------------------

end module IR_CLOUD_TYPE_BAUM_MODULE
