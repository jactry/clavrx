! $Id: awg_cloud_height.f90 583 2014-10-08 03:43:36Z heidinger $
module ACHA_COMP
!---------------------------------------------------------------------
! This module houses the routines associated with...
!
! ACHA Cloud Optical Microphysical Properties
!
! Author: 
!
! Reference:
!
!----------------------------------------------------------------------
  use ACHA_SERVICES_MOD !acha_services_mod.f90 in akh_clavrx_src

  implicit none

  public:: ACHA_COMP_ALGORITHM

  private:: COMPUTE_TAU_REFF_ACHA 

  !--- include the non-system specific variables
  include 'acha_parameters.inc'

  real, private, parameter:: MISSING_VALUE_REAL = -999.0
  integer, private, parameter:: MISSING_VALUE_INTEGER = -999

  contains 

!------------------------------------------------------------------------------
! AWG Cloud Optical and Microphysical Algorithm (ACHA-COMP)
!
! Author: Andrew Heidinger, NOAA
!
! Assumptions
!
! Limitations
!
! NOTE.  This algorithm use the same input and output structures as 
!        the AWG_CLOUD_HEIGHT_ALGORITHM.
!        Do not overwrite elements of the Output structure expect those
!        generated here.
!
!      Output%Tau
!      Output%Ec  (modified from ACHA)
!      Output%Reff
!      Output%Beta (modified from ACHA)
!
!----------------------------------------------------------------------
! modification history
!
!------------------------------------------------------------------------------
  subroutine ACHA_COMP_ALGORITHM(Input, Symbol, Output)

  !===============================================================================
  !  Argument Declaration
  !==============================================================================

  type(symbol_acha), intent(inout) :: Symbol
  type(acha_input_struct), intent(inout) :: Input
  type(acha_output_struct), intent(inout) :: Output

  !===============================================================================
  !  Local Variable Declaration
  !===============================================================================

  integer:: Elem_Idx
  integer:: Line_Idx

  integer:: Cloud_Type
  integer:: Cloud_Phase

  !--- scalar local variables
  real (kind=real4):: Ec_Slant

!-----------------------------------------------------------------------
! BEGIN EXECUTABLE CODE
!-----------------------------------------------------------------------

   !-------------------------------------------------------------------------
   ! Initialization
   !-------------------------------------------------------------------------

   !--------------------------------------------------------------------------
   ! loop over pixels in scanlines
   !--------------------------------------------------------------------------
   Line_loop: do Line_Idx = Line_Idx_min,Input%Number_of_Lines + Line_Idx_min - 1

    Element_Loop:   do Elem_Idx = 1, Input%Number_of_Elements

    !--- for convenience, save nwp indices to local variables
    Cloud_Type = Input%Cloud_Type(Elem_Idx,Line_Idx)
    
    !----------------------------------------------------------------------
    ! determine cloud phase from cloud type for convenience
    !----------------------------------------------------------------------
    Cloud_Phase = symbol%UNKNOWN_PHASE

    if ( (Cloud_Type  == symbol%FOG_TYPE) .or. &
       (Cloud_Type  == symbol%WATER_TYPE) .or. &
       (Cloud_Type  == symbol%SUPERCOOLED_TYPE)) then
       Cloud_Phase = symbol%WATER_PHASE
    endif

    if ( (Cloud_Type  == symbol%CIRRUS_TYPE) .or. &
        (Cloud_Type  == symbol%OVERLAP_TYPE) .or. &
        (Cloud_Type  == symbol%TICE_TYPE) .or. &
        (Cloud_Type  == symbol%OVERSHOOTING_TYPE)) then
        Cloud_Phase = symbol%ICE_PHASE
    endif

    !-----------------------------------------------------------------------------
    ! Estimate Cloud Optical and Microphysical Properties
    !-----------------------------------------------------------------------------
    Ec_Slant =  Output%Ec(Elem_Idx,Line_Idx)

    if (Output%Zc(Elem_Idx,Line_Idx) /= MISSING_VALUE_REAL .and. &
        Output%Ec(Elem_Idx,Line_Idx) /= MISSING_VALUE_REAL .and. &
        Output%Beta(Elem_Idx,Line_Idx) /= MISSING_VALUE_REAL) then


     !--- save nadir adjusted emissivity and optical depth
     if (Output%Ec(Elem_Idx,Line_Idx) < 1.00) then

       call COMPUTE_TAU_REFF_ACHA(symbol, &
                              Output%Beta(Elem_Idx,Line_Idx), &
                              Input%Cosine_Zenith_Angle(Elem_Idx,Line_Idx), &
                              Cloud_Phase, &
                              Ec_Slant, &
                              Output%Ec(Elem_Idx,Line_Idx), &
                              Output%Tau(Elem_Idx,Line_Idx), &
                              Output%Reff(Elem_Idx,Line_Idx))

     else

       Output%Tau(Elem_Idx,Line_Idx) = 20.0
       Output%Ec(Elem_Idx,Line_Idx) = 1.0
       Output%Reff(Elem_Idx,Line_Idx) = 20.0

       if( Cloud_Phase == symbol%ICE_PHASE) then
          Output%Beta(Elem_Idx,Line_Idx) = Beta_Ap_Water
       else
          Output%Beta(Elem_Idx,Line_Idx) = Beta_Ap_Ice
       endif

     endif

    endif
   
 end do Element_Loop

end do Line_Loop

end subroutine  ACHA_COMP_ALGORITHM

!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
subroutine COMPUTE_TAU_REFF_ACHA(symbol,  &
                                 Beta, &
                                 Cosine_Zenith_Angle, &
                                 Cloud_Phase, &
                                 Ec_Slant, & 
                                 Ec, &
                                 Tau, &
                                 Reff)

   type(symbol_acha), intent(in) :: symbol
   real(kind=real4), intent(in):: Beta
   real(kind=real4), intent(in):: Cosine_Zenith_Angle
   real(kind=real4), intent(in):: Ec_Slant
   integer(kind=int4), intent(in):: Cloud_Phase
   real(kind=real4), intent(out):: Ec
   real(kind=real4), intent(out):: Tau
   real(kind=real4), intent(out):: Reff
   real(kind=real4):: Qe_vis
   real(kind=real4):: Qe_11um
   real(kind=real4):: wo_11um
   real(kind=real4):: g_11um
   real(kind=real4):: Tau_Abs_11um
   real(kind=real4):: Temp_R4
   real(kind=real4):: log10_Reff

   Tau = MISSING_VALUE_REAL
   Reff = MISSING_VALUE_REAL
   Ec = MISSING_VALUE_REAL

   if (Cloud_Phase == symbol%ICE_PHASE) then
    Temp_R4 = A_Re_Beta_FIT_ICE +  &
              B_Re_Beta_FIT_ICE*Beta + &
              C_Re_Beta_FIT_ICE*Beta**2 + &
              D_Re_Beta_FIT_ICE*Beta**3
   else
    Temp_R4 = A_Re_Beta_FIT_WATER +  &
              B_Re_Beta_FIT_WATER*Beta + &
              C_Re_Beta_FIT_WATER*Beta**2 + &
              D_Re_Beta_FIT_WATER*Beta**3
   endif

   Reff = 10.0**(1.0/Temp_R4)

   if (Reff > 0.0) then
     log10_Reff = alog10(Reff)
   else
     return
   endif

   if (Cloud_Phase == symbol%ICE_PHASE) then

     Qe_Vis = A_Qe_065um_FIT_ICE +  &
              B_Qe_065um_FIT_ICE*log10_Reff + &
              C_Qe_065um_FIT_ICE*log10_Reff**2

     Qe_11um = A_Qe_11um_FIT_ICE +  &
               B_Qe_11um_FIT_ICE*log10_Reff + &
               C_Qe_11um_FIT_ICE*log10_Reff**2

     wo_11um = A_wo_11um_FIT_ICE +  &
               B_wo_11um_FIT_ICE*log10_Reff + &
               C_wo_11um_FIT_ICE*log10_Reff**2

     g_11um =  A_g_11um_FIT_ICE +  &
               B_g_11um_FIT_ICE*log10_Reff + &
               C_g_11um_FIT_ICE*log10_Reff**2

   else

     Qe_Vis = A_Qe_065um_FIT_WATER +  &
              B_Qe_065um_FIT_WATER*log10_Reff + &
              C_Qe_065um_FIT_WATER*log10_Reff**2

     Qe_11um = A_Qe_11um_FIT_WATER +  &
               B_Qe_11um_FIT_WATER*log10_Reff + &
               C_Qe_11um_FIT_WATER*log10_Reff**2

     wo_11um = A_wo_11um_FIT_WATER +  &
               B_wo_11um_FIT_WATER*log10_Reff + &
               C_wo_11um_FIT_WATER*log10_Reff**2

     g_11um =  A_g_11um_FIT_WATER +  &
               B_g_11um_FIT_WATER*log10_Reff + &
               C_g_11um_FIT_WATER*log10_Reff**2

   endif

   Tau_Abs_11um = -Cosine_Zenith_Angle*alog(1.0 - Ec_Slant) 

   Ec = 1.0 - exp(-Tau_Abs_11um)
   
   Tau = (Qe_vis / Qe_11um) * Tau_Abs_11um / (1.0 - wo_11um * g_11um)

   !--- set negative values to be missing - added by Y Li
   if (Tau < 0) then
      Tau = MISSING_VALUE_REAL
      Reff= MISSING_VALUE_REAL
   endif

end subroutine COMPUTE_TAU_REFF_ACHA 

!----------------------------------------------------------------------
! End of Module
!----------------------------------------------------------------------

end module ACHA_COMP
