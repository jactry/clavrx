module eumetset_overlap

! This module tests the ovelap alogirithm

!=============
   use ACHA_SERVICES_MOD, only : &
           real4, int1, int4, &
           acha_rtm_nwp_struct, ACHA_FETCH_PIXEL_NWP_RTM, &
           acha_output_struct,acha_symbol_struct, &
           acha_input_struct, acha_diag_struct
 
   use pixel_common, only : Covar_Ch27_Ch31_5x5

  implicit none
 
  public :: ovelap_co2_h2o_test

  private :: CO2IRW_CLOUD_HEIGHT
  private :: HEIGHT_H2O_CHANNEL
  private:: NULL_PIX_POINTERS

  real, private, parameter:: MISSING_VALUE_REAL = -999.0
  integer, private, parameter:: MISSING_VALUE_INTEGER = -999
  type(acha_symbol_struct), private :: symbol
  type(acha_rtm_nwp_struct),private :: ACHA_RTM_NWP

  integer, private:: Sfc_Level_RTM
  integer, private:: Tropo_Level_RTM
  integer(kind=INT4), parameter, private:: Line_Idx_Min = 1

  contains

!------------------------------------------------------------------------------

!------------------------------------------------------------------------------

  subroutine ovelap_co2_h2o_test(Input,Symbol, Output, Diag)

  !===============================================================================
  !  Argument Declaration
  !==============================================================================
  implicit none

  type(acha_input_struct), intent(inout) :: Input
  type(acha_symbol_struct), intent(inout) :: Symbol
  type(acha_output_struct), intent(in) :: Output
  type(acha_diag_struct), intent(inout), optional :: Diag

  !===============================================================================
  !  Local Variable Declaration
  !===============================================================================

  integer:: Elem_Idx
  integer:: Line_Idx
!  integer, private:: Sfc_Level_RTM
!  integer, private:: Tropo_Level_RTM

  real :: t_cld_h2o
  real :: z_cld_h2o
  real :: p_cld_h2o
  real :: t_cld_h2o_73
  real :: z_cld_h2o_73
  real :: p_cld_h2o_73
  real :: t_cld_co2
  real :: z_cld_co2
  real :: p_cld_co2
  integer :: overlap_flag
  real :: effec_emiss

  !--- initialize diagnostic output
  !if (present(Diag)) Diag%Array_1 = Missing_Value_Real
  if (present(Diag)) Diag%Array_2 = Missing_Value_Real
  if (present(Diag)) Diag%Array_3 = Missing_Value_Real

      t_cld_h2o = MISSING_VALUE_REAL
      z_cld_h2o = MISSING_VALUE_REAL
      p_cld_h2o = MISSING_VALUE_REAL
      t_cld_h2o_73 = MISSING_VALUE_REAL
      z_cld_h2o_73 = MISSING_VALUE_REAL
      p_cld_h2o_73 = MISSING_VALUE_REAL

      t_cld_co2 = MISSING_VALUE_REAL
      z_cld_co2 = MISSING_VALUE_REAL
      p_cld_co2 = MISSING_VALUE_REAL
      effec_emiss = MISSING_VALUE_REAL
    
   Line_loop: do Line_Idx = Line_Idx_min,Input%Number_of_Lines + Line_Idx_min - 1

    Element_Loop:   do Elem_Idx = 1, Input%Number_of_Elements

      overlap_flag = 0  !MISSING_VALUE_INTEGER

    !---- null profile pointers each time 
    call NULL_PIX_POINTERS(Input, ACHA_RTM_NWP)

    !- populate rtm data
    call ACHA_FETCH_PIXEL_NWP_RTM(Input, symbol, &
                                Elem_Idx,Line_Idx, ACHA_RTM_NWP)

      if      (Input%Cloud_Mask(Elem_Idx,Line_Idx) == symbol%CLEAR &
          .or. Input%Cloud_Mask(Elem_Idx,Line_Idx) == symbol%PROB_CLEAR &
          .or. Input%Invalid_Data_Mask(Elem_Idx,Line_Idx) == symbol%YES) then
          cycle
      end if

  Tropo_Level_RTM = ACHA_RTM_NWP%Tropo_Level
  Sfc_Level_RTM = ACHA_RTM_NWP%Sfc_Level

      if (Input%Chan_On_67um == symbol%YES) then

      if (0) then 
        call CO2IRW_CLOUD_HEIGHT ( &
            Input%Rad_11um(Elem_Idx,Line_Idx)  &                        !inp % sat % rad_ch31 &
           ,Acha_RTM_NWP%Black_Body_Rad_Prof_11um &        !inp % rtm % rad_ch31_bb_prof &
           ,Input%Rad_Clear_11um(Elem_Idx,Line_Idx)  &                        !inp % rtm % rad_ch31_atm_sfc &
           ,Input%Rad_67um(Elem_Idx,Line_Idx)  &                              !inp % sat % rad_ch27 &
           ,Acha_RTM_NWP%Black_Body_Rad_Prof_67um &        !inp % rtm % rad_ch27_bb_prof &
           ,Input%Rad_Clear_67um(Elem_Idx,Line_Idx)  &  !inp % rtm % rad_ch27_atm_sfc &
           ,Tropo_Level_RTM      &   !inp % rtm % tropo_lev &
           ,Sfc_Level_RTM        &   !inp % rtm % sfc_lev &
           ,ACHA_RTM_NWP % p_prof  &  !inp % rtm % p_prof &
           ,ACHA_RTM_NWP % t_prof  &  !inp % rtm % t_prof &
           ,ACHA_RTM_NWP % z_prof  &  !inp % rtm % z_prof &
           , p_cld_h2o &
           , t_cld_h2o &
           , z_cld_h2o )
        end if

        if (1) then
        call HEIGHT_H2O_CHANNEL ( &
            Input%Rad_11um(Elem_Idx,Line_Idx)  &                        !inp % sat % rad_ch31 &
           ,Acha_RTM_NWP%Black_Body_Rad_Prof_11um &        !inp % rtm % rad_ch31_bb_prof &
           ,Input%Rad_Clear_11um(Elem_Idx,Line_Idx)  & !inp % rtm % rad_ch31_atm_sfc &
           ,Input%Rad_67um(Elem_Idx,Line_Idx)  & !inp % sat % rad_ch27 &
           ,Acha_RTM_NWP%Black_Body_Rad_Prof_67um &        !inp % rtm % rad_ch27_bb_prof &
           ,Input%Rad_Clear_67um(Elem_Idx,Line_Idx)  & 
           ,Covar_Ch27_Ch31_5x5(Elem_Idx,Line_Idx)  &
           ,Tropo_Level_RTM      &   !inp % rtm % tropo_lev &
           ,Sfc_Level_RTM        &   !inp % rtm % sfc_lev &
           ,ACHA_RTM_NWP % t_prof  &  !inp % rtm % p_prof &
           ,ACHA_RTM_NWP % z_prof  &  !inp % rtm % t_prof &
           ,ACHA_RTM_NWP % p_prof  &  !inp % rtm % z_prof &
           , t_cld_h2o &
           , z_cld_h2o &
           , p_cld_h2o )
         end if

      end if

      if (Input%Chan_On_73um == symbol%YES) then
        call HEIGHT_H2O_CHANNEL ( &
            Input%Rad_11um(Elem_Idx,Line_Idx)  &                       
           ,Acha_RTM_NWP%Black_Body_Rad_Prof_11um &
           ,Input%Rad_Clear_11um(Elem_Idx,Line_Idx)  & 
           ,Input%Rad_73um(Elem_Idx,Line_Idx)  &
           ,Acha_RTM_NWP%Black_Body_Rad_Prof_73um &
           ,Input%Rad_Clear_73um(Elem_Idx,Line_Idx)  &
           ,Covar_Ch27_Ch31_5x5(Elem_Idx,Line_Idx)  &
           ,Tropo_Level_RTM      &   !inp % rtm % tropo_lev &
           ,Sfc_Level_RTM        &   !inp % rtm % sfc_lev &
           ,ACHA_RTM_NWP % t_prof  &  !inp % rtm % p_prof &
           ,ACHA_RTM_NWP % z_prof  &  !inp % rtm % t_prof &
           ,ACHA_RTM_NWP % p_prof  &  !inp % rtm % z_prof &
           , t_cld_h2o_73 &
           , z_cld_h2o_73 &
           , p_cld_h2o_73 )
         end if

      if (Input%Chan_On_133um == symbol%YES) then

        call CO2IRW_CLOUD_HEIGHT ( &
            Input%Rad_11um(Elem_Idx,Line_Idx)  &                        !inp % sat % rad_ch31 &
           ,Acha_RTM_NWP%Black_Body_Rad_Prof_11um &        !inp % rtm % rad_ch31_bb_prof &
           ,Input%Rad_Clear_11um(Elem_Idx,Line_Idx)  &                        !inp % rtm % rad_ch31_atm_sfc & 
           ,Input%Rad_133um(Elem_Idx,Line_Idx)  &                              !inp % sat % rad_ch27 &
           ,Acha_RTM_NWP%Black_Body_Rad_Prof_133um &        !inp % rtm % rad_ch27_bb_prof &
           ,Input%Rad_Clear_133um(Elem_Idx,Line_Idx)  &  !inp % rtm % rad_ch27_atm_sfc &
           ,Tropo_Level_RTM      &   !inp % rtm % tropo_lev &
           ,Sfc_Level_RTM        &   !inp % rtm % sfc_lev &
           ,ACHA_RTM_NWP % p_prof  &  !inp % rtm % p_prof &
           ,ACHA_RTM_NWP % t_prof  &  !inp % rtm % t_prof &
           ,ACHA_RTM_NWP % z_prof  &  !inp % rtm % z_prof &
           , p_cld_co2 &
           , t_cld_co2 &
           , z_cld_co2 )

      end if

      effec_emiss = EMISSIVITY(Input%Rad_11um(Elem_Idx,Line_Idx),&
                        Input%Rad_Clear_11um(Elem_Idx,Line_Idx),  &
                        Acha_RTM_NWP%Black_Body_Rad_Prof_11um(Tropo_Level_RTM))

     !emissivity threshold
!      if (Output%Ec(Elem_Idx,Line_Idx) < 0.95 .and. &
      if (effec_emiss < 0.9 .and. &
          p_cld_h2o > 0 .and. p_cld_co2 > 0 .and. &
          p_cld_h2o_73 > 0 .and. abs(p_cld_co2-p_cld_h2o_73) > 50. .and. &
          abs(p_cld_co2-p_cld_h2o) > 50. ) then
!          print *,p_cld_h2o, p_cld_co2
          overlap_flag = 1
      end if
     !

     !      Diag%Array_1(Elem_Idx,Line_Idx) = Input%Rad_Clear_11um(Elem_Idx,Line_Idx) !p_cld
     !      Diag%Array_1(Elem_Idx,Line_Idx) = p_cld_co2  !overlap_flag 
           Diag%Array_2(Elem_Idx,Line_Idx) = p_cld_h2o
           Diag%Array_3(Elem_Idx,Line_Idx) = overlap_flag  !p_cld_h2o_73

     end do Element_Loop

    end do Line_Loop

 !---- null profile pointers each time 
 call NULL_PIX_POINTERS(Input, ACHA_RTM_NWP)

  end subroutine ovelap_co2_h2o_test

!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------

!====================================================================
subroutine  CO2IRW_CLOUD_HEIGHT(Rad_11um, &
                             Rad_11um_BB_Profile, &
                             Rad_11um_Clear, &
                             Rad_CO2, &
                             Rad_CO2_BB_Profile,  &
                             Rad_CO2_Clear,  &
                             Tropo_Level, &
                             Sfc_Level, &
                             P_Prof, &
                             T_Prof, &
                             Z_Prof, &
                             Pc,  &
                             Tc,  &
                             Zc)

  real, intent(in):: Rad_11um
  real, intent(in):: Rad_CO2
  real, intent(in), dimension(:):: Rad_11um_BB_Profile
  real, intent(in):: Rad_11um_Clear
  real, intent(in):: Rad_CO2_Clear
  real, intent(in), dimension(:):: Rad_CO2_BB_Profile
  integer, intent(in):: Sfc_Level
  integer, intent(in):: Tropo_Level
  real, intent(in), dimension(:):: P_Prof
  real, intent(in), dimension(:):: T_Prof
  real, intent(in), dimension(:):: Z_Prof
  real, intent(out) :: Pc
  real, intent(out) :: Tc
  real, intent(out) :: Zc

  real:: Rad_CO2_BB_Prediction
  real:: Slope
  real:: Intercept
  real:: Denominator
  integer:: ilev
  integer:: ilev_co2

  real, parameter:: Rad_11um_Thresh = 2.0
  real, parameter:: Rad_CO2_Thresh = 2.0

  !--- initialize
  Pc = MISSING_VALUE_REAL
  Tc = MISSING_VALUE_REAL
  Zc = MISSING_VALUE_REAL

  !--- determine if a solution should be attempted
  if (Rad_11um_Clear - Rad_11um < Rad_11um_Thresh) return
  if (Rad_11um == MISSING_VALUE_REAL) return
  if (Rad_11um_Clear == MISSING_VALUE_REAL) return

  if (Rad_CO2_Clear - Rad_CO2 < Rad_CO2_Thresh) return
  if (Rad_CO2 == MISSING_VALUE_REAL) return
  if (Rad_CO2_Clear == MISSING_VALUE_REAL) return

 !--- attempt a solution

 !--- colder than tropo
 if (Rad_11um < Rad_11um_BB_Profile(Tropo_Level)) then

     ilev_co2 = Tropo_Level

 else   !if not, attempt solution

     !--- determine linear regress of co2 (y)  as a function of window (x)
      Denominator =  Rad_11um - Rad_11um_Clear

      if (Denominator < 0.0) then
             Slope = (Rad_CO2 - Rad_CO2_Clear) / (Denominator)
             Intercept = Rad_CO2 - Slope*Rad_11um
      else
            return
      endif

      !--- brute force solution
      ilev_co2 = 0

      do ilev = Tropo_Level+1, Sfc_Level
          Rad_CO2_BB_Prediction = Slope*Rad_11um_BB_Profile(ilev) + Intercept

          if (Rad_CO2_BB_Prediction < 0) cycle

          if ((Rad_CO2_BB_Prediction > Rad_CO2_BB_Profile(ilev-1)) .and. &
               (Rad_CO2_BB_Prediction <= Rad_CO2_BB_Profile(ilev))) then
               ilev_co2 = ilev
               exit
          endif

      enddo

 endif    !tropopause check

 !--- adjust back to full Rtm profile indices
 if (ilev_co2 > 0) then
       Pc = P_Prof(ilev_co2)
       Tc = T_Prof(ilev_co2)
       Zc = Z_Prof(ilev_co2)
 endif

end subroutine CO2IRW_CLOUD_HEIGHT

   !====================================================================
   ! Function Name: HEIGHT_H2O_CHANNEL
   !
   ! Function: estimate the cloud temperature/height/pressure
   !
   ! Description: Use the 11um and 6.7um obs and the RTM cloud BB profiles
   !              to perform h2o intercept on a pixel level. Filters
   !              restrict this to high clouds only
   !              
   ! Dependencies: 
   !
   ! Restrictions: 
   !
   ! Reference: 
   !
   ! Author: Andrew Heidinger, NOAA/NESDIS
   !
   !====================================================================   
   subroutine HEIGHT_H2O_CHANNEL ( &
           rad_11um &
         , rad_11um_rtm_prof &
         , rad_11um_clear &
         , rad_h2o &
         , rad_h2o_rtm_prof &
         , rad_h2o_clear &
         , Covar_h2o_11 &
         , tropo_lev &
         , sfc_lev &
         , t_prof &
         , z_prof &
         , p_prof &
         , t_cld &
         , z_cld &
         , p_cld)

      implicit none

      real, intent(in) :: rad_11um
      real, intent(in) :: rad_11um_clear
      real , intent(in) :: rad_11um_rtm_prof(:)
      real , intent(in) :: rad_h2o_rtm_prof(:)
      real, intent(in) :: rad_h2o
      real, intent(in) :: rad_h2o_clear
      real, intent(in) :: Covar_h2o_11
      integer, intent(in) :: tropo_lev
      integer, intent(in) :: sfc_lev
      real, intent(in) :: t_prof(:)
      real, intent(in) :: z_prof(:)
      real, intent(in), dimension(:) :: p_prof
      real, intent(out) :: t_cld
      real, intent(out) :: z_cld
      real, intent(out) :: p_cld

      real, parameter :: RAD_11UM_THRESH = 2.0
      real, parameter :: RAD_67UM_THRESH = 0.25
      real, parameter :: BT_CH27_CH31_COVAR_CIRRUS_THRESH = 1.0

      integer :: cld_lev
      integer :: idx_lev

      real :: denominator
      real :: slope
      real :: intercept

      real :: rad_H2O_BB_prediction

      ! -------------------------------            
      t_cld = MISSING_VALUE_REAL
      z_cld = MISSING_VALUE_REAL
      p_cld = MISSING_VALUE_REAL

      ! some tests
      if ( rad_h2o < 0 ) return

      if ( rad_h2o_clear - rad_h2o < RAD_67UM_THRESH ) return

      if ( rad_11um_clear - rad_11um < RAD_11UM_THRESH ) return
      if ( Covar_h2o_11 /= MISSING_VALUE_REAL .and.  Covar_h2o_11 < BT_CH27_CH31_COVAR_CIRRUS_THRESH ) return

      ! - colder than tropopause
      if ( rad_11um < rad_11um_rtm_prof ( tropo_lev ) ) then
         cld_lev = tropo_lev
      else
!--- determine linear regress of h2o (y)  as a function of window (x)
         Denominator =  rad_11um - rad_11um_clear

         slope = (Rad_h2o - Rad_H2o_Clear) / (Denominator)
         intercept = Rad_h2o - Slope * Rad_11um

         cld_lev = 0

         do idx_lev = tropo_lev + 1, sfc_lev

            rad_H2O_BB_prediction = slope * rad_11um_rtm_prof (idx_lev) + Intercept

            if ( rad_H2O_BB_prediction < 0 ) cycle

            if ((rad_H2O_BB_Prediction > rad_h2o_rtm_prof(idx_lev-1)) .and. &
               ( rad_H2O_BB_Prediction <= rad_h2o_rtm_prof(idx_lev))) then
               cld_lev = idx_lev
               exit
          endif

         end do

      end if

      !--- adjust back to full Rtm profile indices
      if (cld_lev > 0 ) then
         t_cld = t_prof( cld_lev )
         z_cld = z_prof( cld_lev )
         p_cld = p_prof( cld_lev )
      end if

   end subroutine HEIGHT_H2O_CHANNEL


!------------------------------------------------------------------------------
! Null Pixel Level Pointers 
!------------------------------------------------------------------------------
subroutine NULL_PIX_POINTERS(Input, ACHA_RTM_NWP)

   type(acha_input_struct), intent(inout) :: Input
   type(acha_rtm_nwp_struct), intent(inout) :: ACHA_RTM_NWP

   ACHA_RTM_NWP%T_Prof => NULL()
   ACHA_RTM_NWP%T_Prof_1 => NULL()
   ACHA_RTM_NWP%T_Prof_2 => NULL()
   ACHA_RTM_NWP%T_Prof_3 => NULL()

   ACHA_RTM_NWP%Z_Prof => NULL()
   ACHA_RTM_NWP%Z_Prof_1 => NULL()
   ACHA_RTM_NWP%Z_Prof_2 => NULL()
   ACHA_RTM_NWP%Z_Prof_3 => NULL()

   if (Input%Chan_On_67um == symbol%YES) then
     ACHA_RTM_NWP%Atm_Rad_Prof_67um =>  NULL()
     ACHA_RTM_NWP%Atm_Trans_Prof_67um =>  NULL()
     ACHA_RTM_NWP%Black_Body_Rad_Prof_67um => NULL()
   endif
   if (Input%Chan_On_85um == symbol%YES) then
     ACHA_RTM_NWP%Atm_Rad_Prof_85um =>  NULL()
     ACHA_RTM_NWP%Atm_Trans_Prof_85um =>  NULL()
   endif

   if (Input%Chan_On_11um == symbol%YES) then
      ACHA_RTM_NWP%Atm_Rad_Prof_11um => NULL()
      ACHA_RTM_NWP%Atm_Trans_Prof_11um => NULL()
      ACHA_RTM_NWP%Black_Body_Rad_Prof_11um => NULL()
   endif
   if (Input%Chan_On_12um == symbol%YES) then
      ACHA_RTM_NWP%Atm_Rad_Prof_12um => NULL()
      ACHA_RTM_NWP%Atm_Trans_Prof_12um => NULL()
   endif
   if (Input%Chan_On_133um == symbol%YES) then
      ACHA_RTM_NWP%Atm_Rad_Prof_133um => NULL()
      ACHA_RTM_NWP%Atm_Trans_Prof_133um => NULL()
   endif


end subroutine NULL_PIX_POINTERS

   !====================================================================
   ! Function Name: EMISSIVITY
   !
   ! Function:
   !  Computes the  effective emissivity
   !
   ! Input:  Rad_Toa - channel radiance at top of atmosphere(toa)
   !         Rad_Clear_Tau - channel radiance at toa for clear skies
   !         Rad_Cloud_BB_Toa - channel radiance at TOA if cloud were a
   !         Black-Body
   !
   ! Output: Emiss - the effective cloud emissivity
   !
   !====================================================================
   function EMISSIVITY(Radiance_Toa, Radiance_Clear_Toa, Radiance_Cloud_BB_Toa) result(Emiss)
      real, intent(in) :: Radiance_Toa
      real, intent(in) :: Radiance_Clear_Toa
      real, intent(in) :: Radiance_Cloud_BB_Toa
      real :: Emiss

      Emiss = Missing_Value_Real

      if (Radiance_Cloud_BB_Toa /= Radiance_Clear_Toa) then
          Emiss = (Radiance_Toa - Radiance_Clear_Toa) / &
            (Radiance_Cloud_BB_Toa - Radiance_Clear_Toa)
       end if

      return

   end function EMISSIVITY


!----------------------------------------------------------------------
! End of Module
!----------------------------------------------------------------------

end module eumetset_overlap
