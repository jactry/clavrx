module NB_CLOUD_MASK_SOLAR_RTM

 use NB_CLOUD_MASK_SERVICES, only: MISSING_VALUE_REAL4

 implicit none

 public:: CLEAR_SKY_TOA_RTM_065UM

 !---------------------------------------------------------------------
 !
 !---------------------------------------------------------------------
 integer, parameter, private:: int1 = selected_int_kind(1)
 integer, parameter, private:: int2 = selected_int_kind(3)
 integer, parameter, private:: int4 = selected_int_kind(8)
 integer, parameter, private:: int8 = selected_int_kind(10)
 integer, parameter, private:: real4 = selected_real_kind(6,37)
 integer, parameter, private:: real8 = selected_real_kind(15,307)
 integer, parameter, private:: ipre = real4

 real, parameter, private:: Rayleigh_Optical_Depth_065um = 0.055
 real, parameter, private:: Aerosol_Optical_Depth_065um = 0.12
 real, parameter, private:: Aerosol_Single_Scatter_Albedo_065um = 0.8
 real, parameter, private:: Aerosol_Asymmetry_Parameter = 0.6
 real, dimension(3), parameter, private:: OD_ozone_coef = (/0.000566454,8.25224e-05,1.94007e-08/)
 real, dimension(3), parameter, private:: OD_h2o_coef = (/  0.000044758, 0.00264790,-0.0000713698/)
 real, parameter, private:: dtor = 0.0174533

 real, dimension(0:13),parameter, private:: Ch1_Sfc_Alb_Umd = &
                                                      (/0.0500, & !0
                                                        0.0406,  &  !1
                                                        0.0444,  &  !2
                                                        0.0406,  &  !3
                                                        0.0542,  &  !4
                                                        0.0423,  &  !5
                                                        0.0423,  &  !6
                                                        0.0576,  &  !7
                                                        0.1242,  &  !8
                                                        0.1781,  &  !9
                                                        0.0788,  &  !10
                                                        0.0677,  &  !11
                                                        0.3124,  &  !12
                                                        0.0918/)    !13

  real, dimension(0:13),parameter, private:: Ch1_Snow_Sfc_Alb_Umd = &
                                                           (/0.66, & !0
                                                             0.23 , &  !1
                                                             0.50,  &  !2
                                                             0.42,  &  !3
                                                             0.25,  &  !4
                                                             0.21,  &  !5
                                                             0.49,  &  !6
                                                             0.59,  &  !7
                                                             0.72,  &  !8
                                                             0.78,  &  !9
                                                             0.70,  &  !10
                                                             0.72,  &  !11
                                                             0.76,  &  !12
                                                             0.65/)    !13



 contains

!====================================================================
! Function Name: CLEAR_SKY_TOA_RTM_065UM
!
! Function:
!   Computes the single scater and aerosol reflectance assuming that the
!   gas is mixed in with scattering
!
! Description: Apply a simple single scattering RTM to compute
!              top-of-atmosphere clear-sky reflectance since the
!              SAPF is unable to do at this time
!  
!
! Calling Sequence:
! call  CLEAR_SKY_TOA_RTM_065UM(TPW, &
!                               TOzone, &
!                               Scat_Zen, &
!                               Sat_Zen, &
!                               Sol_Zen, &
!                               Surface_Reflectance, &
!                               Sfc_Type, &
!                               Snow_Class, &
!                               Toa_Clear_Sky_Refl)
!
! Inputs:
!    TPW - Total Precipitable Water
!    TOzone - Total Ozone
!    Scat_Zen - Scattering Angle
!    Sat_Zen - Sensor Zenith Angle
!    Sol_Zen - Solar Zenith Angle
!    Surface_Reflectance - Surface Relfectance 
!    Sfc_Type - UMD Vegetation Type
!    Snow_Class - Snow Classification
!
! Outputs: 
!    Toa_Clear_Sky_Refl - Top of Atmosphere Clear-sky Reflectance at 0.65 um
!
! Dependencies:  None
!
! Restrictions:  None
!
!====================================================================
 subroutine CLEAR_SKY_TOA_RTM_065UM(Bad_Pixel_Mask, &
                                    TPW, &
                                    TOzone, &
                                    Scat_Zen, &
                                    Sat_Zen, &
                                    Sol_Zen, &
                                    Surface_Reflectance, &
                                    Sfc_Type, &
                                    Snow_Class, &
                                    Toa_Clear_Sky_Refl)

   integer(kind=int1), dimension(:,:), intent(in):: Bad_Pixel_Mask
   real, dimension(:,:), intent(in):: TPW
   real, dimension(:,:), intent(in):: TOzone
   real, dimension(:,:), intent(in):: Scat_Zen
   real, dimension(:,:), intent(in):: Sat_Zen
   real, dimension(:,:), intent(in):: Sol_Zen
   real, dimension(:,:), intent(in):: Surface_Reflectance
   integer(kind=int1), dimension(:,:), intent(in):: Sfc_Type
   integer(kind=int1), dimension(:,:), intent(in):: Snow_Class
   real, dimension(:,:), intent(out):: TOA_Clear_Sky_Refl
   real:: Cos_Scat_Zen
   real:: Cos_Sat_Zen
   real:: Cos_Sol_Zen
   real:: Sfc_Alb_View
   real:: Sfc_Alb_Sun
   real:: Transmission_Sing_Scat
   real:: Air_Mass_Factor
   real:: Aero_Phase_Funct
   real:: Ray_Phase_Funct
   real:: OD_Gas
   real:: OD_Total
   real:: OD_Scat_Total
   real:: OD_Iso_Total
   real:: Trans_Iso_Total_View
   real:: Trans_Iso_Total_Sun
   real:: OD_Iso_Scat_Total
   real:: Eff_Phase_Funct
   real:: Single_Scat_Alb
   real:: Refl_Sing_Scat_a
   real:: Refl_Sing_Scat_b
   real:: Refl_Sing_Scat_c
   real:: Refl_Sing_Scat
   real:: OD_ozone
   real:: OD_h2o
   integer:: Elem_Idx,Line_Idx, Num_Elem, Num_Line

   Num_Elem = size(Sol_Zen,1)
   Num_Line = size(Sol_Zen,2)

   !--- initialize to missing
   TOA_Clear_Sky_Refl = MISSING_VALUE_REAL4

   do Elem_Idx = 1, Num_Elem
    do Line_Idx = 1, Num_Line

    !--- skip if bad data
    if (Bad_Pixel_Mask(ELem_Idx,Line_Idx) == 1) cycle

    !--- skip if night
    if (Sol_Zen(ELem_Idx,Line_Idx) > 90.0) cycle

    !--- compute cosine of scattering angle
    Cos_Scat_Zen = cos(Scat_Zen(Elem_Idx,Line_Idx) * dtor)
    Cos_Sat_Zen = cos(Sat_Zen(Elem_Idx,Line_Idx) * dtor)
    Cos_Sol_Zen = cos(Sol_Zen(Elem_Idx,Line_Idx) * dtor)

    !-------------------------------------------------------------------------------
    !  Surface Albedo
    !-------------------------------------------------------------------------------
    if (Sfc_Type(Elem_Idx,Line_Idx) < 0) cycle   !check for missing sfc type

    !--- set to surface-type default value
    Sfc_Alb_Sun = Ch1_Sfc_Alb_Umd(Sfc_Type(Elem_Idx,Line_Idx))

    !--- if land and white sky is available, use it
    if (Sfc_Type(Elem_Idx,Line_Idx) > 0 .and. Surface_Reflectance(Elem_Idx,Line_Idx) /= MISSING_VALUE_REAL4) then
       Sfc_Alb_Sun =  Surface_Reflectance(Elem_Idx,Line_Idx) / 100.0    ! must be between 0 and 1
    endif

    !--- if snow, use precomputed value based on surface type
    if (Snow_Class(Elem_Idx,Line_Idx) > 1) then
         Sfc_Alb_Sun = Ch1_Snow_Sfc_Alb_Umd(Sfc_Type(Elem_Idx,Line_Idx))
    endif

    Sfc_Alb_View =  Sfc_Alb_Sun

    !-------------------------------------------------------------------------------
    !  Gas Optical Depths
    !-------------------------------------------------------------------------------
    OD_H2O = OD_h2o_coef(1) + OD_h2o_coef(2)*TPW(Elem_Idx,Line_Idx) +  &
             OD_h2o_coef(3)*TPW(Elem_Idx,Line_Idx)**2
   
    OD_Ozone = OD_ozone_coef(1) + OD_ozone_coef(2) * TOzone(Elem_Idx,Line_Idx) + &
               OD_ozone_coef(3) * TOzone(Elem_Idx,Line_Idx)**2

    OD_Gas = OD_h2o + OD_ozone

    !-- compute Rayleigh phase function
    Air_Mass_Factor = 1.0 / Cos_Sat_Zen + 1.0 / Cos_Sol_Zen
   
    Ray_Phase_Funct = 0.75 * (1.0 + Cos_Scat_Zen**2)

    !--- compute total transmission
    OD_Total = Aerosol_Optical_Depth_065um + Rayleigh_Optical_Depth_065um + OD_Gas
   
    Transmission_Sing_Scat = exp(-OD_Total * Air_Mass_Factor)

    OD_Iso_Total = (1.0 - Aerosol_Asymmetry_Parameter) * Aerosol_Optical_Depth_065um + &
                    Rayleigh_Optical_Depth_065um + OD_Gas
   
    Trans_Iso_Total_View = exp(-OD_Iso_Total / Cos_Sat_Zen)
   
    Trans_Iso_Total_Sun = exp(-OD_Iso_Total / Cos_Sol_Zen)

    !--- compute total scattering optical depth
    OD_Scat_Total = Aerosol_Single_Scatter_Albedo_065um *&
                   Aerosol_Optical_Depth_065um + Rayleigh_Optical_Depth_065um
   
    OD_Iso_Scat_Total = Aerosol_Single_Scatter_Albedo_065um * (1.0 - Aerosol_Asymmetry_Parameter) * &
                    Aerosol_Optical_Depth_065um + Rayleigh_Optical_Depth_065um

    !--- single scatter albedo
    Single_Scat_Alb = (Aerosol_Single_Scatter_Albedo_065um * Aerosol_Optical_Depth_065um + &
                     Rayleigh_Optical_Depth_065um) / ( OD_Total )
 
    !aerosol phase function (Henyey-Greenstein)
    Aero_Phase_Funct = (1.0 - Aerosol_Asymmetry_Parameter**2) / &
                    ( (1.0 + Aerosol_Asymmetry_Parameter**2 -  &
                       2.0 * Aerosol_Asymmetry_Parameter*Cos_Scat_Zen)**(1.5) )

    !--- compute effective phase function
    Eff_Phase_Funct = (Aerosol_Single_Scatter_Albedo_065um * Aerosol_Optical_Depth_065um * &
    Aero_Phase_Funct + Rayleigh_Optical_Depth_065um * Ray_Phase_Funct) / (OD_Scat_Total)

    !--- compute single scatter reflectance (0-100%)
    Refl_Sing_Scat_a = Single_Scat_Alb * Eff_Phase_Funct / (4.0 * Air_Mass_Factor * &
                        Cos_Sat_Zen * Cos_Sol_Zen) * (1.0 - &
                        Transmission_Sing_Scat )
 
    Refl_Sing_Scat_b = (OD_Iso_Scat_Total / (2.0*Cos_Sol_Zen)) * &
                        Trans_Iso_Total_View * Sfc_Alb_View

    Refl_Sing_Scat_c = (OD_Iso_Scat_Total / (2.0 * Cos_Sat_Zen)) * &
                        Trans_Iso_Total_Sun * Sfc_Alb_Sun

    Refl_Sing_Scat = 100.0 * (Refl_Sing_Scat_a + Refl_Sing_Scat_b + Refl_Sing_Scat_c)

    TOA_Clear_Sky_Refl(Elem_Idx,Line_Idx) =  Refl_Sing_Scat + 100.0*Transmission_Sing_Scat*Sfc_Alb_View

    enddo
   enddo

 end subroutine CLEAR_SKY_TOA_RTM_065UM

!----------------------------------------------
end module NB_CLOUD_MASK_SOLAR_RTM
