! $Id:$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: dcomp_derived_products_module.f90 (src)
!       DCOMP_DERIVED_PRODUCTS_MODULE (program)
!
! PURPOSE: this modules hold subroutines that derive additional products
!          primarily from DCOMP
!
! DESCRIPTION: 
!
! AUTHORS:
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
! Public routines used in this MODULE:
! COMPUTE_CLOUD_WATER_PATH - compute cloud water path
! COMPUTE_PRECIPITATION - compute precipitation using KNMI approach
! COMPUTE_ADIABATIC_CLOUD_PROPS - compute precipitation using KNMI approach
!
!--------------------------------------------------------------------------------------
MODULE DCOMP_DERIVED_PRODUCTS_MODULE
 use CONSTANTS
 use ALGORITHM_CONSTANTS
 use PIXEL_COMMON
 use NWP_COMMON

 implicit none
 public:: COMPUTE_CLOUD_WATER_PATH, &
          COMPUTE_PRECIPITATION, &
          COMPUTE_ADIABATIC_CLOUD_PROPS, &
          COMPUTE_DCOMP_INSOLATION

  contains

!-----------------------------------------------------------
! compute cloud water path from the optical depth
! and particle size from the dcomp algorithm
!
! The layer values are computed assuming a linear variation
! in cloud water path from the top to the base of the cloud.
! Note CWP = CWP_Ice_Layer + CWP_Water_Layer and 
!      CWP_Scwater is a component of the Water_Layer
! 
!-----------------------------------------------------------
subroutine COMPUTE_CLOUD_WATER_PATH(jmin,jmax)

  integer, intent(in):: jmin
  integer, intent(in):: jmax

  integer:: Elem_Idx
  integer:: Line_Idx
  integer:: Iphase

  real(kind=real4), parameter:: Rho_Water = 1.0    !g/m^3
  real(kind=real4), parameter:: Rho_Ice = 0.917    !g/m^3
  integer:: Lat_NWP_Idx
  integer:: Lon_NWP_Idx
  real:: Cloud_Geometrical_Thickness
  real:: Ice_Layer_Fraction
  real:: Water_Layer_Fraction
  real:: Scwater_Layer_Fraction
  real:: Tau
  real:: Reff

  Cwp_Dcomp = Missing_Value_Real4
  Iwp_Dcomp = Missing_Value_Real4
  Lwp_Dcomp = Missing_Value_Real4
  Cwp_Ice_Layer_Dcomp = Missing_Value_Real4
  Cwp_Water_Layer_Dcomp = Missing_Value_Real4
  Cwp_Scwater_Layer_Dcomp = Missing_Value_Real4
  Tau = Missing_Value_Real4
  Reff = Missing_Value_Real4

  line_loop: DO Line_Idx = jmin, jmax - jmin + 1
    element_loop: DO Elem_Idx = 1, Image%Number_Of_Elements

     !--- assign optical depth and particle size
     if (Geo%Solzen(Elem_Idx,Line_Idx) < 90.0) then 
       Tau = Tau_Dcomp(Elem_Idx,Line_Idx)
       Reff = Reff_Dcomp(Elem_Idx,Line_Idx)
     else
       Tau = Tau_Nlcomp(Elem_Idx,Line_Idx)
       Reff = Reff_Nlcomp(Elem_Idx,Line_Idx)
     endif

     if (Tau == Missing_Value_Real4) cycle
     if (Reff == Missing_Value_Real4) cycle

     !------------------------------------------------
     ! determine phase from cloud type
     ! -1 = undetermined, 0 = water, 1 = ice
     !------------------------------------------------
      Iphase = -1
      if (Cld_Type(Elem_Idx,Line_Idx) == sym%CLEAR_TYPE) then 
              Iphase = -1
      elseif (Cld_Type(Elem_Idx,Line_Idx) == sym%PROB_CLEAR_TYPE)  then
              Iphase = -1 
      elseif (Cld_Type(Elem_Idx,Line_Idx) == sym%FOG_TYPE)  then
              Iphase = 0
      elseif (Cld_Type(Elem_Idx,Line_Idx) == sym%WATER_TYPE)  then
              Iphase = 0
      elseif (Cld_Type(Elem_Idx,Line_Idx) == sym%SUPERCOOLED_TYPE)  then
              Iphase = 0
      elseif (Cld_Type(Elem_Idx,Line_Idx) == sym%MIXED_TYPE)  then
              Iphase = 0
      elseif (Cld_Type(Elem_Idx,Line_Idx) == sym%OPAQUE_ICE_TYPE)  then
              Iphase = 1
      elseif (Cld_Type(Elem_Idx,Line_Idx) == sym%TICE_TYPE)  then
              Iphase = 1
      elseif (Cld_Type(Elem_Idx,Line_Idx) == sym%CIRRUS_TYPE)  then
              Iphase = 1
      elseif (Cld_Type(Elem_Idx,Line_Idx) == sym%OVERLAP_TYPE)  then
              Iphase = 1
      elseif (Cld_Type(Elem_Idx,Line_Idx) == sym%OVERSHOOTING_TYPE)  then
              Iphase = 1
      elseif (Cld_Type(Elem_Idx,Line_Idx) == sym%UNKNOWN_TYPE)  then
              Iphase = -1
      elseif (Cld_Type(Elem_Idx,Line_Idx) == sym%DUST_TYPE)  then
              Iphase = -1
      elseif (Cld_Type(Elem_Idx,Line_Idx) == sym%SMOKE_TYPE)  then
              Iphase = -1
     endif 

     !--- check conditions where this calc should be skipped
     if (Iphase == -1) cycle

     !--- compute cloud water path
     if (Iphase == 0) then
      Cwp_Dcomp(Elem_Idx,Line_Idx) = 0.55*Tau*Reff*Rho_Water
      Lwp_Dcomp(Elem_Idx,Line_Idx) = 0.55*Tau*Reff*Rho_Water
     else
      Cwp_Dcomp(Elem_Idx,Line_Idx) = 0.667*Tau*Reff*Rho_Ice
      Iwp_Dcomp(Elem_Idx,Line_Idx) = 0.667*Tau*Reff*Rho_Ice
     endif

     !--- Partition into Ice, Water and Scwater Layers
     Lon_NWP_Idx = I_Nwp(Elem_Idx,Line_Idx)
     Lat_NWP_Idx = J_Nwp(Elem_Idx,Line_Idx)

     !--- skip if invalid nwp indices
     if (Lat_Nwp_Idx <= 0 .or. Lon_Nwp_Idx <= 0) cycle

     Cloud_Geometrical_Thickness = ACHA%Zc_Top(Elem_Idx,Line_Idx) - ACHA%Zc_Base(Elem_Idx,Line_Idx)
     !--- skip if failed cloud boundares
     if (Cloud_Geometrical_Thickness <= 0.00 .or. ACHA%Zc_Top(Elem_Idx,Line_Idx) <= 0.00) cycle

     Ice_Layer_Fraction = 0.0
     Water_Layer_Fraction = 0.0
     Scwater_Layer_Fraction = 0.0

     if (ACHA%Zc_Base(Elem_Idx,Line_Idx) >= Upper_Limit_Water_Height_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx)) then

         Ice_Layer_Fraction  = 1.0
         Water_Layer_Fraction  = 0.0
         Scwater_Layer_Fraction  = 0.0

     else

         Ice_Layer_Fraction = (ACHA%Zc_Top(Elem_Idx,Line_Idx) - Upper_Limit_Water_Height_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx)) / &
                               Cloud_Geometrical_Thickness

     endif

     if (Ice_Layer_Fraction /= 1.0) then

         if (ACHA%Zc_Top(Elem_Idx,Line_Idx) < Upper_Limit_Water_Height_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx)) then

           Water_Layer_Fraction = 1.0

         else
           Water_Layer_Fraction = (Upper_Limit_Water_Height_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx)-ACHA%Zc_Base(Elem_Idx,Line_Idx)) / &
                                  Cloud_Geometrical_Thickness
         endif

         if ((ACHA%Zc_Top(Elem_Idx,Line_Idx) > Upper_Limit_Water_Height_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx)) .and. &
             (ACHA%Zc_Base(Elem_Idx,Line_Idx) < Freezing_Level_Height_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx))) then

            Scwater_Layer_Fraction = (Upper_Limit_Water_Height_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx) - &
                                     Freezing_Level_Height_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx)) / &
                                     Cloud_Geometrical_Thickness

         elseif ((ACHA%Zc_Top(Elem_Idx,Line_Idx) > Upper_Limit_Water_Height_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx)) .and. &
                   (ACHA%Zc_Base(Elem_Idx,Line_Idx) < Upper_Limit_Water_Height_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx))) then
  
              Scwater_Layer_Fraction = (Upper_Limit_Water_Height_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx) - &
                                       ACHA%Zc_Base(Elem_Idx,Line_Idx)) / &
                                       Cloud_Geometrical_Thickness

         elseif ((ACHA%Zc_Top(Elem_Idx,Line_Idx) > Freezing_Level_Height_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx)) .and. &
                   (ACHA%Zc_Base(Elem_Idx,Line_Idx) < Freezing_Level_Height_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx))) then

              Scwater_Layer_Fraction = (ACHA%Zc_Top(Elem_Idx,Line_Idx) - Freezing_Level_Height_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx)) / &
                                        Cloud_Geometrical_Thickness

          endif

     endif

     Ice_Layer_Fraction = max(0.0,Ice_Layer_Fraction)
     Water_Layer_Fraction = max(0.0,Water_Layer_Fraction)
     Scwater_Layer_Fraction = max(0.0,Scwater_Layer_Fraction)

     Cwp_Ice_Layer_Dcomp(Elem_Idx,Line_Idx) = Ice_Layer_Fraction * Cwp_Dcomp(Elem_Idx,Line_Idx)
     Cwp_Water_Layer_Dcomp(Elem_Idx,Line_Idx) = Water_Layer_Fraction * Cwp_Dcomp(Elem_Idx,Line_Idx)
     Cwp_Scwater_Layer_Dcomp(Elem_Idx,Line_Idx) = Scwater_Layer_Fraction * Cwp_Dcomp(Elem_Idx,Line_Idx)

    enddo element_loop
  enddo line_loop


end subroutine COMPUTE_CLOUD_WATER_PATH
!-----------------------------------------------------------------------------
!--- compute Number concentration and Geometrical Height
!-----------------------------------------------------------------------------
subroutine COMPUTE_ADIABATIC_CLOUD_PROPS(Line_Idx_Min,Num_Lines)
  integer, intent(in):: Line_Idx_Min
  integer, intent(in):: Num_Lines

  !--- local parameters used in lwp, iwp, N and H computation
  real(kind=real4), parameter:: Rho_Water = 1000.0    !kg/m^3
  real(kind=real4), parameter:: Q_Eff_Sca  = 2.0
  real(kind=real4), parameter:: Drop_Dis_Width  = 0.8
  real(kind=real4):: Condensation_Rate
  real(kind=real4):: Water_Path_Cloud
  integer:: Elem_Idx, Line_Idx, Elem_Idx_Min, Elem_Idx_Max, Line_Idx_Max, Num_Elements

  real(kind=real4)::  T_Cloud, Reff_Cloud, Tau_Cloud

  Elem_Idx_Min = 1
  Num_Elements = Image%Number_Of_Elements
  Elem_Idx_Max = Elem_Idx_Min + Num_Elements - 1
  Line_Idx_Max = Line_Idx_Min + Num_Lines - 1

  Hcld_Dcomp = Missing_Value_Real4
  Cdnc_Dcomp = Missing_Value_Real4

  line_loop: do Line_Idx = Line_Idx_Min, Line_Idx_Max

    element_loop: do Elem_Idx = Elem_Idx_Min, Elem_Idx_Max

      !--- skip bad pixels
      if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) cycle

      !--- skip non cloud pixels
      if (Cld_Type(Elem_Idx,Line_Idx) == sym%CLEAR_TYPE .or. &
          Cld_Type(Elem_Idx,Line_Idx) == sym%PROB_CLEAR_TYPE) then
          Hcld_Dcomp(Elem_Idx,Line_Idx) = 0.0
          Cdnc_Dcomp(Elem_Idx,Line_Idx) = 0.0
          cycle
      endif

      !--- skip ice
      if (Cld_Type(Elem_Idx,Line_Idx) /= sym%FOG_TYPE .and. &
          Cld_Type(Elem_Idx,Line_Idx) /= sym%WATER_TYPE .and. &
          Cld_Type(Elem_Idx,Line_Idx) /= sym%SUPERCOOLED_TYPE) then
          cycle
      endif

      !--- make local aliases of global variables
      Water_Path_Cloud = Cwp_Dcomp(Elem_Idx,Line_Idx) / 1000.0   !kg/m^2
      T_Cloud = ACHA%Tc(Elem_Idx,Line_Idx)
      Reff_Cloud = Reff_Dcomp(Elem_Idx,Line_Idx)
      Tau_Cloud = Tau_Dcomp(Elem_Idx,Line_Idx)

     !--- compute Number concentration and Geometrical Height
     !--- filter for the clouds where this is applicable
     if (T_cloud > 268.0 .and. T_cloud < 300.0 .and. &
         Reff_Cloud > 5.0 .and. Reff_Cloud < 35.0 .and. &
         Tau_Cloud > 5 .and. Tau_Cloud < 50.0) then

         !-- condensation rate (kg/m^3/m)
         Condensation_Rate = exp(-21.0553+T_cloud*0.0536887) / 1000.0

         !geometrical height (meters)
         Hcld_Dcomp(Elem_Idx,Line_Idx)  =  &
                    ( 2.0 / Condensation_Rate * Water_Path_Cloud )**0.5

         !Number concentration (cm-3)
         Cdnc_Dcomp(Elem_Idx,Line_Idx) = 2.0**(-5.0/2.0)/Drop_Dis_Width *   &
                             Tau_Cloud**3.0 * Water_Path_Cloud**(-5.0/2.0) *  &
                            (3.0/5.0*Q_eff_sca*pi)**(-3.0) *  &
                            (3.0*Condensation_Rate/4.0/pi/rho_water)**(-2.0) * &
                             Condensation_Rate**(5.0/2.0)/ 1.0E6
     endif

    end do element_loop
  end do line_loop

end subroutine COMPUTE_ADIABATIC_CLOUD_PROPS

!---------------------------------------------------------------------------------------
! compute precipitation from the KNMI approach
!
! Citation: Roebeling, R. A., and I. Holleman (2009), 
!           SEVIRI rainfall retrieval and validation using weather radar observations, 
!           J. Geophys. Res., 114, D21202, doi:10.1029/2009JD012102.
!
! /***************** PROCEDURE Calculate_Precip
!! *********************************
! * Procedure    : Calculate_Precip
! * Description  : Retrieves Precip
!
! * Author       : Rob Roebeling
! * Calls                :
!
!******************************************************************************/
!void Calculate_PRECIP( float *precip, float tau, float reff, int phase, int cch
!)
!{
!  /* - Declarations */
!  float precip_temp, precip_max, lwp_p,dprecip;
!
!  /* Start Sub-Routine */
!
!  /* Parameterization based on Akos Horvarh and Roger Davies, 2007 */
!  lwp_p = tau * reff * 2.0/3.0;
!
!  if (lwp_p >= 150.0 && ((reff >= 15.0 && phase == J_LIQ) || phase == J_ICE) )
!  {
!
!    dprecip    = (float)cch;
!    if ((float)cch > 7000.0)  dprecip  = 7000.0;
!
!    precip_temp = pow((lwp_p/120.0-1.0),1.6)/(dprecip/1000.0);
!    precip_max  = 5.0+pow((float)cch/1000.0,1.6);
!
!    if (precip_temp >= precip_max) precip_temp = precip_max;
!
!    *precip     = precip_temp;
!  }
!  else
!  {
!      *precip  = J_NCL;
!  }
!}
!/* ------------------------------------------------------------------------- */
!---------------------------------------------------------------------------------------
subroutine COMPUTE_PRECIPITATION(Line_Idx_Min,Num_Lines)

  integer, intent(in):: Line_Idx_Min
  integer, intent(in):: Num_Lines

  real (kind=real4), parameter:: dH_0 = 0.6 !km
  real (kind=real4), parameter:: CWP_0 = 120.0 !g/m^2
  real (kind=real4), parameter:: Lapse_Rate = 6.5 !K/km
  real (kind=real4), parameter:: Alpha = 1.6 !dimensionless
  real (kind=real4), parameter:: C = 1.0 !mm / hour
  real (kind=real4), parameter:: CWP_T = 150.0 !g/m^2
  real (kind=real4), parameter:: Ceps_T = 15.0 !micron
  integer, parameter:: N_box = 100
  real (kind=real4), parameter:: dH_Max = 7.0  !km
  real (kind=real4) :: CTT_Max
  real (kind=real4) :: CTT_Pix
  real (kind=real4) :: Reff_Pix
  real (kind=real4) :: CWP_Pix
  real (kind=real4) :: dH
  real (kind=real4) :: Rain_Rate_Max
  integer:: Line_Idx_Max
  integer:: Elem_Idx_Max
  integer:: Elem_Idx_Min
  integer:: Num_Elements
  integer:: Line_Idx
  integer:: Line_Idx_1
  integer:: Line_Idx_2
  integer:: Elem_Idx
  integer:: Elem_Idx_1
  integer:: Elem_Idx_2


  Elem_Idx_Min = 1
  Num_Elements = Image%Number_Of_Elements
  Elem_Idx_Max = Elem_Idx_Min + Num_Elements - 1
  Line_Idx_Max = Line_Idx_Min + Num_Lines - 1

  Rain_Rate_Dcomp = Missing_Value_Real4

  line_loop: DO Line_Idx = Line_Idx_Min, Line_Idx_Max

    Line_Idx_1  = min(Line_Idx_Max-1,max(1,Line_Idx - N_box /2))
    Line_Idx_2  = min(Line_Idx_Max,max(2,Line_Idx + N_box /2))

    element_loop: DO Elem_Idx = Elem_Idx_Min, Elem_Idx_Max

      !--- skip bad pixels
      if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) cycle

      !--- skip non cloud pixels
      if (Cld_Type(Elem_Idx,Line_Idx) == sym%CLEAR_TYPE .or. &
          Cld_Type(Elem_Idx,Line_Idx) == sym%PROB_CLEAR_TYPE) then
          Rain_Rate_Dcomp(Elem_Idx,Line_Idx) = 0.0
          cycle
      endif

      CWP_Pix = Cwp_Dcomp(Elem_Idx,Line_Idx)
      CTT_Pix = ACHA%Tc(Elem_Idx,Line_Idx)
      Reff_Pix = Reff_Dcomp(Elem_Idx,Line_Idx)

      if (Geo%Solzen(Elem_Idx,Line_Idx) > 90.0) then
        Reff_Pix = Reff_Nlcomp(Elem_Idx,Line_Idx)
      endif

      !--- skip bad pixels
      if (Reff_Pix <= 0.0) cycle

      !--- screen low water path clouds
      if (CWP_Pix < CWP_T) then
        Rain_Rate_Dcomp(Elem_Idx,Line_Idx) = 0.0
        cycle
      endif

      !--- screen small particle size liquid water clouds
      if (Cld_Type(Elem_Idx,Line_Idx) == sym%FOG_TYPE .or. &
          Cld_Type(Elem_Idx,Line_Idx) == sym%WATER_TYPE .or. &
          Cld_Type(Elem_Idx,Line_Idx) == sym%SUPERCOOLED_TYPE) then
         if (Reff_Pix < Ceps_T) then
            Rain_Rate_Dcomp(Elem_Idx,Line_Idx) = 0.0
            cycle
         endif
      endif

     !--- define box to look for maximum cloud temperature
      Line_Idx_1  = min(Line_Idx_Max-1,max(1,Line_Idx - N_box /2))
      Line_Idx_2  = min(Line_Idx_Max,max(2,Line_Idx + N_box /2))
      Elem_Idx_1  = min(Elem_Idx_Max-1,max(1,Elem_Idx - N_box /2))
      Elem_Idx_2  = min(Elem_Idx_Max,max(2,Elem_Idx + N_box /2))
      CTT_Max = maxval(ACHA%Tc(Elem_Idx_1:Elem_Idx_2,Line_Idx_1:Line_Idx_2))  


      !--- compute precip height
      dH = min(dH_Max,(CTT_Max - CTT_Pix) / Lapse_Rate + dH_0)

      !--- compute rain rate
      Rain_Rate_Max = 5.0 + dH**Alpha
      Rain_Rate_Dcomp(Elem_Idx,Line_Idx) = (((CWP_Pix - CWP_0)/CWP_0)**Alpha) / dH
      Rain_Rate_Dcomp(Elem_Idx,Line_Idx) = min(Rain_Rate_Max,Rain_Rate_Dcomp(Elem_Idx,Line_Idx))


    enddo element_loop
  enddo line_loop

end subroutine COMPUTE_PRECIPITATION
!---------------------------------------------------------------------------------------
! compute insolation using the cloud properties
!
! Citation: J.A. Coakley (2003),  Reflectance and Albedo, Surface
!
! Progam takes the cloud tranmission and the spherical albedo from DCOMP to estimate
! solar insolation. 
! 
! Note, the cloud tranmission from DCOMP is a total transmission (diffuse +
! direct).  This routines separates them.
!
! Regression for broad-band solar transmission are taken from (ref here)
!
! Currently, 0.65 um MODIS white sky albedoes are used for the surface albedo
! over land. 
!
! Weaknesses
! 1) We need to develop appropriate direct and diffuse broad-band values of surface albedo
! 2) We need to add aerosol impacts
!
! 
!---------------------------------------------------------------------------------------
subroutine COMPUTE_DCOMP_INSOLATION(Line_Idx_Min,Num_Lines,Sun_Earth_Distance)

  integer, intent(in):: Line_Idx_Min
  integer, intent(in):: Num_Lines
  real, intent(in):: Sun_Earth_Distance

  real (kind=real4), parameter:: SOLAR_CONSTANT = 1356.0 !W/m^2
  real (kind=real4) :: Cloud_Spherical_Albedo
  real (kind=real4) :: Cloud_Optical_Depth
  real (kind=real4) :: Solar_Zenith_Angle
  real (kind=real4) :: Cosine_Solar_Zenith_Angle
  real (kind=real4) :: Cloud_Transmission_Diffuse
  real (kind=real4) :: Cloud_Transmission_Direct
  real (kind=real4) :: Surface_Albedo_Direct
  real (kind=real4) :: Surface_Albedo_Diffuse
  real (kind=real4) :: Insolation_Dcomp_Diffuse_Black_Surface
  real (kind=real4) :: Insolation_Dcomp_Direct_Black_Surface
  real (kind=real4) :: Insolation_Dcomp_Diffuse
  real (kind=real4) :: Insolation_Dcomp_Direct
  real (kind=real4) :: Fo_Toa
  real (kind=real4) :: Fo
  real (kind=real4) :: Tpw
  real (kind=real4) :: Tozone
  integer:: Line_Idx_Max
  integer:: Elem_Idx_Max
  integer:: Elem_Idx_Min
  integer:: Num_Elements
  integer:: Line_Idx
  integer:: Elem_Idx
  integer:: Land_Class
  integer:: Lat_Nwp_Idx
  integer:: Lon_Nwp_Idx
  real:: tau_h2o
  real:: tau_o3 
  real:: tau_co2
  real:: tau_ray 
  real:: tau_total
  real:: atm_trans 
  real:: Surface_Pressure


  Elem_Idx_Min = 1
  Num_Elements = Image%Number_Of_Elements
  Elem_Idx_Max = Elem_Idx_Min + Num_Elements - 1
  Line_Idx_Max = Line_Idx_Min + Num_Lines - 1

  Insolation_Dcomp = Missing_Value_Real4

  Fo_Toa = SOLAR_CONSTANT / (Sun_Earth_Distance**2)

  line_loop: DO Line_Idx = Line_Idx_Min, Line_Idx_Max
    element_loop: DO Elem_Idx = Elem_Idx_Min, Elem_Idx_Max

      Lon_Nwp_Idx = I_Nwp(Elem_Idx,Line_Idx)
      Lat_Nwp_Idx = J_Nwp(Elem_Idx,Line_Idx)

      Cloud_Optical_Depth = Tau_Dcomp(Elem_Idx,Line_Idx)  
      Solar_Zenith_Angle = Geo%Solzen(Elem_Idx,Line_Idx)
      Land_Class = Sfc%Land(Elem_Idx,Line_Idx)
      TPW = Tpw_Nwp_Pix(Elem_Idx,Line_Idx)
      Tozone = 0.0
      Surface_Pressure = 1010.00
      if (Lon_Nwp_Idx > 0 .and. Lat_Nwp_Idx > 0) then
          Tozone = 0.001*Ozone_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx)   !atm-cm
          Surface_Pressure = Psfc_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx) !hPa
      endif

      !--- adjust gases for slant path
      Cosine_Solar_Zenith_Angle = cos(Geo%Solzen(Elem_Idx,Line_Idx)*DTOR)

!     H2O_Trans_Direct = 1.0 - 2.9*TPW / ((1.0 + 141.5*TPW)**(0.635) + 5.925*TPW)
!     Ozone_Trans_Direct = 1.0 - 0.02118*Ozone / (1.0 + 0.042*Ozone + 0.000323*(Ozone**2))
!     Fo = Fo_Toa * H2O_Trans_Direct * Ozone_Trans_Direct


      tau_h2o = 0.104*(Tpw**0.30)
      tau_o3 = 0.038*(Tozone**0.44)
      tau_co2 = 0.0076*((Surface_Pressure/1013.25)**0.29)
      tau_ray = 0.038*(Surface_Pressure/1013.25)
      tau_total = tau_ray + tau_h2o + tau_co2 + tau_ray
      atm_trans = exp(-1.0*tau_total / Cosine_Solar_Zenith_Angle)

      Fo = Fo_Toa * atm_trans

      !--- skip data that can not be processed
      if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) cycle
      if (Solar_Zenith_Angle > 70.0) cycle

      !--- determine surface albedo
      if (Land_Class == sym%LAND) then
         Surface_Albedo_Diffuse = ch(1)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx) / 100.0
         Surface_Albedo_Direct = Surface_Albedo_Diffuse
      else
         Surface_Albedo_Diffuse = 0.06
         Surface_Albedo_Direct = 0.026/(Cosine_Solar_Zenith_Angle**1.7+0.065) + &
                                 0.15*(Cosine_Solar_Zenith_Angle - 0.1)* &
                                      (Cosine_Solar_Zenith_Angle-0.5)*   &
                                      (Cosine_Solar_Zenith_Angle-1.0)
      endif

      !-- set cloud trans and albedo, if clear, set to transparent values
      if (Cloud_Optical_Depth > 0.0) then
        Cloud_Spherical_Albedo = Cloud_063um_Spherical_Albedo(Elem_Idx,Line_Idx)
        Cloud_Transmission_Direct = exp( -1.0 * Cloud_Optical_Depth / Cosine_Solar_Zenith_Angle)
        Cloud_Transmission_Diffuse = Cloud_063um_Transmission_Solar(Elem_Idx,Line_Idx)  - Cloud_Transmission_Direct
      else
        Cloud_Spherical_Albedo = 0.0
        Cloud_Transmission_Direct =  1.0
        Cloud_Transmission_Diffuse =  0.0
      endif
       
      Insolation_Dcomp_Direct_Black_Surface = Fo * Cloud_Transmission_Direct * Cosine_Solar_Zenith_Angle

      Insolation_Dcomp_Diffuse_Black_Surface = Fo * Cloud_Transmission_Diffuse * Cosine_Solar_Zenith_Angle

!     !-- Coakley
!     Insolation_Dcomp_Direct = Insolation_Dcomp_Direct_Black_Surface * (1.0 +  &
!                 Surface_Albedo_Direct * Cloud_Spherical_Albedo / (1.0 - Surface_Albedo_Diffuse * Cloud_Spherical_Albedo))

!     Insolation_Dcomp_Diffuse = Insolation_Dcomp_Diffuse_Black_Surface / &
!                                (1.0 - Surface_Albedo_Diffuse * Cloud_Spherical_Albedo) 

      !-- Heidinger Formulation
      Insolation_Dcomp_Direct = Insolation_Dcomp_Direct_Black_Surface * (1.0 +  &
                  Surface_Albedo_Direct * Cloud_Spherical_Albedo / (1.0 - Surface_Albedo_Diffuse * Cloud_Spherical_Albedo))

      Insolation_Dcomp_Diffuse = Insolation_Dcomp_Diffuse_Black_Surface *  &
                  (1.0 + Surface_Albedo_Diffuse * Cloud_Spherical_Albedo / (1.0 - Surface_Albedo_Diffuse * Cloud_Spherical_Albedo))

      !--- combine
      Insolation_Dcomp(Elem_Idx,Line_Idx) = Insolation_Dcomp_Direct + Insolation_Dcomp_Diffuse

    enddo element_loop
  enddo line_loop

end subroutine COMPUTE_DCOMP_INSOLATION

!-----------------------------------------------------------
! end of MODULE
!-----------------------------------------------------------
end module DCOMP_DERIVED_PRODUCTS_MODULE
