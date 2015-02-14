!$Id$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: nwp_common.f90 (src)
!       NWP_COMMON (program)
!
! PURPOSE:  This module holds the radiative transfer quantities needed for
!           the algorithms
!
! DESCRIPTION: 
!           note, there two type of nwp data
!            1- the pressure level data
!            2- the data on different surface grid
!
!            the only data assumed to be a on the surface grid are
!              - surface temperature
!              - Weasd depth
!              - u and v wind speed at 10m
!
!           the surface and pressure level grid may be different
!           i_Nwp, j_Nwp points to a cell in the pressure level data
!
!           In the GFS data, the pressure and surface grids are the same, in the
!           NCEP reanalysis, they differ.
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
! public::
!   CREATE_NWP_ARRAYS
!   DESTROY_NWP_ARRAYS
!   FIND_NWP_GRID_CELL
!   MAP_PIXEL_NWP
!   KNOWING_P_COMPUTE_T_Z_NWP
!   KNOWING_Z_COMPUTE_T_P_NWP
!   KNOWING_T_COMPUTE_P_Z_NWP
!   FIND_NWP_LEVELS
!   INTERPOLATE_NWP
!   INTERPOLATE_PROFILE
!   INTERPOLATE_NWP_TZ_PROFILES
!   COMPUTE_COAST_MASK_NWP
!   QC_NWP
!   COMPUTE_PIXEL_NWP_PARAMETERS
!   MODIFY_TSFC_NWP_PIX
!   COMPUTE_NWP_PARAMETERS
!   PROF_LOOKUP_USING_P
!   PROF_LOOKUP_USING_T
!   PROF_LOOKUP_USING_T_LAPSE
!   PROF_LOOKUP_USING_T_PROF
!   PROF_LOOKUP_USING_Z
!   TEMPORAL_INTERP_TMPSFC_NWP
!   CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY
!--------------------------------------------------------------------------------------
module NWP_COMMON
   
  use CONSTANTS
  use PIXEL_COMMON
  use NUMERICAL_ROUTINES

  implicit none
  private
  private:: FIND_NWP_LEVELS, &
            COMPUTE_NWP_CLOUD_PARAMETERS, &
            CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY

  public:: CREATE_NWP_ARRAYS,  &
           DESTROY_NWP_ARRAYS,  &
           FIND_NWP_GRID_CELL,  &
           MAP_PIXEL_NWP, &
           KNOWING_P_COMPUTE_T_Z_NWP,  &
           KNOWING_Z_COMPUTE_T_P_NWP,  &
           KNOWING_Z_COMPUTE_T_P_NWP_ARBITRARY_LEVELS,  &
           KNOWING_T_COMPUTE_P_Z_NWP,  &
           COMPUTE_NWP_LEVELS_SEGMENT, &
           INTERPOLATE_NWP,  &
           INTERPOLATE_PROFILE,  &
           INTERPOLATE_NWP_TZ_PROFILES,  &
           COMPUTE_COAST_MASK_NWP, &
           QC_NWP, &
           MODIFY_TSFC_NWP_PIX, &
           COMPUTE_SEGMENT_NWP_CLOUD_PARAMETERS, &
           PROF_LOOKUP_USING_P, &
           PROF_LOOKUP_USING_T, &
           PROF_LOOKUP_USING_T_LAPSE, &
           PROF_LOOKUP_USING_T_PROF, &
           PROF_LOOKUP_USING_Z, &
           TEMPORAL_INTERP_TMPSFC_NWP, &
           COMPUTE_PIXEL_NWP_PARAMETERS

 interface CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY
     module procedure  &
         CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY_I1, &
         CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY_I2, &
         CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY_I4, &
         CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY_R4
 end interface

 interface INTERPOLATE_NWP
     module procedure  &
         INTERPOLATE_NWP_I1, &
         INTERPOLATE_NWP_I2, &
         INTERPOLATE_NWP_I4, &
         INTERPOLATE_NWP_R4
 end interface
!----------------------------------------------------------------------
!--- set this parameter to 1 when reading GFS hdf files that have
!--- x as the first index, not z
!----------------------------------------------------------------------
  integer, public, save:: REFORMAT_GFS_ZXY

! NWP array declarations
  integer (kind=int4), save, public :: npoints, Nlevels_Nwp
  integer (kind=int4), save, public :: Nlat_Nwp, Nlon_Nwp
  real (kind=real4),   save, public :: dLat_Nwp, dLon_Nwp
  real (kind=real4), public :: lat1_Nwp
  real (kind=real4), public :: lon1_Nwp
  real (kind=real4), save, public :: missing_Nwp
  real, public, parameter :: Psfc_max_Nwp = 1100.0
  real (kind=real8), public, save :: ncep_time_Before
  real (kind=real8), public, save :: ncep_time_After
  integer,           dimension(:,:), allocatable, public, save :: Mask_Nwp
  integer,           dimension(:,:), allocatable, public, save :: bad_Nwp_mask
  integer,           dimension(:,:), allocatable, public, save :: sfc_type_Nwp
  real (kind=real4), dimension(:,:), allocatable, public, save :: Satzen_Nwp
  real (kind=real4), dimension(:,:), allocatable, public, save :: Solzen_Nwp

  real (kind=real4), dimension(:,:), allocatable, public, save, target :: Psfc_Nwp
  real (kind=real4), dimension(:,:), allocatable, public, save :: Pmsl_Nwp
  real (kind=real4), dimension(:,:), allocatable, public, save :: Zsfc_Nwp
  real (kind=real4), dimension(:,:), allocatable, public, save :: Tmpsfc_Nwp
  real (kind=real4), dimension(:,:), allocatable, public, save :: Tmpsfc_Nwp_Before
  real (kind=real4), dimension(:,:), allocatable, public, save :: Tmpsfc_Nwp_After
  real (kind=real4), dimension(:,:), allocatable, public, save :: Tmpair_Nwp
  real (kind=real4), dimension(:,:), allocatable, public, save :: Tmpair_uni_Nwp
  real (kind=real4), dimension(:,:), allocatable, public, save :: Tmpsfc_uni_Nwp
  real (kind=real4), dimension(:,:), allocatable, public, save :: T_Trop_Nwp
  real (kind=real4), dimension(:,:), allocatable, public, save :: P_Trop_Nwp
  real (kind=real4), dimension(:,:), allocatable, public, save :: Rhsfc_Nwp
  real (kind=real4), dimension(:,:), allocatable, public, save :: Tpw_Nwp
  real (kind=real4), dimension(:,:), allocatable, public, save :: Uth_Nwp
  real (kind=real4), dimension(:,:), allocatable, public, save :: Hght500_Nwp
  real (kind=real4), dimension(:,:), allocatable, public, save :: Ozone_Nwp
  real (kind=real4), dimension(:,:), allocatable, public, save :: Weasd_Nwp
  real (kind=real4), dimension(:,:), allocatable, public, save :: U_Wnd_10m_Nwp
  real (kind=real4), dimension(:,:), allocatable, public, save :: V_Wnd_10m_Nwp
  real (kind=real4), dimension(:,:), allocatable, public, save :: Wnd_Spd_10m_Nwp
  real (kind=real4), dimension(:,:), allocatable, public, save :: Wnd_Dir_10m_Nwp
  integer (kind=int1), dimension(:,:), allocatable, public, save :: Land_Nwp
  real (kind=real4), dimension(:,:), allocatable, public, save :: Ice_Nwp
  integer (kind=int1), dimension(:,:), allocatable, public, save :: Tropo_Level_Nwp
  integer (kind=int4), dimension(:,:), allocatable, public, save :: Level850_Nwp
  integer (kind=int4), dimension(:,:), allocatable, public, save :: Level700_Nwp
  integer (kind=int4), dimension(:,:), allocatable, public, save :: Level500_Nwp
  integer (kind=int1), dimension(:,:), allocatable, public, save :: Sfc_Level_Nwp
  integer (kind=int1), dimension(:,:), allocatable, public, save :: Inversion_Level_Nwp
  integer (kind=int4), dimension(:,:,:), allocatable, public, save :: Inversion_Level_Profile_Nwp
  real (kind=real4), dimension(:,:), allocatable, public, save :: Lifting_Condensation_Level_Height_Nwp !km
  real (kind=real4), dimension(:,:), allocatable, public, save :: Convective_Condensation_Level_Height_Nwp !km
  real (kind=real4), dimension(:,:), allocatable, public, save :: Freezing_Level_Height_Nwp !km
  real (kind=real4), dimension(:,:), allocatable, public, save :: Upper_Limit_Water_Height_Nwp !km
  real (kind=real4), dimension(:,:), allocatable, public, save :: K_Index_Nwp !K
  real (kind=real4), dimension(:,:), allocatable, public, save :: Pc_Nwp
  real (kind=real4), dimension(:,:), allocatable, public, save :: Sc_Lwp_Nwp
  real (kind=real4), dimension(:,:), allocatable, public, save :: Lwp_Nwp
  real (kind=real4), dimension(:,:), allocatable, public, save :: Iwp_Nwp
  real (kind=real4), dimension(:,:), allocatable, public, save :: Cwp_Nwp
  real (kind=real4), dimension(:,:), allocatable, public, save :: Cloud_Fraction_Satellite_Nwp
  real (kind=real4), dimension(:,:), allocatable, public, save :: High_Cloud_Fraction_Satellite_Nwp
  real (kind=real4), dimension(:,:), allocatable, public, save :: Mid_Cloud_Fraction_Satellite_Nwp
  real (kind=real4), dimension(:,:), allocatable, public, save :: Low_Cloud_Fraction_Satellite_Nwp
  integer (kind=int1), dimension(:,:), allocatable, public, save :: Ncld_Layers_Nwp
  integer (kind=int1), dimension(:,:), allocatable, public, save :: Cld_Type_Nwp
  real (kind=real4), dimension(:), allocatable, public, save :: Lat_Nwp
  real (kind=real4), dimension(:), allocatable, public, save :: Lon_Nwp
  real (kind=real4), dimension(:), allocatable, target, public, save :: P_Std_Nwp
  real (kind=real4), dimension(:,:,:), allocatable, target, save, public :: Z_Prof_Nwp
  real (kind=real4), dimension(:,:,:), allocatable, target, save, public :: T_Prof_Nwp
  real (kind=real4), dimension(:,:,:), allocatable, save, public :: Rh_Prof_Nwp
  real (kind=real4), dimension(:,:,:), allocatable, save, public :: Ozone_Prof_Nwp
  real (kind=real4), dimension(:,:,:), allocatable, target, save, public :: Tpw_Prof_Nwp
  real (kind=real4), dimension(:,:,:), allocatable, target, save, public :: Clwmr_Prof_Nwp
  real (kind=real4), dimension(:,:,:), allocatable, target, save, public :: U_Wnd_Prof_Nwp
  real (kind=real4), dimension(:,:,:), allocatable, target, save, public :: V_Wnd_Prof_Nwp
  real (kind=real4), dimension(:,:), allocatable, save, public :: temp2d_Nwp_1
  real (kind=real4), dimension(:,:), allocatable, save, public :: temp2d_Nwp_2
  real (kind=real4), dimension(:,:,:), allocatable, save, public :: temp3d_Nwp_1
  real (kind=real4), dimension(:,:,:), allocatable, save, public :: temp3d_Nwp_2
  real (kind=real4), dimension(:,:,:), allocatable, save, public :: temp3d

  integer(kind=int4), save, public:: nwp_start_hour
  integer(kind=int4), save, public:: nwp_end_hour

!--- nwp profiles interpolated to the pixel level
  real (kind=real4), dimension(:), allocatable, public, save :: T_Prof_Nwp_pix
  real (kind=real4), dimension(:), allocatable, public, save :: Z_Prof_Nwp_pix

!--- local parameters
  real(kind=real4), public, parameter :: P_Trop_Max = 300.0
  real(kind=real4), public, parameter :: P_Trop_Min = 25.0
  real(kind=real4), public, parameter :: P_Inversion_Min = 700.0
  real(kind=real4), public, parameter :: Delta_T_Inversion = 0.0  !05

contains
!-------------------------------------------------------------
! subroutine QC_NWP()
!
! Subroutine to quality control NWP data
!
! Check the values of some fields and set Bad_Nwp_Mask
! accordingly
!
! The tests run here are arbitrary but are based on known
! failures
!
!-------------------------------------------------------------
 subroutine QC_NWP()

  integer:: Lon_Nwp_Idx, Lat_Nwp_Idx

  
  do Lon_Nwp_Idx = 1, Nlon_Nwp
     do Lat_Nwp_Idx = 1, Nlat_Nwp

        if ((P_Trop_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx) <= 0.0) .or. &
            (T_Trop_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx) <= 0.0) .or. &
            (Zsfc_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx) > 10000.0) .or. &
            (Psfc_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx) > 1500.0) .or. &
            (Tmpsfc_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx) > 400.0) .or. &
            (Tmpsfc_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx) <= 0.0)) then

            Bad_Nwp_Mask(Lon_Nwp_Idx,Lat_Nwp_Idx) = sym%YES

        endif

     enddo
  enddo

end subroutine QC_NWP

!----------------------------------------------------------------------
! subroutine COMPUTE_TSFC_NWP(i1,nx,j1,ny,Smooth_nwp_opt)
!
! compute a pixel level surface temperature from the NWP fields
! and smooth if option chosen
!
! i1 = first element index
! nx = number of element indices
! j1 = first element index
! ny = number of element indices
! Smooth_nwp_opt = flag to smooth nwp
!
! This must be called After MAP_PIXEL_NWP
!----------------------------------------------------------------------
subroutine MODIFY_TSFC_NWP_PIX(Elem_Idx_Start,Num_Elements,Line_Idx_Start,Num_Lines)

  integer(kind=int4), intent(in):: Elem_Idx_Start
  integer(kind=int4), intent(in):: Num_Elements
  integer(kind=int4), intent(in):: Line_Idx_Start
  integer(kind=int4), intent(in):: Num_Lines
  integer(kind=int4) :: Elem_Idx_End
  integer(kind=int4) :: Line_Idx_End
  integer(kind=int4) :: Elem_Idx
  integer(kind=int4) :: Line_Idx
  real (kind=real4) :: Delta_Zsfc
  real (kind=real4) :: Delta_Tsfc
  real (kind=real4) :: Delta_Lapse_Rate
  real(kind=real4) :: Zsfc_Nwp_Pix
  integer(kind=int4) :: Ilev_start
  integer(kind=int4) :: Ilev_end
  integer(kind=int4) :: Lon_Nwp_Idx
  integer(kind=int4) :: Lat_Nwp_Idx
  integer(kind=int4) :: Lon_Nwp_Idx_x
  integer(kind=int4) :: Lat_Nwp_Idx_x
  integer(kind=int4) :: Sfc_Level_Idx

  Elem_Idx_End = Elem_Idx_Start + Num_Elements - 1
  Line_Idx_End = Line_Idx_Start + Num_Lines - 1

  do Elem_Idx = Elem_Idx_Start, Elem_Idx_End
    do Line_Idx = Line_Idx_Start,Line_Idx_End

     if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) then
             Tsfc_Nwp_Pix(Elem_Idx,Line_Idx) = Missing_Value_Real4
             cycle
     endif

     Lon_Nwp_Idx = i_Nwp(Elem_Idx,Line_Idx)
     Lat_Nwp_Idx = j_Nwp(Elem_Idx,Line_Idx)

     if (Lon_Nwp_Idx == 0 .or. Lat_Nwp_Idx == 0) then
             Tsfc_Nwp_Pix(Elem_Idx,Line_Idx) = Missing_Value_Real4
             cycle
     endif

     Lon_Nwp_Idx_x = i_Nwp_x(Elem_Idx,Line_Idx)
     Lat_Nwp_Idx_x = j_Nwp_x(Elem_Idx,Line_Idx)
     Sfc_Level_Idx = Sfc_Level_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx)

     !----------------------------------------------------------------------------------
     ! modify Tsfc_Nwp_Pix for sub-nwp elevation
     !
     !  Zsfc = pixel level elevation in meters
     !  Zsfc_Nwp = nwp level elevation in km
     !
     !----------------------------------------------------------------------------------
     if (Sfc%Land(Elem_Idx,Line_Idx) == sym%LAND) then

        !--- assume all surface features are in lowest half of profile
        Ilev_end = Nlevels_Nwp
        Ilev_start = Nlevels_Nwp/2

        if ((Zsfc_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx) /= Missing_Value_Real4) .and. &
           (Sfc%Zsfc(Elem_Idx,Line_Idx) /= Missing_Value_Real4) .and. &
           (Lon_Nwp_Idx > 0) .and. (Lat_Nwp_Idx > 0) .and. (Lon_Nwp_Idx_x > 0) .and. (Lat_Nwp_Idx_x > 0) .and. &
           (Sfc_Level_Idx > 1)) then

          !--- compute a smooth surface elevation from NWP 
          Zsfc_Nwp_Pix = INTERPOLATE_NWP( &
                       Zsfc_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx), &
                       Zsfc_Nwp(Lon_Nwp_Idx_x,Lat_Nwp_Idx), &
                       Zsfc_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx_x), &
                       Zsfc_Nwp(Lon_Nwp_Idx_x,Lat_Nwp_Idx_x), &
                       Lon_Nwp_fac(Elem_Idx,Line_Idx), &
                       Lat_Nwp_fac(Elem_Idx,Line_Idx))
        
          !--- compute the near surface lapse rate (K/m) 
          Delta_Lapse_Rate = (T_Prof_Nwp(Sfc_Level_Idx,Lon_Nwp_Idx,Lat_Nwp_Idx) &
                              - T_Prof_Nwp(Sfc_Level_Idx-1,Lon_Nwp_Idx,Lat_Nwp_Idx)) / &
                            (Z_Prof_Nwp(Sfc_Level_Idx,Lon_Nwp_Idx,Lat_Nwp_Idx) &
                            - Z_Prof_Nwp(Sfc_Level_Idx-1,Lon_Nwp_Idx,Lat_Nwp_Idx))
        else
          Delta_Lapse_Rate = 0
        endif

        !--- compute the pertubation to NWP surface temp to account for sub-grid elevation
        Delta_Zsfc = Sfc%Zsfc(Elem_Idx,Line_Idx) - Zsfc_Nwp_Pix !meters
        Delta_Tsfc = Delta_Lapse_Rate * Delta_Zsfc       !K
        Tsfc_Nwp_Pix(Elem_Idx,Line_Idx) = Tsfc_Nwp_Pix(Elem_Idx,Line_Idx) + Delta_Tsfc   !K

     endif

    enddo
  enddo   

end subroutine MODIFY_TSFC_NWP_PIX

!----------------------------------------------------------------------
! subroutine COMPUTE_NWP_PARAMETERS(Smooth_nwp_opt)
!
! compute parameters from NWP fields and smooth if option chosen
!
! This must be called After MAP_PIXEL_NWP
!----------------------------------------------------------------------
subroutine COMPUTE_PIXEL_NWP_PARAMETERS(Smooth_nwp_opt)

  integer(kind=int4), intent(in):: Smooth_nwp_opt
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(Tmpsfc_Nwp,Tsfc_Nwp_Pix,Smooth_nwp_opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(T_Trop_Nwp,Ttropo_Nwp_Pix,Smooth_nwp_opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(Tmpair_Nwp,Tair_Nwp_Pix,Smooth_nwp_opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(Rhsfc_Nwp,Rh_Nwp_Pix,Smooth_nwp_opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(Psfc_Nwp,Psfc_Nwp_Pix,Smooth_nwp_opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(Pmsl_Nwp,Pmsl_Nwp_Pix,Smooth_nwp_opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(Tpw_Nwp,Tpw_Nwp_Pix,Smooth_nwp_opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(Ozone_Nwp,Ozone_Nwp_Pix,Smooth_nwp_opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(K_Index_Nwp,K_Index_Nwp_Pix,Smooth_nwp_opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(Pmsl_Nwp,Pmsl_Nwp_Pix,Smooth_nwp_opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(Sc_Lwp_Nwp,Sc_Lwp_Nwp_Pix,Smooth_nwp_opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(Lwp_Nwp,Lwp_Nwp_Pix,Smooth_nwp_opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(Iwp_Nwp,Iwp_Nwp_Pix,Smooth_nwp_opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(Cwp_Nwp,Cwp_Nwp_Pix,Smooth_nwp_opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(Pc_Nwp,Pc_Nwp_Pix,Smooth_nwp_opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(Cloud_Fraction_Satellite_Nwp,Cfrac_Nwp_Pix,Smooth_nwp_opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(Ncld_Layers_Nwp,Ncld_Layers_Nwp_Pix,Smooth_nwp_opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(Cld_Type_Nwp,Cld_Type_Nwp_Pix,Smooth_nwp_opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(Wnd_Spd_10m_Nwp,Wnd_Spd_10m_Nwp_Pix,Smooth_nwp_opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(Wnd_Dir_10m_Nwp,Wnd_Dir_10m_Nwp_Pix,Smooth_nwp_opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(Lifting_Condensation_Level_Height_Nwp,LCL_Height_Nwp_Pix,Smooth_nwp_opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(Convective_Condensation_Level_Height_Nwp,CCL_Height_Nwp_Pix,Smooth_nwp_opt)

end subroutine COMPUTE_PIXEL_NWP_PARAMETERS

!-------------------------------------------------------------
! subroutine MAP_PIXEL_NWP(j1,j2)
!
! Subroutine to find nwp cell where each in a segment lies
!-------------------------------------------------------------
 subroutine MAP_PIXEL_NWP(Number_of_Elements,Number_of_Lines)

  integer, intent(in):: Number_of_Elements,Number_of_Lines
  integer:: Elem_Idx,Line_Idx,Ierr


  i_Nwp = Missing_Value_int1
  j_Nwp = Missing_Value_int1
  i_Nwp_x = Missing_Value_int1
  j_Nwp_x = Missing_Value_int1
  Lon_Nwp_Fac = Missing_Value_Real4
  Lat_Nwp_Fac = Missing_Value_Real4

  do Line_Idx = 1, Number_of_Lines
     do Elem_Idx = 1, Number_of_Elements
                                                                                                                                         
      !--- check for valid geolocation
      if (Nav%Lon(Elem_Idx,Line_Idx) < -180.0 .or. Nav%Lon(Elem_Idx,Line_Idx) > 180.0 .or. &
          Nav%Lat(Elem_Idx,Line_Idx) < -90.0 .or. Nav%Lat(Elem_Idx,Line_Idx) > 90.0) then
          cycle
      endif 
        
      !--- compute NWP cell to pixel mapping
      call FIND_NWP_GRID_CELL(Nav%Lon(Elem_Idx,Line_Idx),Nav%Lat(Elem_Idx,Line_Idx), &
                              i_Nwp(Elem_Idx,Line_Idx),j_Nwp(Elem_Idx,Line_Idx), &
                              i_Nwp_x(Elem_Idx,Line_Idx),j_Nwp_x(Elem_Idx,Line_Idx),  &
                              Lon_Nwp_fac(Elem_Idx,Line_Idx), Lat_Nwp_fac(Elem_Idx,Line_Idx),Ierr)

       !-- if there is an error, flag pixel as bad
      if (Ierr == 1) then
         i_Nwp(Elem_Idx,Line_Idx) = Missing_Value_int1
         j_Nwp(Elem_Idx,Line_Idx) = Missing_Value_int1
         i_Nwp_x(Elem_Idx,Line_Idx) = Missing_Value_int1
         j_Nwp_x(Elem_Idx,Line_Idx) = Missing_Value_int1
         Lon_Nwp_fac(Elem_Idx,Line_Idx) = Missing_Value_Real4
         Lat_Nwp_fac(Elem_Idx,Line_Idx) = Missing_Value_Real4
         cycle
      endif

     enddo
  enddo

 end subroutine MAP_PIXEL_NWP
!------------------------------------------------------------------
! Compute NWP Levels for each NWP Gridcell
!
! must be called aftrer MAP_PIXEL_NWP
!------------------------------------------------------------------
 subroutine COMPUTE_NWP_LEVELS_SEGMENT(Number_of_Elements,Number_of_Lines)

  integer, intent(in):: Number_of_Elements
  integer, intent(in):: Number_of_Lines
  integer:: Elem_Idx
  integer:: Line_Idx
  integer:: Lat_NWP_Idx
  integer:: Lon_NWP_Idx

  !--- intialize levels to missing
  Sfc_Level_Nwp = Missing_Value_Int1
  Tropo_Level_Nwp = Missing_Value_Int1
  Inversion_Level_Nwp = Missing_Value_Int1
  Level850_Nwp = Missing_Value_Int1
  Level700_Nwp = Missing_Value_Int1
  Level500_Nwp = Missing_Value_Int1

  !--- loop through each pixel and if the
  do Line_Idx = 1, Number_of_Lines
     do Elem_Idx = 1, Number_of_Elements

      !--- alias nwp indices for this nwp cell using predetermined global variables
      Lon_NWP_Idx = i_Nwp(Elem_Idx,Line_Idx)   
      Lat_NWP_Idx = j_Nwp(Elem_Idx,Line_Idx)  

      !--- check for valid nwp mapping and data, if not skip
      if (Lon_NWP_Idx < 1 .or. Lat_Nwp_Idx < 1) cycle
      if (Bad_Nwp_Mask(Lon_NWP_Idx,Lat_NWP_Idx) == sym%YES) cycle

      !-- if this populated for this nwp cell, skip this pixel
      if (Sfc_Level_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx) /= Missing_Value_Int1) cycle

      !--- find needed nwp levels for this nwp cell, store in global variables
      call FIND_NWP_LEVELS(Lon_Nwp_Idx,Lat_Nwp_Idx)

     enddo
  enddo

 end subroutine COMPUTE_NWP_LEVELS_SEGMENT
!------------------------------------------------------------------
! subroutine FIND_NWP_GRID_CELL(lon, lat, Lon_Nwp_Idx, Lat_Nwp_Idx, Lon_Nwp_Idxx, Lat_Nwp_Idxx, lonfac, latfac, ierror)
!
! Subroutine to convert lat, lon into NWP grid-cell
!
! input:
!   lon - longitude (-180 to 180)
!   lat - latitude (-90 to 90)
! output:
!   Lat_Nwp_Idx - nwp latitude index of the nearest nwp latitude
!   Lon_Nwp_Idx - nwp longitude index of the nearest nwp longitude
!   Lat_Nwp_Idxx - nwp latitude index of the nearest nwp latitude diagonal
!   Lon_Nwp_Idxx - nwp longitude index of the nearest nwp longitude diagonal
!   latfac - latitude weight between Lat_Nwp_Idx(0.0) and Lat_Nwp_Idxx(1.0)
!   lonfac - longitude weight between Lon_Nwp_Idx(0.0) and Lon_Nwp_Idxx(1.0)
!  
!
! imagine a pixel, x, surrounded by nwp vertices (O)
!
!        O            O            O
!                     -   x
!
!        O            O            O
!                                 --- 
!
!  point O is (Lon_Nwp_Idx,Lat_Nwp_Idx) and O is (Lon_Nwp_Idxx,Lat_Nwp_Idxx)
!        -                   ---
!
! modified to return information needed to spatially interpolate
!------------------------------------------------------------------
  subroutine FIND_NWP_GRID_CELL(lon, lat, Lon_Nwp_Idx, Lat_Nwp_Idx, Lon_Nwp_Idxx, Lat_Nwp_Idxx, lonfac, latfac, ierror)

    real (kind=real4), intent(in) :: lon, lat
    integer (kind=int4), intent(out) :: Lon_Nwp_Idx, Lat_Nwp_Idx, Lon_Nwp_Idxx, Lat_Nwp_Idxx, ierror
    real (kind=real4), intent(out) :: latfac
    real (kind=real4), intent(out) :: lonfac

    real (kind=real4) :: rlat
    real (kind=real4) :: rlon
    real (kind=real4) :: rLon_Nwp
    real (kind=real4) :: rLon_Nwpx
    integer:: Is_Dateline

    ierror = 0
    rlon = lon
    rlat = lat
    Is_Dateline = 0

    !--- convert negative lons to go from 180 to 360 degrees
    if (rlon < 0.0) then
       rlon = rlon + 360.0
    endif

    !--- Find Position in NWP grid
    if (rlon < 0.0 .or. rlon > 360.0 .or. rlat < -90.0 .or. rlat > 90.0) then
       ierror = 1
       Lat_Nwp_Idx = 0
       Lon_Nwp_Idx = 0
    else
       Lat_Nwp_Idx = max(1, min(Nlat_Nwp, nint( (rlat-lat1_Nwp) / dLat_Nwp + 1.0) ))
       Lon_Nwp_Idx = max(1, min(Nlon_Nwp, nint( (rlon-lon1_Nwp) / dLon_Nwp + 1.0) ))
    endif

    rLon_Nwp = Lon_Nwp(Lon_Nwp_Idx)
    if (Lon_Nwp(Lon_Nwp_Idx) < 0.0) then
       rLon_Nwp = rLon_Nwp + 360.0
    endif

    !---  latitude interpolation information
    if (Lat_Nwp_Idx > 1 .and. Lat_Nwp_Idx < Nlat_Nwp) then
      if (sign(1.0,lat-Lat_Nwp(Lat_Nwp_Idx)) == sign(1.0,dLat_Nwp)) then
         Lat_Nwp_Idxx = Lat_Nwp_Idx + 1
      else
         Lat_Nwp_Idxx = Lat_Nwp_Idx - 1
      endif
      Lat_Nwp_Idxx = min(Nlat_Nwp,max(1,Lat_Nwp_Idxx))

      !--- compute latitude interpolation factor
      if (Lat_Nwp(Lat_Nwp_Idxx) /= Lat_Nwp(Lat_Nwp_Idx)) then
         latfac = (lat - Lat_Nwp(Lat_Nwp_Idx))/(Lat_Nwp(Lat_Nwp_Idxx)-Lat_Nwp(Lat_Nwp_Idx))
      else
         latfac = 0.0
      endif
      
    endif

    !---- determine dateline flag
    if (abs(lon - Lon_Nwp(Lon_Nwp_Idx)) > abs(dLon_Nwp)) then
       Is_Dateline = 1
    endif

    !---  longitude interpolation information
    if (Lon_Nwp_Idx > 1 .and. Lon_Nwp_Idx < Nlon_Nwp) then
      if (sign(1.0,lon-Lon_Nwp(Lon_Nwp_Idx)) == sign(1.0,dLon_Nwp)) then
          if (Is_Dateline == 0) then
            Lon_Nwp_Idxx = Lon_Nwp_Idx + 1
          else
            Lon_Nwp_Idxx = Lon_Nwp_Idx - 1
          endif
      else
          if (Is_Dateline == 0) then
            Lon_Nwp_Idxx = Lon_Nwp_Idx - 1
          else
            Lon_Nwp_Idxx = Lon_Nwp_Idx + 1
          endif
      endif
    endif
    Lon_Nwp_Idxx = min(Nlon_Nwp,max(1,Lon_Nwp_Idxx))

    !--- make a positive definite value of lon at Lon_Nwp_Idxx
    rLon_Nwpx = Lon_Nwp(Lon_Nwp_Idxx)
    if (Lon_Nwp(Lon_Nwp_Idxx) < 0.0) then
       rLon_Nwpx = rLon_Nwpx + 360.0
    endif

    !--- recompute date line flag including Lon_Nwp_Idxx point
    if (abs(lon - Lon_Nwp(Lon_Nwp_Idx)) > abs(dLon_Nwp)) then
       Is_Dateline = 1
    endif
    if (abs(lon - Lon_Nwp(Lon_Nwp_Idxx)) > abs(dLon_Nwp)) then
       Is_Dateline = 1
    endif

    !--- compute latitude interpolation factor
    if (Is_Dateline == 0) then
       if (Lon_Nwp(Lon_Nwp_Idxx)/=Lon_Nwp(Lon_Nwp_Idx)) then
          lonfac = (lon - Lon_Nwp(Lon_Nwp_Idx))/(Lon_Nwp(Lon_Nwp_Idxx)-Lon_Nwp(Lon_Nwp_Idx))
       else
          lonfac = 0.0
       endif
    else
       if (rLon_Nwp /= rLon_Nwpx) then
          lonfac = abs((rlon - rLon_Nwp)/(rLon_Nwp-rLon_Nwpx))
       else
          lonfac = 0.0
       endif
    endif

    !--- constrain
    lonfac = min(0.5,max(0.0,lonfac))
    latfac = min(0.5,max(0.0,latfac))
    Lon_Nwp_Idxx = min(Nlon_Nwp,max(1,Lon_Nwp_Idxx))
    Lat_Nwp_Idxx = min(Nlat_Nwp,max(1,Lat_Nwp_Idxx))

  end subroutine FIND_NWP_GRID_CELL


!----------------------------------------------------------------------------
! Function INTERPOLATE_NWP
! 
! general interpolation routine for nwp fields
!
! description of arguments
! Lon_Nwp_Idx, Lat_Nwp_Idx - nwp indices of closest nwp cell
! Lon_Nwp_Idxx,Lat_Nwp_Idxx - nwp indices of nwp cells of diagnoal of bounding box
! lonx - longitude weighting factor 
! latx = latitude weighting factor
! z1 = data(Lon_Nwp_Idx, Lat_Nwp_Idx)
! z2 = data(Lon_Nwp_Idxx,Lat_Nwp_Idx)
! z3 = data(Lon_Nwp_Idx,Lat_Nwp_Idxx)
! z4 = data(Lon_Nwp_Idxx,Lat_Nwp_Idxx)
!---------------------------------------------------------------------------
 function INTERPOLATE_NWP_R4(z1,z2,z3,z4,lonx,latx) result(z)
  real, intent(in):: z1
  real, intent(in):: z2
  real, intent(in):: z3
  real, intent(in):: z4
  real, intent(in):: lonx
  real, intent(in):: latx
  real:: z
  !--- linear inteprpolation scheme
  if (minval((/z1,z2,z3,z4/)) == Missing_Value_Real4) then
    z = z1
  else
    z =  (1.0-lonx) * ((1.0-latx) * z1 + (latx)* z3) + &
         (lonx) * ((1.0-latx) * z2 + (latx)* z4)
  endif
 end function INTERPOLATE_NWP_R4
 function INTERPOLATE_NWP_I4(z1,z2,z3,z4,lonx,latx) result(z)
  integer(kind=int4), intent(in):: z1
  integer(kind=int4), intent(in):: z2
  integer(kind=int4), intent(in):: z3
  integer(kind=int4), intent(in):: z4
  real, intent(in):: lonx
  real, intent(in):: latx
  integer(kind=int4):: z
  !--- linear inteprpolation scheme
  if (minval((/z1,z2,z3,z4/)) == Missing_Value_Int4) then
    z = z1
  else
    z =  (1.0-lonx) * ((1.0-latx) * z1 + (latx)* z3) + &
         (lonx) * ((1.0-latx) * z2 + (latx)* z4)
  endif
 end function INTERPOLATE_NWP_I4
 function INTERPOLATE_NWP_I2(z1,z2,z3,z4,lonx,latx) result(z)
  integer(kind=int2), intent(in):: z1
  integer(kind=int2), intent(in):: z2
  integer(kind=int2), intent(in):: z3
  integer(kind=int2), intent(in):: z4
  real, intent(in):: lonx
  real, intent(in):: latx
  integer(kind=int2):: z
  !--- linear inteprpolation scheme
  if (minval((/z1,z2,z3,z4/)) == Missing_Value_Int2) then
    z = z1
  else
    z =  (1.0-lonx) * ((1.0-latx) * z1 + (latx)* z3) + &
         (lonx) * ((1.0-latx) * z2 + (latx)* z4)
  endif
 end function INTERPOLATE_NWP_I2
 function INTERPOLATE_NWP_I1(z1,z2,z3,z4,lonx,latx) result(z)
  integer(kind=int1), intent(in):: z1
  integer(kind=int1), intent(in):: z2
  integer(kind=int1), intent(in):: z3
  integer(kind=int1), intent(in):: z4
  real, intent(in):: lonx
  real, intent(in):: latx
  integer(kind=int1):: z
  !--- linear inteprpolation scheme
  if (minval((/z1,z2,z3,z4/)) == Missing_Value_Int1) then
    z = z1
  else
    z =  (1.0-lonx) * ((1.0-latx) * z1 + (latx)* z3) + &
         (lonx) * ((1.0-latx) * z2 + (latx)* z4)
  endif
 end function INTERPOLATE_NWP_I1

!---------------------------------------------------------------------------
! subroutine INTERPOLATE_PROFILE(z1,z2,z3,z4,lonx,latx,z)
!
! description of arguments
! Lon_Nwp_Idx, Lat_Nwp_Idx - nwp indices of closest nwp cell
! Lon_Nwp_Idxx,Lat_Nwp_Idxx - nwp indices of nwp cells of diagonal of bounding box
! lonx - longitude weighting factor 
! latx = latitude weighting factor
! z1 = data(Lon_Nwp_Idx, Lat_Nwp_Idx)
! z2 = data(Lon_Nwp_Idxx,Lat_Nwp_Idx)
! z3 = data(Lon_Nwp_Idx,Lat_Nwp_Idxx)
! z4 = data(Lon_Nwp_Idxx,Lat_Nwp_Idxx)
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

 subroutine INTERPOLATE_PROFILE(z1,z2,z3,z4,lonx,latx,z)

  real, dimension(:), intent(in):: z1
  real, dimension(:), intent(in):: z2
  real, dimension(:), intent(in):: z3
  real, dimension(:), intent(in):: z4
  real, intent(in):: lonx
  real, intent(in):: latx
  real, dimension(:), intent(out):: z

  !--- linear inteprpolation scheme
  z =  (1.0-lonx) * ((1.0-latx) * z1 + (latx)* z3) + &
           (lonx) * ((1.0-latx) * z2 + (latx)* z4)

 end subroutine INTERPOLATE_PROFILE


!--------------------------------------------------------------------------
! subroutine INTERPOLATE_NWP_TZ_PROFILES(Elem_Idx,Line_Idx)
!
! spatially interpolate TZ profiles
!
!--------------------------------------------------------------------------
 subroutine INTERPOLATE_NWP_TZ_PROFILES(Elem_Idx,Line_Idx)

   integer, intent(in):: Elem_Idx,Line_Idx
   integer:: Ilev



   do Ilev = 1,Nlevels_Nwp

      T_Prof_Nwp_pix(Ilev) = INTERPOLATE_NWP(T_Prof_Nwp(Ilev,i_Nwp(Elem_Idx,Line_Idx),j_Nwp(Elem_Idx,Line_Idx)), &
                                             T_Prof_Nwp(Ilev,i_Nwp(Elem_Idx,Line_Idx),j_Nwp_x(Elem_Idx,Line_Idx)), &
                                             T_Prof_Nwp(Ilev,i_Nwp_x(Elem_Idx,Line_Idx),j_Nwp(Elem_Idx,Line_Idx)), &
                                             T_Prof_Nwp(Ilev,i_Nwp_x(Elem_Idx,Line_Idx),j_Nwp_x(Elem_Idx,Line_Idx)), &
                                             Lon_Nwp_fac(Elem_Idx,Line_Idx), Lat_Nwp_fac(Elem_Idx,Line_Idx))

      Z_Prof_Nwp_pix(Ilev) = INTERPOLATE_NWP(Z_Prof_Nwp(Ilev,i_Nwp(Elem_Idx,Line_Idx),j_Nwp(Elem_Idx,Line_Idx)), &
                                             Z_Prof_Nwp(Ilev,i_Nwp(Elem_Idx,Line_Idx),j_Nwp_x(Elem_Idx,Line_Idx)), &
                                             Z_Prof_Nwp(Ilev,i_Nwp_x(Elem_Idx,Line_Idx),j_Nwp(Elem_Idx,Line_Idx)), &
                                             Z_Prof_Nwp(Ilev,i_Nwp_x(Elem_Idx,Line_Idx),j_Nwp_x(Elem_Idx,Line_Idx)), &
                                             Lon_Nwp_fac(Elem_Idx,Line_Idx), Lat_Nwp_fac(Elem_Idx,Line_Idx))
   enddo


 end subroutine INTERPOLATE_NWP_TZ_PROFILES

!-----------------------------------------------------------------------------
! subroutine CREATE_NWP_ARRAYS()
!
! allocate and initialize the memory needed for the nwp arrays
!-----------------------------------------------------------------------------
  subroutine CREATE_NWP_ARRAYS()

    integer:: Lon_Nwp_Idx, Lat_Nwp_Idx

!   Allocate arrays

    allocate(Mask_Nwp(Nlon_Nwp, Nlat_Nwp))

    allocate(bad_Nwp_mask(Nlon_Nwp, Nlat_Nwp))

    allocate(sfc_type_Nwp(Nlon_Nwp, Nlat_Nwp))
    allocate(Satzen_Nwp(Nlon_Nwp, Nlat_Nwp))
    allocate(Solzen_Nwp(Nlon_Nwp, Nlat_Nwp))

    allocate(Pmsl_Nwp(Nlon_Nwp, Nlat_Nwp))
    allocate(Psfc_Nwp(Nlon_Nwp, Nlat_Nwp))
    allocate(Zsfc_Nwp(Nlon_Nwp, Nlat_Nwp))
    allocate(Tmpsfc_Nwp(Nlon_Nwp, Nlat_Nwp))
    allocate(Tmpsfc_Nwp_Before(Nlon_Nwp, Nlat_Nwp))
    allocate(Tmpsfc_Nwp_After(Nlon_Nwp, Nlat_Nwp))
    allocate(Tmpair_Nwp(Nlon_Nwp, Nlat_Nwp))
    allocate(Tmpair_uni_Nwp(Nlon_Nwp, Nlat_Nwp))
    allocate(Tmpsfc_uni_Nwp(Nlon_Nwp, Nlat_Nwp))
    allocate(T_Trop_Nwp(Nlon_Nwp, Nlat_Nwp))
    allocate(P_Trop_Nwp(Nlon_Nwp, Nlat_Nwp))
    allocate(Rhsfc_Nwp(Nlon_Nwp, Nlat_Nwp))
    allocate(Uth_Nwp(Nlon_Nwp, Nlat_Nwp))
    allocate(hght500_Nwp(Nlon_Nwp, Nlat_Nwp))
    allocate(Tpw_Nwp(Nlon_Nwp, Nlat_Nwp))
    allocate(Ozone_Nwp(Nlon_Nwp, Nlat_Nwp))
    allocate(Weasd_Nwp(Nlon_Nwp, Nlat_Nwp))
    allocate(U_Wnd_10m_Nwp(Nlon_Nwp, Nlat_Nwp))
    allocate(V_Wnd_10m_Nwp(Nlon_Nwp, Nlat_Nwp))
    allocate(Wnd_Spd_10m_Nwp(Nlon_Nwp, Nlat_Nwp))
    allocate(Wnd_Dir_10m_Nwp(Nlon_Nwp, Nlat_Nwp))
    allocate(P_Std_Nwp(Nlevels_Nwp))
    allocate(Lat_Nwp(Nlat_Nwp))
    allocate(Lon_Nwp(Nlon_Nwp))
    allocate(Z_Prof_Nwp(Nlevels_Nwp, Nlon_Nwp, Nlat_Nwp))
    allocate(T_Prof_Nwp(Nlevels_Nwp, Nlon_Nwp, Nlat_Nwp))
    allocate(Ozone_Prof_Nwp(Nlevels_Nwp, Nlon_Nwp, Nlat_Nwp))
    allocate(Rh_Prof_Nwp(Nlevels_Nwp, Nlon_Nwp, Nlat_Nwp))
    allocate(Tpw_Prof_Nwp(Nlevels_Nwp, Nlon_Nwp, Nlat_Nwp))
    allocate(Clwmr_Prof_Nwp(Nlevels_Nwp, Nlon_Nwp, Nlat_Nwp))
    allocate(U_Wnd_Prof_Nwp(Nlevels_Nwp, Nlon_Nwp, Nlat_Nwp))
    allocate(V_Wnd_Prof_Nwp(Nlevels_Nwp, Nlon_Nwp, Nlat_Nwp))
    allocate(Land_Nwp(Nlon_Nwp, Nlat_Nwp))
    allocate(Ice_Nwp(Nlon_Nwp, Nlat_Nwp))
    allocate(Sfc_Level_Nwp(Nlon_Nwp, Nlat_Nwp))
    allocate(Tropo_Level_Nwp(Nlon_Nwp, Nlat_Nwp))
    allocate(Inversion_Level_Nwp(Nlon_Nwp, Nlat_Nwp))
    allocate(Inversion_Level_Profile_Nwp(Nlevels_Nwp,Nlon_Nwp, Nlat_Nwp))
    allocate(Lifting_Condensation_Level_Height_Nwp(Nlon_Nwp, Nlat_Nwp))
    allocate(Convective_Condensation_Level_Height_Nwp(Nlon_Nwp, Nlat_Nwp))
    allocate(Freezing_Level_Height_Nwp(Nlon_Nwp, Nlat_Nwp))
    allocate(Upper_Limit_Water_Height_Nwp(Nlon_Nwp, Nlat_Nwp))
    allocate(K_Index_Nwp(Nlon_Nwp, Nlat_Nwp))
    allocate(Level850_Nwp(Nlon_Nwp, Nlat_Nwp))
    allocate(Level500_Nwp(Nlon_Nwp, Nlat_Nwp))
    allocate(Level700_Nwp(Nlon_Nwp, Nlat_Nwp))
    allocate(Pc_Nwp(Nlon_Nwp, Nlat_Nwp))
    allocate(Sc_Lwp_Nwp(Nlon_Nwp, Nlat_Nwp))
    allocate(Lwp_Nwp(Nlon_Nwp, Nlat_Nwp))
    allocate(Iwp_Nwp(Nlon_Nwp, Nlat_Nwp))
    allocate(Cwp_Nwp(Nlon_Nwp, Nlat_Nwp))
    allocate(Cloud_Fraction_Satellite_Nwp(Nlon_Nwp, Nlat_Nwp))
    allocate(High_Cloud_Fraction_Satellite_Nwp(Nlon_Nwp, Nlat_Nwp))
    allocate(Mid_Cloud_Fraction_Satellite_Nwp(Nlon_Nwp, Nlat_Nwp))
    allocate(Low_Cloud_Fraction_Satellite_Nwp(Nlon_Nwp, Nlat_Nwp))
    allocate(Ncld_Layers_Nwp(Nlon_Nwp, Nlat_Nwp))
    allocate(Cld_Type_Nwp(Nlon_Nwp, Nlat_Nwp))

    allocate(temp2d_Nwp_1(Nlon_Nwp, Nlat_Nwp))
    allocate(temp2d_Nwp_2(Nlon_Nwp, Nlat_Nwp))
    allocate(temp3d_Nwp_1(Nlevels_Nwp, Nlon_Nwp, Nlat_Nwp))
    allocate(temp3d_Nwp_2(Nlevels_Nwp, Nlon_Nwp, Nlat_Nwp))

    allocate(T_Prof_Nwp_pix(Nlevels_Nwp))
    allocate(Z_Prof_Nwp_pix(Nlevels_Nwp))

    if (REFORMAT_GFS_ZXY == 1) then
     allocate(temp3d(Nlon_Nwp, Nlat_Nwp,Nlevels_Nwp))
    else
     allocate(temp3d(Nlevels_Nwp,Nlon_Nwp, Nlat_Nwp))
    endif

    Mask_Nwp = 0
    Bad_Nwp_Mask = sym%NO
    Sfc_Type_Nwp = 0
    Satzen_Nwp = 0
    Solzen_Nwp = 0
    Pmsl_Nwp = 0
    Psfc_Nwp = 0
    Zsfc_Nwp = 0
    Tmpsfc_Nwp = 0
    Tmpsfc_Nwp_Before = 0
    Tmpsfc_Nwp_After = 0
    Tmpair_Nwp = 0
    Tmpair_Uni_Nwp = 0
    Tmpsfc_Uni_Nwp = 0
    T_Trop_Nwp = 0
    P_Trop_Nwp = 0
    Rhsfc_Nwp = 0
    Uth_Nwp = 0
    Hght500_Nwp = 0
    Tpw_Nwp = 0
    Ozone_Nwp = 0
    Weasd_Nwp = 0
    U_Wnd_10m_Nwp = 0
    V_Wnd_10m_Nwp = 0
    Wnd_Spd_10m_Nwp = 0
    Wnd_Dir_10m_Nwp = 0
    Lon_Nwp = 0
    Lat_Nwp = 0
    P_Std_Nwp = 0
    Z_Prof_Nwp = 0
    T_Prof_Nwp = 0
    Rh_Prof_Nwp = 0.0
    Ozone_Prof_Nwp = 0.0
    Tpw_Prof_Nwp = 0.0
    Clwmr_Prof_Nwp = 0.0
    U_Wnd_Prof_Nwp = 0.0
    V_Wnd_Prof_Nwp = 0.0
    Land_Nwp = 0
    ice_Nwp = 0.0
    temp3d_Nwp_1 = 0
    temp3d_Nwp_2 = 0
    temp3d = 0
    Sfc_Level_Nwp = 0
    Tropo_Level_Nwp = 0
    Level850_Nwp = 0
    Level700_Nwp = 0
    Level500_Nwp = 0
    Inversion_Level_Nwp = 0
    Inversion_Level_Profile_Nwp = 0
    Lifting_Condensation_Level_Height_Nwp = 0
    Convective_Condensation_Level_Height_Nwp = 0
    Freezing_Level_Height_Nwp = 0
    Upper_Limit_Water_Height_Nwp = 0
    K_Index_Nwp = 0
    Pc_Nwp = 0.0
    Sc_Lwp_Nwp = 0.0
    Iwp_Nwp = 0.0
    Lwp_Nwp = 0.0
    Cwp_Nwp = 0.0
    Cloud_Fraction_Satellite_Nwp = 0.0
    High_Cloud_Fraction_Satellite_Nwp = 0.0
    Mid_Cloud_Fraction_Satellite_Nwp = 0.0
    Low_Cloud_Fraction_Satellite_Nwp = 0.0
    Ncld_Layers_Nwp = 0
    Cld_Type_Nwp = 0

    T_Prof_Nwp_pix = 0
    Z_Prof_Nwp_pix = 0

!--- create nwp lat and lon vectors
    do Lat_Nwp_Idx = 1, Nlat_Nwp
       Lat_Nwp(Lat_Nwp_Idx) = lat1_Nwp + (Lat_Nwp_Idx-1) * dLat_Nwp
    end do

    do Lon_Nwp_Idx = 1, Nlon_Nwp
       Lon_Nwp(Lon_Nwp_Idx) = lon1_Nwp + (Lon_Nwp_Idx-1) * dLon_Nwp
       if (Lon_Nwp(Lon_Nwp_Idx) > 180.0) then
         Lon_Nwp(Lon_Nwp_Idx) = Lon_Nwp(Lon_Nwp_Idx) - 360.0
       endif
    end do

end subroutine CREATE_NWP_ARRAYS


!-----------------------------------------------------------------------------
! subroutine DESTROY_NWP_ARRAYS
!
! deallocate the memory needed for the nwp arrays
!-----------------------------------------------------------------------------
subroutine DESTROY_NWP_ARRAYS

    if (allocated(Mask_Nwp))          deallocate(Mask_Nwp)
    if (allocated(bad_Nwp_mask))      deallocate(bad_Nwp_mask)
    if (allocated(sfc_type_Nwp))      deallocate(sfc_type_Nwp)
    if (allocated(Satzen_Nwp))        deallocate(Satzen_Nwp)
    if (allocated(Solzen_Nwp))        deallocate(Solzen_Nwp)

    if (allocated(Land_Nwp))           deallocate(Land_Nwp)
    if (allocated(ice_Nwp))            deallocate(ice_Nwp)
    if (allocated(Sfc_Level_Nwp))      deallocate(Sfc_Level_Nwp)
    if (allocated(Tropo_Level_Nwp))    deallocate(Tropo_Level_Nwp)
    if (allocated(Level850_Nwp))       deallocate(Level850_Nwp)
    if (allocated(Level700_Nwp))       deallocate(Level700_Nwp)
    if (allocated(Level500_Nwp))       deallocate(Level500_Nwp)
    if (allocated(Inversion_Level_Nwp)) deallocate(Inversion_Level_Nwp)
    if (allocated(Inversion_Level_Profile_Nwp)) deallocate(Inversion_Level_Profile_Nwp)
    if (allocated(Lifting_Condensation_Level_Height_Nwp))    deallocate(Lifting_Condensation_Level_Height_Nwp)
    if (allocated(Convective_Condensation_Level_Height_Nwp))    deallocate(Convective_Condensation_Level_Height_Nwp)
    if (allocated(Freezing_Level_Height_Nwp))    deallocate(Freezing_Level_Height_Nwp)
    if (allocated(Upper_Limit_Water_Height_Nwp))    deallocate(Upper_Limit_Water_Height_Nwp)
    if (allocated(K_index_Nwp))       deallocate(K_Index_Nwp)
    if (allocated(Pmsl_Nwp))          deallocate(Pmsl_Nwp)
    if (allocated(Psfc_Nwp))          deallocate(Psfc_Nwp)
    if (allocated(Zsfc_Nwp))          deallocate(Zsfc_Nwp)
    if (allocated(Tmpsfc_Nwp))        deallocate(Tmpsfc_Nwp)
    if (allocated(Tmpsfc_Nwp_Before)) deallocate(Tmpsfc_Nwp_Before)
    if (allocated(Tmpsfc_Nwp_After))  deallocate(Tmpsfc_Nwp_After)
    if (allocated(Tmpair_Nwp))        deallocate(Tmpair_Nwp)
    if (allocated(Tmpair_uni_Nwp))    deallocate(Tmpair_uni_Nwp)
    if (allocated(Tmpsfc_uni_Nwp))    deallocate(Tmpsfc_uni_Nwp)
    if (allocated(T_Trop_Nwp))        deallocate(T_Trop_Nwp)
    if (allocated(P_Trop_Nwp))        deallocate(P_Trop_Nwp)
    if (allocated(Rhsfc_Nwp))         deallocate(Rhsfc_Nwp)
    if (allocated(hght500_Nwp))       deallocate(hght500_Nwp)
    if (allocated(Uth_Nwp))           deallocate(Uth_Nwp)
    if (allocated(Tpw_Nwp))           deallocate(Tpw_Nwp)
    if (allocated(Ozone_Nwp))         deallocate(Ozone_Nwp)
    if (allocated(Weasd_Nwp))         deallocate(Weasd_Nwp)
    if (allocated(U_Wnd_10m_Nwp))     deallocate(U_Wnd_10m_Nwp)
    if (allocated(V_Wnd_10m_Nwp))     deallocate(V_Wnd_10m_Nwp)
    if (allocated(Wnd_Spd_10m_Nwp))   deallocate(Wnd_Spd_10m_Nwp)
    if (allocated(Wnd_Dir_10m_Nwp))   deallocate(Wnd_Dir_10m_Nwp)
    if (allocated(Lat_Nwp))           deallocate(Lat_Nwp)
    if (allocated(Lon_Nwp))           deallocate(Lon_Nwp)
    if (allocated(P_Std_Nwp))         deallocate(P_Std_Nwp)
    if (allocated(Z_Prof_Nwp))        deallocate(Z_Prof_Nwp)
    if (allocated(T_Prof_Nwp))        deallocate(T_Prof_Nwp)
    if (allocated(Rh_Prof_Nwp))       deallocate(Rh_Prof_Nwp)
    if (allocated(Ozone_Prof_Nwp))    deallocate(Ozone_Prof_Nwp)
    if (allocated(Tpw_Prof_Nwp))      deallocate(Tpw_Prof_Nwp)
    if (allocated(Clwmr_Prof_Nwp))    deallocate(Clwmr_Prof_Nwp)
    if (allocated(U_Wnd_Prof_Nwp))    deallocate(U_Wnd_Prof_Nwp)
    if (allocated(V_Wnd_Prof_Nwp))    deallocate(V_Wnd_Prof_Nwp)
    if (allocated(temp2d_Nwp_1))      deallocate(temp2d_Nwp_1)
    if (allocated(temp2d_Nwp_2))      deallocate(temp2d_Nwp_2)
    if (allocated(temp3d_Nwp_1))      deallocate(temp3d_Nwp_1)
    if (allocated(temp3d_Nwp_2))      deallocate(temp3d_Nwp_2)
    if (allocated(temp3d))            deallocate(temp3d)
    if (allocated(T_Prof_Nwp_pix))    deallocate(T_Prof_Nwp_pix)
    if (allocated(Z_Prof_Nwp_pix))    deallocate(Z_Prof_Nwp_pix)
    if (allocated(Pc_Nwp))            deallocate(Pc_Nwp)
    if (allocated(Sc_Lwp_Nwp))        deallocate(Sc_Lwp_Nwp)
    if (allocated(Lwp_Nwp))           deallocate(Lwp_Nwp)
    if (allocated(Iwp_Nwp))           deallocate(Iwp_Nwp)
    if (allocated(Cwp_Nwp))           deallocate(Cwp_Nwp)
    if (allocated(Cloud_Fraction_Satellite_Nwp))       deallocate(Cloud_Fraction_Satellite_Nwp)
    if (allocated(High_Cloud_Fraction_Satellite_Nwp))  deallocate(High_Cloud_Fraction_Satellite_Nwp)
    if (allocated(Mid_Cloud_Fraction_Satellite_Nwp))   deallocate(Mid_Cloud_Fraction_Satellite_Nwp)
    if (allocated(Low_Cloud_Fraction_Satellite_Nwp))   deallocate(Low_Cloud_Fraction_Satellite_Nwp)
    if (allocated(Ncld_Layers_Nwp))   deallocate(Ncld_Layers_Nwp)
    if (allocated(Cld_Type_Nwp))   deallocate(Cld_Type_Nwp)

  end subroutine DESTROY_NWP_ARRAYS

!----------------------------------------------------------------------
! subroutine KNOWING_P_COMPUTE_T_Z_NWP(Lon_Nwp_Idx,Lat_Nwp_Idx,P,T,Z,Ilev)
!
! derive height and temperature from a profile knowing pressure
!----------------------------------------------------------------------
 subroutine KNOWING_P_COMPUTE_T_Z_NWP(Lon_Nwp_Idx,Lat_Nwp_Idx,P,T,Z,Ilev)
  integer, intent(in):: Lon_Nwp_Idx, Lat_Nwp_Idx
  real, intent(in):: P
  real, intent(out):: T,Z
  integer, intent(out):: Ilev
  real:: dp, dt, dz


!--- interpolate pressure profile
  call LOCATE(P_Std_Nwp,Nlevels_Nwp,P,Ilev)
  Ilev = max(1,min(Nlevels_Nwp-1,Ilev))

  dp = P_Std_Nwp(Ilev+1) - P_Std_Nwp(Ilev)
  dt = T_Prof_Nwp(Ilev+1,Lon_Nwp_Idx,Lat_Nwp_Idx) - T_Prof_Nwp(Ilev,Lon_Nwp_Idx,Lat_Nwp_Idx)
  dz = Z_Prof_Nwp(Ilev+1,Lon_Nwp_Idx,Lat_Nwp_Idx) - Z_Prof_Nwp(Ilev,Lon_Nwp_Idx,Lat_Nwp_Idx)

!--- perform interpolation
  if (dp /= 0.0) then
   T = T_Prof_Nwp(Ilev, Lon_Nwp_Idx, Lat_Nwp_Idx) + dt/dp * (P - P_Std_Nwp(Ilev))
   Z = Z_Prof_Nwp(Ilev, Lon_Nwp_Idx, Lat_Nwp_Idx) + dz/dp * (P - P_Std_Nwp(Ilev))
  else
   T = T_Prof_Nwp(Ilev, Lon_Nwp_Idx, Lat_Nwp_Idx) 
   Z = Z_Prof_Nwp(Ilev, Lon_Nwp_Idx, Lat_Nwp_Idx)
  endif
 end subroutine KNOWING_P_COMPUTE_T_Z_NWP

!----------------------------------------------------------------------
! subroutine KNOWING_Z_COMPUTE_T_P_NWP(Lon_Nwp_Idx,Lat_Nwp_Idx,P,T,Z,Ilev)
!
! derive pressure and temperature from a profile knowing height
!----------------------------------------------------------------------
 subroutine KNOWING_Z_COMPUTE_T_P_NWP(Lon_Nwp_Idx,Lat_Nwp_Idx,P,T,Z,Ilev)
  integer, intent(in):: Lon_Nwp_Idx, Lat_Nwp_Idx
  real, intent(in):: Z
  real, intent(out):: T,P
  integer, intent(out):: Ilev
  real:: dp
  real:: dt
  real:: dz
  integer:: kstart
  integer:: kend
  integer:: Nlevels_temp

  !---
  kstart = Tropo_Level_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx)
  kend = Sfc_Level_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx)
  Nlevels_temp = kend - kstart + 1

  !--- interpolate pressure profile
  call LOCATE(Z_Prof_Nwp(kstart:kend,Lon_Nwp_Idx,Lat_Nwp_Idx),Nlevels_temp,Z,Ilev)
  Ilev = max(1,min(Nlevels_temp-1,Ilev))
  Ilev = Ilev + kstart - 1

  dp = P_Std_Nwp(Ilev+1) - P_Std_Nwp(Ilev)
  dt = T_Prof_Nwp(Ilev+1,Lon_Nwp_Idx,Lat_Nwp_Idx) - T_Prof_Nwp(Ilev,Lon_Nwp_Idx,Lat_Nwp_Idx)
  dz = Z_Prof_Nwp(Ilev+1,Lon_Nwp_Idx,Lat_Nwp_Idx) - Z_Prof_Nwp(Ilev,Lon_Nwp_Idx,Lat_Nwp_Idx)

  !--- perform interpolation
  if (dp /= 0.0) then
   T = T_Prof_Nwp(Ilev, Lon_Nwp_Idx, Lat_Nwp_Idx) + dt/dz * (Z - Z_Prof_Nwp(Ilev,Lon_Nwp_Idx,Lat_Nwp_Idx))
   P = P_Std_Nwp(Ilev) + dp/dz * (Z - Z_Prof_Nwp(Ilev,Lon_Nwp_Idx,Lat_Nwp_Idx))
  else
   T = T_Prof_Nwp(Ilev, Lon_Nwp_Idx, Lat_Nwp_Idx) 
   P = P_Std_Nwp(Ilev)
  endif

 end subroutine KNOWING_Z_COMPUTE_T_P_NWP


!----------------------------------------------------------------------
! subroutine KNOWING_Z_COMPUTE_T_P_NWP(Lon_Nwp_Idx,Lat_Nwp_Idx,P,T,Z,Ilev)
!
! derive pressure and temperature from a profile knowing height
!----------------------------------------------------------------------
 subroutine KNOWING_Z_COMPUTE_T_P_NWP_ARBITRARY_LEVELS(Lon_Nwp_Idx,Lat_Nwp_Idx,P,T,Z, &
                                                       kstart,kend,Ilev)
  integer, intent(in):: Lon_Nwp_Idx, Lat_Nwp_Idx
  real, intent(in):: Z
  integer, intent(in):: kstart
  integer, intent(in):: kend
  real, intent(out):: T,P
  integer, intent(out):: Ilev
  real:: dp
  real:: dt
  real:: dz
  integer:: Nlevels_temp

   Nlevels_temp = kend - kstart + 1

!--- interpolate pressure profile
  call LOCATE(Z_Prof_Nwp(kstart:kend,Lon_Nwp_Idx,Lat_Nwp_Idx),Nlevels_temp,Z,Ilev)
  Ilev = max(1,min(Nlevels_temp-1,Ilev))
  Ilev = Ilev + kstart - 1

  dp = P_Std_Nwp(Ilev+1) - P_Std_Nwp(Ilev)
  dt = T_Prof_Nwp(Ilev+1,Lon_Nwp_Idx,Lat_Nwp_Idx) - T_Prof_Nwp(Ilev,Lon_Nwp_Idx,Lat_Nwp_Idx)
  dz = Z_Prof_Nwp(Ilev+1,Lon_Nwp_Idx,Lat_Nwp_Idx) - Z_Prof_Nwp(Ilev,Lon_Nwp_Idx,Lat_Nwp_Idx)

!--- perform interpolation
  if (dp /= 0.0) then
   T = T_Prof_Nwp(Ilev, Lon_Nwp_Idx, Lat_Nwp_Idx) + dt/dz * (Z - Z_Prof_Nwp(Ilev,Lon_Nwp_Idx,Lat_Nwp_Idx))
   P = P_Std_Nwp(Ilev) + dp/dz * (Z - Z_Prof_Nwp(Ilev,Lon_Nwp_Idx,Lat_Nwp_Idx))
  else
   T = T_Prof_Nwp(Ilev, Lon_Nwp_Idx, Lat_Nwp_Idx)
   P = P_Std_Nwp(Ilev)
  endif

 end subroutine KNOWING_Z_COMPUTE_T_P_NWP_ARBITRARY_LEVELS
!----------------------------------------------------------------------------
! subroutine KNOWING_T_COMPUTE_P_Z_NWP(Lon_Nwp_Idx,Lat_Nwp_Idx,P,T,Z,klev,i,j,ierr)
!
! Interpolate P and Z from a profile knowing T
!
! klev = result is between klev and llev + 1
! 
! note, the highest level allowed is the tropopause
!
!----------------------------------------------------------------------------
 subroutine KNOWING_T_COMPUTE_P_Z_NWP(Lon_Nwp_Idx,Lat_Nwp_Idx,P,T,Z,klev,Elem_Idx,Line_Idx,ierr,z_interp_weight)

  integer, intent(in):: Lon_Nwp_Idx, Lat_Nwp_Idx
  real, intent(in):: T
  integer, intent(in):: Elem_Idx,Line_Idx
  real, intent(out):: P,Z
  integer, intent(out):: klev,ierr
  integer:: Nlevels_temp,kstart,kend,Level_within_Inversion_flag
  real:: dp, dt, dz,z1,z2,t1,t2
  real,optional,intent(out):: z_interp_weight

  ierr = sym%NO
  Z = Missing_Value_Real4
  P = Missing_Value_Real4

  !--- test for existence of a valid solution with troposphere
  kstart = Tropo_Level_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx)
  kend = Sfc_Level_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx)
  Nlevels_temp = kend - kstart + 1

  !--- interpolate temperature profile

  !--- check to see if warmer than max, than assume at surface
  if (T > maxval(T_Prof_Nwp(kstart:kend,Lon_Nwp_Idx,Lat_Nwp_Idx))) then
    P = P_Std_Nwp(kend)
    Z = Z_Prof_Nwp(kend,Lon_Nwp_Idx,Lat_Nwp_Idx)
    klev = kend - 1
    ierr = sym%NO
       if (present(z_interp_weight)) z_interp_weight = 0.
    return
  endif

  !--- check to see if colder than min, than assume tropopause
  if (T < minval(T_Prof_Nwp(kstart:kend,Lon_Nwp_Idx,Lat_Nwp_Idx))) then
    P = P_Std_Nwp(kstart)
    Z = Z_Prof_Nwp(kstart,Lon_Nwp_Idx,Lat_Nwp_Idx)
    klev = kstart + 1
    ierr = sym%NO
       if (present(z_interp_weight)) z_interp_weight = 0.
    return
  endif
  
  !--- if there is an Inversion, look below first
  Level_within_Inversion_flag = 0
  if (Inversion_Level_Nwp(Lon_Nwp_Idx, Lat_Nwp_Idx) > 0) then
   kstart = Inversion_Level_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx)
   kend = Sfc_Level_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx)
   Nlevels_temp = kend - kstart + 1
   call LOCATE(T_Prof_Nwp(kstart:kend,Lon_Nwp_Idx,Lat_Nwp_Idx),Nlevels_temp,T,klev)
   if ((klev > 0) .and. (klev < Nlevels_temp -1)) then
     klev = klev + kstart - 1
     Level_within_Inversion_flag = 1
   endif
  endif

  !--- if no solution within an Inversion, look above
  if (Level_within_Inversion_flag == 0) then

    kstart = Tropo_Level_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx)
    kend = Sfc_Level_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx)
    Nlevels_temp = kend - kstart + 1
    call LOCATE(T_Prof_Nwp(kstart:kend,Lon_Nwp_Idx,Lat_Nwp_Idx),Nlevels_temp,T,klev)
    klev = klev + kstart - 1
    klev = max(1,min(Nlevels_Nwp-1,klev))

  endif

  !-- if solution is above trop, set to trop values
  if (klev < Tropo_Level_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx)) then

    P = P_Std_Nwp(Tropo_Level_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx))
    Z = Z_Prof_Nwp(Tropo_Level_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx),Lon_Nwp_Idx,Lat_Nwp_Idx)

  else
    !--- determine derivatives
    dp = P_Std_Nwp(klev+1) - P_Std_Nwp(klev)

    if (Smooth_nwp_flag == sym%NO) then
     z2 = Z_Prof_Nwp(klev+1,Lon_Nwp_Idx,Lat_Nwp_Idx)
     z1 = Z_Prof_Nwp(klev,Lon_Nwp_Idx,Lat_Nwp_Idx)
     t2 = T_Prof_Nwp(klev+1,Lon_Nwp_Idx,Lat_Nwp_Idx)
     t1 = T_Prof_Nwp(klev,Lon_Nwp_Idx,Lat_Nwp_Idx)
    else
     z1 = INTERPOLATE_NWP(Z_Prof_Nwp(klev,i_Nwp(Elem_Idx,Line_Idx),j_Nwp(Elem_Idx,Line_Idx)), &
                          Z_Prof_Nwp(klev,i_Nwp_x(Elem_Idx,Line_Idx),j_Nwp(Elem_Idx,Line_Idx)), &
                          Z_Prof_Nwp(klev,i_Nwp(Elem_Idx,Line_Idx),j_Nwp_x(Elem_Idx,Line_Idx)), &
                          Z_Prof_Nwp(klev,i_Nwp_x(Elem_Idx,Line_Idx),j_Nwp_x(Elem_Idx,Line_Idx)), &
                          Lon_Nwp_fac(Elem_Idx,Line_Idx), Lat_Nwp_fac(Elem_Idx,Line_Idx))
     z2 = INTERPOLATE_NWP(Z_Prof_Nwp(klev+1,i_Nwp(Elem_Idx,Line_Idx),j_Nwp(Elem_Idx,Line_Idx)), &
                          Z_Prof_Nwp(klev+1,i_Nwp_x(Elem_Idx,Line_Idx),j_Nwp(Elem_Idx,Line_Idx)), &
                          Z_Prof_Nwp(klev+1,i_Nwp(Elem_Idx,Line_Idx),j_Nwp_x(Elem_Idx,Line_Idx)), &
                          Z_Prof_Nwp(klev+1,i_Nwp_x(Elem_Idx,Line_Idx),j_Nwp_x(Elem_Idx,Line_Idx)), &
                          Lon_Nwp_fac(Elem_Idx,Line_Idx), Lat_Nwp_fac(Elem_Idx,Line_Idx))
     t1 = INTERPOLATE_NWP(T_Prof_Nwp(klev,i_Nwp(Elem_Idx,Line_Idx),j_Nwp(Elem_Idx,Line_Idx)), &
                          T_Prof_Nwp(klev,i_Nwp_x(Elem_Idx,Line_Idx),j_Nwp(Elem_Idx,Line_Idx)), &
                          T_Prof_Nwp(klev,i_Nwp(Elem_Idx,Line_Idx),j_Nwp_x(Elem_Idx,Line_Idx)), &
                          T_Prof_Nwp(klev,i_Nwp_x(Elem_Idx,Line_Idx),j_Nwp_x(Elem_Idx,Line_Idx)), &
                          Lon_Nwp_fac(Elem_Idx,Line_Idx), Lat_Nwp_fac(Elem_Idx,Line_Idx))
     t2 = INTERPOLATE_NWP(T_Prof_Nwp(klev+1,i_Nwp(Elem_Idx,Line_Idx),j_Nwp(Elem_Idx,Line_Idx)), &
                          T_Prof_Nwp(klev+1,i_Nwp_x(Elem_Idx,Line_Idx),j_Nwp(Elem_Idx,Line_Idx)), &
                          T_Prof_Nwp(klev+1,i_Nwp(Elem_Idx,Line_Idx),j_Nwp_x(Elem_Idx,Line_Idx)), &
                          T_Prof_Nwp(klev+1,i_Nwp_x(Elem_Idx,Line_Idx),j_Nwp_x(Elem_Idx,Line_Idx)), &
                          Lon_Nwp_fac(Elem_Idx,Line_Idx), Lat_Nwp_fac(Elem_Idx,Line_Idx))

    endif

    dz = z2 - z1
    dt = t2 - t1

    !--- perform interpolation
    if (dt /= 0.0) then
     P = P_Std_Nwp(klev) + dp/dt * (T - t1)
     Z = z1 + dz/dt * (T - t1)
    else
     P = P_Std_Nwp(klev) 
     Z = z1
    endif
 

     if (present(z_interp_weight)) z_interp_weight =  (z - z1)/dz

  !--- end above trop check
  endif

 end subroutine KNOWING_T_COMPUTE_P_Z_NWP

!-------------------------------------------------------------
! subroutine FIND_NWP_LEVELS(Lon_Nwp_Idx,Lat_Nwp_Idx)
!
! subroutine to find key levels in the profiles
!-------------------------------------------------------------
subroutine FIND_NWP_LEVELS(Lon_Nwp_Idx,Lat_Nwp_Idx)

 integer (kind=int4), intent(in):: Lon_Nwp_Idx, Lat_Nwp_Idx
 integer (kind=int4):: k
 real:: P_Trop_temp
 real:: dT_dZ
 real:: es
 real:: e
 real:: T
 real:: Td
 real:: T850
 real:: T700
 real:: T500
 real:: Td850
 real:: Td700

   !--------------------------------------------------------------------
   !--- find surface level (standard closest but less than sfc pressure)
   !--------------------------------------------------------------------
   Sfc_Level_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx) = Nlevels_Nwp
   do k = Nlevels_Nwp,1,-1
       if (P_Std_Nwp(k) < Psfc_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx)) then
            Sfc_Level_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx) = k
            exit
       endif
    enddo

   !--------------------------------------------------------------------
   !--- find some standard levels 
   !--------------------------------------------------------------------
   Level850_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx) = 0
   Level700_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx) = 0
   Level500_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx) = 0
   do k = 1, Sfc_Level_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx)-1
      if ((P_Std_Nwp(k) <= 850.0) .and. (P_Std_Nwp(k+1) > 850.0)) then
        Level850_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx) = k
      endif
      if ((P_Std_Nwp(k) <= 700.0) .and. (P_Std_Nwp(k+1) > 700.0)) then
        Level700_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx) = k
      endif
      if ((P_Std_Nwp(k) <= 500.0) .and. (P_Std_Nwp(k+1) > 500.0)) then
        Level500_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx) = k
      endif
   enddo

   !--------------------------------------------------------------------
   !--- find tropopause level  based on tropopause pressure
   !--- tropopause is between tropopause_level and tropopaue_level + 1
   !--------------------------------------------------------------------

   !--- constrain tropopause pressure to be greater than 75 mb
   P_Trop_Temp = max(P_Trop_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx),75.0)

   Tropo_Level_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx) = 1
   do k = 1, Sfc_Level_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx)-1
      if ((P_Std_Nwp(k) <= P_Trop_temp) .and. &
          (P_Std_Nwp(k+1) > P_Trop_temp)) then
        Tropo_Level_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx) = k
      endif
   enddo

    !--- check if tropopause level found
    if (Tropo_Level_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx) == 0) then
        Tropo_Level_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx) = 1   !assume top level if no trop found
    endif         

!--------------------------------------------------------------------
!--- find stratopause level starting at 500 mb by looking for levels
!--- with the minimum temp
!--- Note, no point in doing this until higher vertical resolution
!--- profiles are used.
!--------------------------------------------------------------------

   !---------------------------------------------------------------------
   ! Inversion Level Profile
   !---------------------------------------------------------------------
   Inversion_Level_Profile_Nwp(:,Lon_Nwp_Idx,Lat_Nwp_Idx) = sym%NO
   do k = Tropo_Level_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx),Sfc_Level_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx)-1
      if (T_Prof_Nwp(k,Lon_Nwp_Idx,Lat_Nwp_Idx) - T_Prof_Nwp(k+1,Lon_Nwp_Idx,Lat_Nwp_Idx) > Delta_T_Inversion) then
        Inversion_Level_Profile_Nwp(k,Lon_Nwp_Idx,Lat_Nwp_Idx) = sym%YES
      endif
   enddo

   !---------------------------------------------------------------------
   ! find Inversion level - highest level Inversion below Tropopause
   !---------------------------------------------------------------------
   Inversion_Level_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx) = 0
   do k = Tropo_Level_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx),Sfc_Level_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx)-1
      if ((T_Prof_Nwp(k,Lon_Nwp_Idx,Lat_Nwp_Idx) - T_Prof_Nwp(k+1,Lon_Nwp_Idx,Lat_Nwp_Idx) > Delta_T_Inversion) .and. &
          (P_Std_Nwp(k) >= P_Inversion_Min)) then
          Inversion_Level_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx) = k
          exit 
      endif
   enddo

   !-------------------------------------------------------------------
   ! Find the Height of the Freezing Level in the NWP Profiles
   ! Start at the Tropopause and work down
   !-------------------------------------------------------------------
   Freezing_Level_Height_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx) = Z_Prof_Nwp(Sfc_Level_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx),Lon_Nwp_Idx,Lat_Nwp_Idx)
   do k = Tropo_Level_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx),Sfc_Level_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx)-1

      if (T_Prof_Nwp(k,Lon_Nwp_Idx,Lat_Nwp_Idx) > 273.15) then

       if (Z_Prof_Nwp(k,Lon_Nwp_Idx,Lat_Nwp_Idx) > Z_Prof_Nwp(k-1,Lon_Nwp_Idx,Lat_Nwp_Idx)) then
         dT_dZ = (T_Prof_Nwp(k,Lon_Nwp_Idx,Lat_Nwp_Idx) - T_Prof_Nwp(k-1,Lon_Nwp_Idx,Lat_Nwp_Idx)) /  &
                 (Z_Prof_Nwp(k,Lon_Nwp_Idx,Lat_Nwp_Idx) - Z_Prof_Nwp(k-1,Lon_Nwp_Idx,Lat_Nwp_Idx))
         Freezing_Level_Height_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx) = (T_Prof_Nwp(k,Lon_Nwp_Idx,Lat_Nwp_Idx) - 273.15) /  &
                                                               dT_dZ + Z_Prof_Nwp(k,Lon_Nwp_Idx,Lat_Nwp_Idx)
       else
         Freezing_Level_Height_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx) = Z_Prof_Nwp(k-1,Lon_Nwp_Idx,Lat_Nwp_Idx)
       endif

       exit

      endif
   enddo

   !---------------------------------------------------------------------------
   ! Find the Height of the Lifting(LCL) and Convective Condensation Level(CCL) 
   ! Heights from the NWP Profiles
   !----------------------------------------------------------------------------
    Lifting_Condensation_Level_Height_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx) = Missing_Value_Real4
    Convective_Condensation_Level_Height_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx) = Missing_Value_Real4
    T = Tmpair_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx)    ! K
    if (T > 180.0) then
     if (T > 253.0) then
       es = VAPOR(T)    !saturation vapor pressure wrt water hpa
     else
       es = VAPOR_ICE(T) !saturation vapor pressure wrt ice hpa
     endif
     e = es * Rhsfc_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx) / 100.0  !vapor pressure in hPa
     Td = 273.15 + 243.5 * alog(e / 6.112)  / (17.67 - alog(e/6.112))      !Dewpoint T in K
     Lifting_Condensation_Level_Height_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx) = 0.125*(T - Td)  !km
     Convective_Condensation_Level_Height_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx) = (T - Td)/4.4  !km
   endif

   !---------------------------------------------------------------------------
   ! K Index 
   !----------------------------------------------------------------------------
   K_Index_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx) = Missing_Value_Real4
   if (Sfc_Level_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx) > 0 .and. &
       Level700_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx) > 0 .and. &
       Level500_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx) > 0) then

     if (Level850_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx) == 0) then
         k = Sfc_Level_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx)
     else
         k = Level850_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx)
     endif
     T850 = T_Prof_Nwp(k,Lon_Nwp_Idx,Lat_Nwp_Idx)
     T700 = T_Prof_Nwp(Level700_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx),Lon_Nwp_Idx,Lat_Nwp_Idx)
     T500 = T_Prof_Nwp(Level500_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx),Lon_Nwp_Idx,Lat_Nwp_Idx)
     if (T850 > 253.0) then
       es = VAPOR(T850)    !saturation vapor pressure wrt water hpa
      else
       es = VAPOR_ICE(T850) !saturation vapor pressure wrt ice hpa
     endif
     e = es * Rh_Prof_Nwp(k,Lon_Nwp_Idx,Lat_Nwp_Idx) / 100.0  !vapor pressure in hPa
     Td850 = 273.15 + 243.5 * alog(e / 6.112)  / (17.67 - alog(e/6.112))       !Dewpoint T in K

     if (T700 > 253.0) then
       es = VAPOR(T700)    !saturation vapor pressure wrt water hpa
     else
       es = VAPOR_ICE(T700) !saturation vapor pressure wrt ice hpa
     endif
     e = es * Rh_Prof_Nwp(Level700_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx),Lon_Nwp_Idx,Lat_Nwp_Idx) / 100.0  !vapor pressure in hPa
     Td700 = 273.15 + 243.5 * alog(e / 6.112)  / (17.67 - alog(e/6.112))       !Dewpoint T in K
    
     if (T700 /= Td700) then
        K_Index_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx) = (T850 - T500) + Td850 - (T700 - Td700) - 273.15
     endif
   endif

   !-------------------------------------------------------------------
   ! Find the Height above which to consider CWP to be 100% ice in the NWP Profiles
   ! Start at the Tropopause and work down
   !-------------------------------------------------------------------
   Upper_Limit_Water_Height_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx) &
      = Z_Prof_Nwp(Sfc_Level_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx),Lon_Nwp_Idx,Lat_Nwp_Idx)

   do k = Tropo_Level_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx),Sfc_Level_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx)-1

      if (T_Prof_Nwp(k,Lon_Nwp_Idx,Lat_Nwp_Idx) > 260.0) then

       if (Z_Prof_Nwp(k,Lon_Nwp_Idx,Lat_Nwp_Idx) > Z_Prof_Nwp(k-1,Lon_Nwp_Idx,Lat_Nwp_Idx)) then
         dT_dZ = (T_Prof_Nwp(k,Lon_Nwp_Idx,Lat_Nwp_Idx) - T_Prof_Nwp(k-1,Lon_Nwp_Idx,Lat_Nwp_Idx)) /  &
                 (Z_Prof_Nwp(k,Lon_Nwp_Idx,Lat_Nwp_Idx) - Z_Prof_Nwp(k-1,Lon_Nwp_Idx,Lat_Nwp_Idx))
         Upper_Limit_Water_Height_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx) = (T_Prof_Nwp(k,Lon_Nwp_Idx,Lat_Nwp_Idx) - 260.0) / &
                 dT_dZ + Z_Prof_Nwp(k,Lon_Nwp_Idx,Lat_Nwp_Idx)
       else
         Upper_Limit_Water_Height_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx) = Z_Prof_Nwp(k-1,Lon_Nwp_Idx,Lat_Nwp_Idx)
       endif

       exit

      endif
   enddo

end subroutine FIND_NWP_LEVELS

!======================================================================
! subroutine COMPUTE_COAST_MASK_NWP(j1,j2)
!
! Compute Coast Mask for NWP data
!
! dependencies - Land_Mask and land_Mask_Nwp should be already computed
!======================================================================
 subroutine COMPUTE_COAST_MASK_NWP(j1,j2)

  integer, intent(in):: j1, j2
  integer:: Elem_Idx,Line_Idx,Lon_Nwp_Idx,Lat_Nwp_Idx

  !--- loop over all pixels in segment
  do Line_Idx = j1,j1+j2-1
     do Elem_Idx = 1, Image%Number_Of_Elements

        !--- initialize to no
        Sfc%Coast_Mask_Nwp(Elem_Idx,Line_Idx) = sym%NO
       
        !--- check for valid pixels
        if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES)  then
          cycle
        endif

        !--- save nwp indices
        Lon_Nwp_Idx = i_Nwp(Elem_Idx,Line_Idx)
        Lat_Nwp_Idx = j_Nwp(Elem_Idx,Line_Idx)

        !--- derive nwp coast mask
        if (Sfc%Land_Mask(Elem_Idx,Line_Idx) /= Land_Nwp(Lon_Nwp_Idx,Lat_Nwp_Idx)) then
           Sfc%Coast_Mask_Nwp(Elem_Idx,Line_Idx) = sym%YES
        endif 

     !--- end loop over all pixels in segment
     enddo
   enddo

 end subroutine COMPUTE_COAST_MASK_NWP
!======================================================================
! subroutine COMPUTE_NWP_CLOUD_PARAMETERS(Clwmr_Profile,
!    temperature_Profile,Pressure_Profile,P_Cld_Top,iwp,lwp,cwp,
!    Number_of_Cloud_Layers)
!
! P_Cld_Top = Pressure of Cloud Top (tau = 1)
!======================================================================
 subroutine COMPUTE_NWP_CLOUD_PARAMETERS(Tropopause_Level_Nwp, &
                                         Surface_Level_Nwp, &
                                         Clwmr_Profile, &
                                         Rel_Humid_Profile, &
                                         Temperature_Profile, &
                                         Pressure_Profile, &
                                         P_Cld_Top,Iwp,Lwp,Sc_Lwp,Cwp,&
                                         Cloud_Fraction_Satellite, &
                                         High_Cloud_Fraction_Satellite, &
                                         Mid_Cloud_Fraction_Satellite, &
                                         Low_Cloud_Fraction_Satellite, &
                                         Number_Of_Cloud_Layers, &
                                         Cloud_Type)
   implicit none                                      

  integer(kind=int1), intent(in):: Tropopause_Level_Nwp
  integer(kind=int1), intent(in):: Surface_Level_Nwp
  real, dimension(:), intent(in):: Clwmr_Profile
  real, dimension(:), intent(in):: Rel_Humid_Profile
  real, dimension(:), intent(in):: Temperature_Profile
  real, dimension(:), intent(in):: Pressure_Profile
  real, intent(out):: P_Cld_Top
  real, intent(out):: Iwp
  real, intent(out):: Lwp
  real, intent(out):: Sc_Lwp
  real, intent(out):: Cwp
  integer (kind=int1), intent(out):: Number_Of_Cloud_Layers
  integer (kind=int1), intent(out):: Cloud_Type

  real, dimension(size(Clwmr_Profile)):: Cloud_Fraction_Profile
  real:: Sat_Specific_Humidity
  real:: Cloud_Fraction_Satellite
  real:: High_Cloud_Fraction_Satellite
  real:: Mid_Cloud_Fraction_Satellite
  real:: Low_Cloud_Fraction_Satellite
  real:: T_Cld_Top
 
  real:: Clwmr_Min
  real, parameter:: Lwp_Threshold  = 5.0
  real, parameter:: Iwp_Threshold  = 1.0
  real, parameter:: Frac_Min_Threshold  = 1.0e-06
  real, parameter:: Clwmr_Min_Threshold  = 0.0
  real, parameter:: Cwp_Min_Threshold  = 2.0     !g/m^2 - arbitrary need to investigate
  real, parameter:: Max_Temperature_Ice = 273.15
  real, parameter:: Min_Temperature_Water  = 253.15
  real, parameter:: Max_Temperature_Sc_Water = 273.15
  real, parameter:: Min_Temperature_Sc_Water  = 263.15
  integer:: Number_Of_Levels_In_Profile
  real:: Factor
  real:: Clwmr_Ice_Layer
  real:: Clwmr_Water_Layer
  real:: Clwmr_Sc_Water_Layer
  real:: Ice_Frac_Top
  real:: Ice_Frac_Bot
  real:: Sc_Water_Frac_Top
  real:: Sc_Water_Frac_Bot
  integer:: Ilay
  integer:: Ilev
  real, parameter:: k1 = 0.25
  real, parameter:: k2 = 100
  real, parameter:: k3 = 0.49
  real:: T
  real:: e
  real:: es

  Number_Of_Levels_In_Profile = size(Clwmr_Profile)

  if (Tropopause_Level_Nwp <= 0 .or. Surface_Level_Nwp <= 0) then
     return
  endif

  !---------------------------------------------------------------------
  ! convert clwmr profiles into optical depths
  ! iwp and lwp are in g/m^2
  !---------------------------------------------------------------------
  Lwp = 0.0
  Iwp = 0.0
  Cwp = 0.0
  Sc_Lwp = 0.0
  P_Cld_Top = Missing_Value_Real4
  T_Cld_Top = Missing_Value_Real4

     do Ilay = Tropopause_Level_Nwp-1 , Surface_Level_Nwp-1

        Ice_Frac_Top  = min(1.0,max(0.0,&
                           (Max_Temperature_Ice - Temperature_Profile(ilay))/ &
                           (Max_Temperature_Ice-Min_Temperature_Water)))

        Ice_Frac_Bot  = min(1.0,max(0.0,&
                           (Max_Temperature_Ice - Temperature_Profile(ilay+1))/ &
                           (Max_Temperature_Ice-Min_Temperature_Water)))

        Clwmr_Ice_Layer = 0.5 * (Ice_Frac_Top*Clwmr_Profile(ilay) +  &
                                 Ice_Frac_Bot*Clwmr_Profile(ilay+1))

        Clwmr_Water_Layer = 0.5 * ((1.0-Ice_Frac_Top)*Clwmr_Profile(ilay) + &
                                   (1.0-Ice_Frac_Bot)*Clwmr_Profile(ilay+1))

        Sc_Water_Frac_Top  = min(1.0,max(0.0,&
                           (Max_Temperature_Sc_Water - Temperature_Profile(ilay))/ &
                           (Max_Temperature_Sc_Water - Min_Temperature_Sc_Water)))

        Sc_Water_Frac_Bot  = min(1.0,max(0.0,&
                           (Max_Temperature_Sc_Water - Temperature_Profile(ilay+1))/ &
                           (Max_Temperature_Sc_Water - Min_Temperature_Sc_Water)))

        Clwmr_Sc_Water_Layer = 0.5 * (Sc_Water_Frac_Top*Clwmr_Profile(ilay) + &
                                     Sc_Water_Frac_Bot*Clwmr_Profile(ilay+1))

        Factor = 1000.0 * 100.0 * (Pressure_Profile(ilay+1) - Pressure_Profile(ilay)) / g

        Iwp = Iwp + Clwmr_Ice_Layer * Factor
        Lwp = Lwp + Clwmr_Water_Layer * Factor
        Sc_Lwp = Sc_Lwp + Clwmr_Sc_Water_Layer * Factor

        Cwp = Lwp + Iwp

        !--- compute P_Cld_Top
        if ((Cwp >= Iwp_Threshold) .and. (P_Cld_Top == Missing_Value_Real4)) then
            P_Cld_Top = P_Std_Nwp(ilay+1)
            T_Cld_Top = Temperature_Profile(ilay+1)
        endif

      enddo

  !----------------------------------------------------------------------
  ! test to see if sufficient amount of cloud water was found in the 
  ! column to estimate cloud fractions.  if not, exit
  !----------------------------------------------------------------------
  if (Cwp < Cwp_Min_Threshold) then
     Number_Of_Cloud_Layers = 0
     Cloud_Type = sym%CLEAR_TYPE
     Cloud_Fraction_Satellite = 0.0
     High_Cloud_Fraction_Satellite = 0.0
     Mid_Cloud_Fraction_Satellite = 0.0
     Low_Cloud_Fraction_Satellite = 0.0
     return
  endif

  !----------------------------------------------------------------------
  ! compute number of distinct cloud layers in GFS profile
  !----------------------------------------------------------------------
  Number_Of_Cloud_Layers = 0

  if (P_Cld_Top /= Missing_Value_Real4) then

         do ilay = 3, Surface_Level_Nwp - 1

            Clwmr_Min = 1.0e-05 * Pressure_Profile(ilay)/1000.0    !

            if (Clwmr_Profile(ilay) > Clwmr_Min .and.   &
                Clwmr_Profile(ilay-1) < Clwmr_Min .and. &
                Clwmr_Profile(ilay-2) < Clwmr_Min) then

                Number_Of_Cloud_Layers = Number_Of_Cloud_Layers + 1

            endif
        enddo
  endif 

  !----------------------------------------------------------------------
  ! compute cloud type
  !----------------------------------------------------------------------
  Cloud_Type = sym%CLEAR_TYPE

  if (T_Cld_Top /= Missing_Value_Real4) then

      Cloud_Type = sym%CIRRUS_TYPE
      if (T_Cld_Top > 258.0) Cloud_Type = sym%SUPERCOOLED_TYPE
      if (T_Cld_Top > 273.0) Cloud_Type = sym%WATER_TYPE

      if (Cloud_Type == sym%CIRRUS_TYPE) then
        if (Number_Of_Cloud_Layers > 1) Cloud_Type = sym%OVERLAP_TYPE
        if (Cwp > 100) Cloud_Type = sym%OPAQUE_ICE_TYPE
      endif

  endif

  !----------------------------------------------------------------------
  ! compute cloud fraction profile
  ! procedure described at http://www.emc.ncep.noaa.gov/GFS/doc.php
  ! Note, Clwmr_Profile has units of g/g
  !----------------------------------------------------------------------
  Cloud_Fraction_Profile = 0
  do Ilev = Tropopause_Level_Nwp-1, Surface_Level_Nwp

    if (Clwmr_Profile(Ilev) > Clwmr_Min_Threshold) then

     T = Temperature_Profile(ilev)
     if (T > 180.0) then
      if (T > 253.0) then
        es = VAPOR(T)    !saturation vapor pressure wrt water hpa
      else
        es = VAPOR_ICE(T) !saturation vapor pressure wrt ice hpa
      endif
      e = es * Rel_Humid_Profile(ilev) / 100.0  !vapor pressure in hPa

      Sat_Specific_Humidity = (0.622 * es) / (Pressure_Profile(ilev) - es)
  
      Cloud_Fraction_Profile(Ilev) = ((Rel_Humid_Profile(ilev)/100.0)**k1) * (1.0 - &
                              exp( (-1.0*k2*1000.0*Clwmr_Profile(ilev)) /  &
                                   (((1.0-Rel_Humid_Profile(ilev)/100.0)*Sat_Specific_Humidity)**k3)))

      if (Rel_Humid_Profile(Ilev) > 100.0) Cloud_Fraction_Profile(Ilev) = 1.0

     endif

    endif

  enddo

  !----------------------------------------------------------------------
  ! compute cloud fraction as seen by a satellite using Max-Random Overlap
  !----------------------------------------------------------------------
  Cloud_Fraction_Satellite = 0.0
  High_Cloud_Fraction_Satellite = 0.0
  Mid_Cloud_Fraction_Satellite = 0.0
  Low_Cloud_Fraction_Satellite = 0.0

  do Ilev = Tropopause_Level_Nwp-1, Surface_Level_Nwp

       if (Cloud_Fraction_Profile(Ilev) > Frac_Min_Threshold)  then

          if (Cloud_Fraction_Satellite == 0.0) then

               Cloud_Fraction_Satellite = Cloud_Fraction_Profile(Ilev)

          else
!           if (Cloud_Fraction_Profile(Ilev-1) > Frac_Min_Threshold .and. &
!             Cloud_Fraction_Profile(Ilev) > Frac_Min_Threshold) then       !max overlap

!             Cloud_Fraction_Satellite_Temp = max(Cloud_Fraction_Profile(Ilev-1),Cloud_Fraction_Profile(Ilev))
!             Cloud_Fraction_Satellite = max(Cloud_Fraction_Satellite,Cloud_Fraction_Satellite_Temp)

!           else    !random overlap

!             Cloud_Fraction_Satellite_Temp = 1.0 - (1.0-Cloud_Fraction_Satellite)*(1.0-Cloud_Fraction_Profile(Ilev))

!           endif

            Cloud_Fraction_Satellite = max(Cloud_Fraction_Satellite, Cloud_Fraction_Profile(Ilev))

          endif

          if (Pressure_Profile(Ilev) <= 440.0) then
               High_Cloud_Fraction_Satellite = Cloud_Fraction_Satellite
          elseif (Pressure_Profile(Ilev) > 440.0 .and. Pressure_Profile(Ilev) <= 680.0) then
               Mid_Cloud_Fraction_Satellite = Cloud_Fraction_Satellite - High_Cloud_Fraction_Satellite
          else
               Low_Cloud_Fraction_Satellite = Cloud_Fraction_Satellite &
                     - High_Cloud_Fraction_Satellite - Mid_Cloud_Fraction_Satellite
          endif

          if (Cloud_Fraction_Satellite == 1.0) then
             exit
          endif
               
       endif

!write (unit=6,fmt="(I4,10f8.5)") Ilev, Cloud_Fraction_Profile(Ilev), Cloud_Fraction_Satellite

  enddo

 end subroutine COMPUTE_NWP_CLOUD_PARAMETERS

!----------------------------------------------------------------------
! SUBROUTINE prof_lookup_using_p(Zlev, plev, Tlev, z, p, t, Ilev, a) 
!
!----------------------------------------------------------------------
SUBROUTINE prof_lookup_using_p(Zlev, plev, Tlev, z, p, t, Ilev, a)
   
  REAL (kind=real4), intent(in), dimension(:) :: Zlev, plev, Tlev
  REAL (kind=real4), intent(in) :: p
  REAL (kind=real4), intent(out) :: z, t
  INTEGER (kind=int4), intent(out) :: Ilev
  REAL (kind=real4), intent(out) :: a
  INTEGER (kind=int4) :: Nlev
  REAL (kind=real4) :: dp, dt, dz

  Nlev = size(plev,1)
   
  call LOCATE(plev,Nlev,p,Ilev)
  Ilev = max(1,min(Nlev-1,Ilev))

  dp = Plev(Ilev+1) - Plev(Ilev)
  dt = Tlev(Ilev+1) - Tlev(Ilev)
  dz = Zlev(Ilev+1) - Zlev(Ilev)

  if (dp /= 0.0) then
    a = (p - Plev(Ilev))/dp
    t = Tlev(Ilev) + a*dt
    z = Zlev(Ilev) + a*dz
  else
    a = 0.0
    t = Tlev(Ilev) 
    z = Zlev(Ilev)
  endif
 
END SUBROUTINE prof_lookup_using_p

!----------------------------------------------------------------------
! SUBROUTINE prof_lookup_using_t(Zlev, plev, Tlev, z, p, t, Ilev, a)
!----------------------------------------------------------------------
SUBROUTINE prof_lookup_using_t(Zlev, plev, Tlev, z, p, t, Ilev, a)
   
  REAL (kind=real4), intent(in), dimension(:) :: Zlev, plev, Tlev
  REAL (kind=real4), intent(in) :: t
  REAL (kind=real4), intent(out) :: p, z
  INTEGER (kind=int4), intent(out) :: Ilev
  REAL (kind=real4), intent(out) :: a
  INTEGER (kind=int4) :: Nlev

  Nlev = size(plev,1)
   
  if (t < Tlev(1)) then
    CALL prof_lookup_using_t_lapse(Zlev, plev, Tlev, z, p, t, Ilev, a)
  else
    CALL prof_lookup_using_T_Prof(Zlev, plev, Tlev, z, p, t, Ilev, a)
  endif
 
END SUBROUTINE prof_lookup_using_t

!----------------------------------------------------------------------
! SUBROUTINE prof_lookup_using_t_lapse(Zlev, plev, Tlev, z, p, t, Ilev, a)
!----------------------------------------------------------------------
SUBROUTINE prof_lookup_using_t_lapse(Zlev, plev, Tlev, z, p, t, Ilev, a)
   
  REAL (kind=real4), intent(in), dimension(:) :: Zlev, plev, Tlev
  REAL (kind=real4), intent(in) :: t
  REAL (kind=real4), intent(out) :: p, z
  INTEGER (kind=int4), intent(out) :: Ilev
  REAL (kind=real4), intent(out) :: a
  REAL (kind=real4) :: dt, dz, t_tmp
  REAL (kind=real4) :: DZ_STRATO_LAPSE_RATE
  REAL (kind=real4), PARAMETER :: DT_STRATO_LAPSE_RATE = -6.5 !k
  REAL (kind=real4), PARAMETER :: DZ_STRATO_LAPSE_RATE_M = 1000.0 !m
  REAL (kind=real4), PARAMETER :: DZ_STRATO_LAPSE_RATE_KM = 1.0 !km
 
 
  IF (zlev(1) > 1000.0) THEN
    DZ_STRATO_LAPSE_RATE = DZ_STRATO_LAPSE_RATE_M
  ELSE
    DZ_STRATO_LAPSE_RATE = DZ_STRATO_LAPSE_RATE_KM
  ENDIF


  dt = DT_STRATO_LAPSE_RATE
  dz = DZ_STRATO_LAPSE_RATE

  a = (t - Tlev(1))/dt
  z = Zlev(1) + a*dz
  
  !This is a hack, really should be using full profile
  CALL prof_lookup_using_z(Zlev, plev, Tlev, z, p, t_tmp, Ilev, a)   

END SUBROUTINE prof_lookup_using_t_lapse

!----------------------------------------------------------------------
! SUBROUTINE prof_lookup_using_T_Prof(Zlev, plev, Tlev, z, p, t, Ilev, a)
!----------------------------------------------------------------------

SUBROUTINE prof_lookup_using_T_Prof(Zlev, plev, Tlev, z, p, t, Ilev, a)
   
  REAL (kind=real4), intent(in), dimension(:) :: Zlev, plev, Tlev
  REAL (kind=real4), intent(in) :: t
  REAL (kind=real4), intent(out) :: p, z
  INTEGER (kind=int4), intent(out) :: Ilev
  REAL (kind=real4), intent(out) :: a
  INTEGER (kind=int4) :: Nlev
  REAL (kind=real4) :: dp, dt, dz

  Nlev = size(plev,1)
   
  call LOCATE(Tlev,Nlev,t,Ilev)
  Ilev = max(1,min(Nlev-1,Ilev))

  dp = Plev(Ilev+1) - Plev(Ilev)
  dt = Tlev(Ilev+1) - Tlev(Ilev)
  dz = Zlev(Ilev+1) - Zlev(Ilev)

  if (dp /= 0.0 .and. dz /= 0) then
    a = (t - Tlev(Ilev))/dt
    p = Plev(Ilev) + a*dp
    z = Zlev(Ilev) + a*dz
  else
    a = 0.0
    p = Plev(Ilev) 
    z = Zlev(Ilev)
  endif
 
END SUBROUTINE prof_lookup_using_T_Prof

!----------------------------------------------------------------------
! SUBROUTINE prof_lookup_using_z(Zlev, plev, Tlev, z, p, t, Ilev, a)
!----------------------------------------------------------------------

SUBROUTINE prof_lookup_using_z(Zlev, plev, Tlev, z, p, t, Ilev, a)
   
  REAL (kind=real4), intent(in), dimension(:) :: Zlev, plev, Tlev
  REAL (kind=real4), intent(in) :: z
  REAL (kind=real4), intent(out) :: p, t
  INTEGER (kind=int4), intent(out) :: Ilev
  REAL (kind=real4), intent(out) :: a
  INTEGER (kind=int4) :: Nlev
  REAL (kind=real4) :: dp, dt, dz

  Nlev = size(plev,1)
   
  call LOCATE(Zlev,Nlev,z,Ilev)
  Ilev = max(1,min(Nlev-1,Ilev))

  dp = Plev(Ilev+1) - Plev(Ilev)
  dt = Tlev(Ilev+1) - Tlev(Ilev)
  dz = Zlev(Ilev+1) - Zlev(Ilev)

  if (dz /= 0.0) then
    a = (z - Zlev(Ilev))/dz
    p = Plev(Ilev) + a*dp
    t = Tlev(Ilev) + a*dt
  else
    a = 0.0
    p = Plev(Ilev) 
    t = Tlev(Ilev)
  endif
 
END SUBROUTINE prof_lookup_using_z

!----------------------------------------------------------------------
SUBROUTINE COMPUTE_SEGMENT_NWP_CLOUD_PARAMETERS()

  INTEGER:: Elem_Idx 
  INTEGER:: Line_Idx 
  INTEGER:: Lon_Idx 
  INTEGER:: Lat_Idx 

  !--- initialize global arrays which are output of this routine
  Cld_Type_Nwp = Missing_Value_Int1
  Sc_Lwp_Nwp = Missing_Value_Real4
  Lwp_Nwp = Missing_Value_Real4
  Iwp_Nwp = Missing_Value_Real4
  Cwp_Nwp = Missing_Value_Real4
  Pc_Nwp = Missing_Value_Real4
  Ncld_Layers_Nwp = Missing_Value_Int1
  Cloud_Fraction_Satellite_Nwp = Missing_Value_Real4
  High_Cloud_Fraction_Satellite_Nwp = Missing_Value_Real4
  Mid_Cloud_Fraction_Satellite_Nwp = Missing_Value_Real4
  Low_Cloud_Fraction_Satellite_Nwp = Missing_Value_Real4

  !--- loop over pixels in segment
   line_loop: do Line_Idx = 1, Image%Number_Of_Lines_Read_This_Segment
     element_loop: do Elem_Idx = 1, Image%Number_Of_Elements

      !--- check for bad pixels
      if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) then
        cycle
      endif

      !--- check for space views
      if (Space_Mask(Elem_Idx,Line_Idx) == sym%YES) then
        cycle
      endif

      Lon_Idx = I_Nwp(Elem_Idx,Line_Idx)
      Lat_Idx = J_Nwp(Elem_Idx,Line_Idx)

      if (Lon_Idx <= 0 .or. Lat_Idx <= 0) then
        print *, "Missing lat,lon ", Elem_Idx,Line_Idx,Lon_Idx, Lat_Idx
        cycle
      endif

      !--- check to see if this cell was already done
      if (Cld_Type_Nwp(Lon_Idx,Lat_Idx) >= 0) then
        cycle
      endif

      !--- compute cloud parameters from nwp profiles
      call COMPUTE_NWP_CLOUD_PARAMETERS(Tropo_Level_Nwp(Lon_Idx,Lat_Idx), &
                                        Sfc_Level_Nwp(Lon_Idx,Lat_Idx), &
                                        Clwmr_Prof_Nwp(:,Lon_Idx,Lat_Idx), &
                                        Rh_Prof_Nwp(:,Lon_Idx,Lat_Idx), &
                                        T_Prof_Nwp(:,Lon_Idx,Lat_Idx), &
                                        P_Std_Nwp, &
                                        Pc_Nwp(Lon_Idx,Lat_Idx), &
                                        Iwp_Nwp(Lon_Idx,Lat_Idx), &
                                        Lwp_Nwp(Lon_Idx,Lat_Idx), &
                                        Sc_Lwp_Nwp(Lon_Idx,Lat_Idx), &
                                        Cwp_Nwp(Lon_Idx,Lat_Idx), &
                                        Cloud_Fraction_Satellite_Nwp(Lon_Idx,Lat_Idx), &
                                        High_Cloud_Fraction_Satellite_Nwp(Lon_Idx,Lat_Idx), &
                                        Mid_Cloud_Fraction_Satellite_Nwp(Lon_Idx,Lat_Idx), &
                                        Low_Cloud_Fraction_Satellite_Nwp(Lon_Idx,Lat_Idx), &
                                        Ncld_Layers_Nwp(Lon_Idx,Lat_Idx), &
                                        Cld_Type_Nwp(Lon_Idx,Lat_Idx))

     end do element_loop
  end do line_loop

END SUBROUTINE COMPUTE_SEGMENT_NWP_CLOUD_PARAMETERS

!----------------------------------------------------------------------
! subroutine TEMPORAL_INTERP_TMPSFC_NWP(start_itime, end_itime)
!----------------------------------------------------------------------
subroutine TEMPORAL_INTERP_TMPSFC_NWP(start_itime, end_itime)
integer(kind=int4), intent(in):: start_itime
integer(kind=int4), intent(in):: end_itime
real(kind=real4):: x_interp
real(kind=real4):: start_hour
real(kind=real4)::  mean_hours
real(kind=real4):: segment_start_time
real(kind=real4):: segment_end_time
real(kind=real4):: segment_mean_time
integer(kind=int4):: nwp_start_hour_segment
integer(kind=int4):: nwp_end_hour_segment

!--- check for valid times (assume positive)
if (start_itime < 0 .or. end_itime < 0) then
  return
endif

!--- determine the bound day values for this orbit
start_hour = start_itime/60.0/60.0/1000.0

nwp_start_hour_segment = 6*int(start_hour/6)
nwp_end_hour_segment = nwp_start_hour_segment + 6

if (nwp_start_hour == 0 .and. nwp_end_hour_segment == 24) then
   nwp_end_hour_segment = 0
endif
if (nwp_end_hour == 24 .and. nwp_start_hour_segment == 0) then
   nwp_start_hour_segment = 24
endif

!--- initialize interpolation weight
x_interp = Missing_Value_Real4

!--- 
if (nwp_end_hour_segment == nwp_start_hour) then
  x_interp = 0.0
elseif (nwp_start_hour_segment == nwp_end_hour) then
  x_interp = 1.0
endif

if ((nwp_start_hour_segment == nwp_start_hour) .and. & 
    (nwp_end_hour_segment == nwp_end_hour)) then

!----  determine mean time for this segment
segment_start_time = start_itime/86400000.0_real4
segment_end_time = end_itime/86400000.0_real4
segment_mean_time = 0.5*(segment_start_time + segment_end_time)
mean_hours = segment_mean_time*24.0

 !--- determine interpolation weight for 6 hrly ncep fields
 !
 x_interp = (mean_hours  - nwp_start_hour)/  &
            (nwp_end_hour - nwp_start_hour)
 x_interp = max(0.0,min(1.0,x_interp))

endif

 !--- apply interpolation
if (x_interp /= Missing_Value_Real4) then 
 where ( (Tmpsfc_Nwp_Before /= Missing_Value_Real4) .and. (Tmpsfc_Nwp_After /= Missing_Value_Real4) )
    Tmpsfc_Nwp = (1.0 - x_interp) * Tmpsfc_Nwp_Before + x_interp * Tmpsfc_Nwp_After
 end where
endif

end subroutine TEMPORAL_INTERP_TMPSFC_NWP

!--------------------------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------------------------
subroutine CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY_R4(xNwp,xPix,interp_method_in)

    real, dimension(:,:), intent(in):: xNwp
    real, dimension(:,:), intent(out):: xPix
    integer, intent(in), optional:: interp_method_in 
    integer:: interp_method
    integer:: Elem_Idx
    integer:: Line_Idx
    xPix = Missing_Value_Real4
    interp_method = 0
    if (present(interp_method_in)) interp_method = interp_method_in
    do Elem_Idx = 1, size(I_Nwp(:,1))
     do Line_Idx = 1, size(I_Nwp(1,:))
      if ((I_Nwp(Elem_Idx,Line_Idx) > 0).and.(J_Nwp(Elem_Idx,Line_Idx) > 0)) then
        xPix(Elem_Idx,Line_Idx) = xNwp(I_Nwp(Elem_Idx,Line_Idx),J_Nwp(Elem_Idx,Line_Idx))
        if (interp_method /= 0) then
         if ((I_Nwp_x(Elem_Idx,Line_Idx) > 0).and.(J_Nwp_x(Elem_Idx,Line_Idx) > 0)) then
            xPix(Elem_Idx,Line_Idx) = INTERPOLATE_NWP( &
                       xNwp(I_Nwp(Elem_Idx,Line_Idx),J_Nwp(Elem_Idx,Line_Idx)), &
                       xNwp(I_Nwp_x(Elem_Idx,Line_Idx),J_Nwp(Elem_Idx,Line_Idx)), &
                       xNwp(I_Nwp(Elem_Idx,Line_Idx),J_Nwp_x(Elem_Idx,Line_Idx)), &
                       xNwp(I_Nwp_x(Elem_Idx,Line_Idx),J_Nwp_x(Elem_Idx,Line_Idx)), &
                       Lon_Nwp_Fac(Elem_Idx,Line_Idx), &
                       Lat_Nwp_Fac(Elem_Idx,Line_Idx))
         endif
        endif
      endif
     enddo
    enddo
end subroutine CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY_R4
!--------------------------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------------------------
subroutine CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY_I1(xNwp,xPix,interp_method_in)
    integer(kind=int1), dimension(:,:), intent(in):: xNwp
    integer(kind=int1), dimension(:,:), intent(out):: xPix
    integer, intent(in), optional:: interp_method_in 
    integer:: interp_method
    integer:: Elem_Idx
    integer:: Line_Idx
    xPix = Missing_Value_Int1
    interp_method = 0
    if (present(interp_method_in)) interp_method = interp_method_in
    do Elem_Idx = 1, size(I_Nwp(:,1))
     do Line_Idx = 1, size(I_Nwp(1,:))
      if ((I_Nwp(Elem_Idx,Line_Idx) > 0).and.(J_Nwp(Elem_Idx,Line_Idx) > 0)) then
        xPix(Elem_Idx,Line_Idx) = xNwp(I_Nwp(Elem_Idx,Line_Idx),J_Nwp(Elem_Idx,Line_Idx))
        if (interp_method /= 0) then
         if ((I_Nwp_x(Elem_Idx,Line_Idx) > 0).and.(J_Nwp_x(Elem_Idx,Line_Idx) > 0)) then

            xPix(Elem_Idx,Line_Idx) = INTERPOLATE_NWP( &
                       xNwp(I_Nwp(Elem_Idx,Line_Idx),J_Nwp(Elem_Idx,Line_Idx)), &
                       xNwp(I_Nwp_x(Elem_Idx,Line_Idx),J_Nwp(Elem_Idx,Line_Idx)), &
                       xNwp(I_Nwp(Elem_Idx,Line_Idx),J_Nwp_x(Elem_Idx,Line_Idx)), &
                       xNwp(I_Nwp_x(Elem_Idx,Line_Idx),J_Nwp_x(Elem_Idx,Line_Idx)), &
                       Lon_Nwp_Fac(Elem_Idx,Line_Idx), &
                       Lat_Nwp_Fac(Elem_Idx,Line_Idx))
         endif
        endif
      endif
     enddo
    enddo
end subroutine CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY_I1

subroutine CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY_I2(xNwp,xPix,interp_method_in)
    integer(kind=int2), dimension(:,:), intent(in):: xNwp
    integer(kind=int2), dimension(:,:), intent(out):: xPix
    integer, intent(in), optional:: interp_method_in 
    integer:: interp_method
    integer:: Elem_Idx
    integer:: Line_Idx
    xPix = Missing_Value_Int2
    interp_method = 0
    if (present(interp_method_in)) interp_method = interp_method_in
    do Elem_Idx = 1, size(I_Nwp(:,1))
     do Line_Idx = 1, size(I_Nwp(1,:))
      if ((I_Nwp(Elem_Idx,Line_Idx) > 0).and.(J_Nwp(Elem_Idx,Line_Idx) > 0)) then
        xPix(Elem_Idx,Line_Idx) = xNwp(I_Nwp(Elem_Idx,Line_Idx),J_Nwp(Elem_Idx,Line_Idx))
        if (interp_method /= 0) then
         if ((I_Nwp_x(Elem_Idx,Line_Idx) > 0).and.(J_Nwp_x(Elem_Idx,Line_Idx) > 0)) then
            xPix(Elem_Idx,Line_Idx) = INTERPOLATE_NWP( &
                       xNwp(I_Nwp(Elem_Idx,Line_Idx),J_Nwp(Elem_Idx,Line_Idx)), &
                       xNwp(I_Nwp_x(Elem_Idx,Line_Idx),J_Nwp(Elem_Idx,Line_Idx)), &
                       xNwp(I_Nwp(Elem_Idx,Line_Idx),J_Nwp_x(Elem_Idx,Line_Idx)), &
                       xNwp(I_Nwp_x(Elem_Idx,Line_Idx),J_Nwp_x(Elem_Idx,Line_Idx)), &
                       Lon_Nwp_Fac(Elem_Idx,Line_Idx), &
                       Lat_Nwp_Fac(Elem_Idx,Line_Idx))
         endif
        endif
      endif
     enddo
    enddo
end subroutine CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY_I2

subroutine CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY_I4(xNwp,xPix,interp_method_in)
    integer(kind=int4), dimension(:,:), intent(in):: xNwp
    integer(kind=int4), dimension(:,:), intent(out):: xPix
    integer, intent(in), optional:: interp_method_in 
    integer:: interp_method
    integer:: Elem_Idx
    integer:: Line_Idx
    xPix = Missing_Value_Int4
    interp_method = 0
    if (present(interp_method_in)) interp_method = interp_method_in
    do Elem_Idx = 1, size(I_Nwp(:,1))
     do Line_Idx = 1, size(I_Nwp(1,:))
      if ((I_Nwp(Elem_Idx,Line_Idx) > 0).and.(J_Nwp(Elem_Idx,Line_Idx) > 0)) then
        xPix(Elem_Idx,Line_Idx) = xNwp(I_Nwp(Elem_Idx,Line_Idx),J_Nwp(Elem_Idx,Line_Idx))
        if (interp_method /= 0) then
         if ((I_Nwp_x(Elem_Idx,Line_Idx) > 0).and.(J_Nwp_x(Elem_Idx,Line_Idx) > 0)) then
            xPix(Elem_Idx,Line_Idx) = INTERPOLATE_NWP( &
                       xNwp(I_Nwp(Elem_Idx,Line_Idx),J_Nwp(Elem_Idx,Line_Idx)), &
                       xNwp(I_Nwp_x(Elem_Idx,Line_Idx),J_Nwp(Elem_Idx,Line_Idx)), &
                       xNwp(I_Nwp(Elem_Idx,Line_Idx),J_Nwp_x(Elem_Idx,Line_Idx)), &
                       xNwp(I_Nwp_x(Elem_Idx,Line_Idx),J_Nwp_x(Elem_Idx,Line_Idx)), &
                       Lon_Nwp_Fac(Elem_Idx,Line_Idx), &
                       Lat_Nwp_Fac(Elem_Idx,Line_Idx))
         endif
        endif
      endif
     enddo
    enddo
end subroutine CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY_I4


end module NWP_COMMON
