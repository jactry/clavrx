! $Id$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: sfc_prop_umd.f90 (src)
!       SURFACE_PROPERTIES (program)
!
! PURPOSE: 
!          this is module of surface properties (reflectance, emissivity) for
!          the UMD surface type classification
!
! DESCRIPTION: 
!     0: Water
!     1: Evergreen Needleleaf Forests
!     2: Evergreen Broadleaf Forests
!     3: Deciduous Needleleaf Forests
!     4: Deciduous Broadleaf Forests
!     5: Mixed Forests
!     6: Woodlands
!     7: Wooded Grasslands/Shrubs
!     8: Closed Bushlands or Shrublands
!     9: Open Shrublands
!     10: Grasses
!     11: Cropland
!     12: Bare
!     13: Urban and Built
!
!     Note, surface albedoes range from 0 to 1.0
!
! AUTHORS:
!  Andrew Heidinger, Andrew.Heidinger@noaa.gov
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
! File I/O: None
!
! Public Routines
!  SETUP_UMD_PROPS - assign values to radiative properties for each land type
!
! Private Routines: None
!
!--------------------------------------------------------------------------------------
module SURFACE_PROPERTIES
use CONSTANTS
use PIXEL_COMMON

implicit none

private
public:: SETUP_UMD_PROPS,  &
         GET_PIXEL_SFC_EMISS_FROM_SFC_TYPE, &
         COMPUTE_BINARY_LAND_COAST_MASKS, &
         COMPUTE_COAST_MASK_FROM_LAND_MASK

integer, parameter, public:: ntype_sfc=14
real(kind=real4), dimension(0:ntype_sfc-1), public,save:: Ch1_Sfc_Alb_Umd
real(kind=real4), dimension(0:ntype_sfc-1), public,save:: Ch2_Sfc_Alb_Umd
real(kind=real4), dimension(0:ntype_sfc-1), public,save:: Ch6_Sfc_Alb_Umd
real(kind=real4), dimension(0:ntype_sfc-1), public,save:: Ch20_Sfc_Alb_Umd
real(kind=real4), dimension(0:ntype_sfc-1), public,save:: Ch20_Sfc_Emiss_Umd
real(kind=real4), dimension(0:ntype_sfc-1), public,save:: Ch27_Sfc_Emiss_Umd
real(kind=real4), dimension(0:ntype_sfc-1), public,save:: Ch28_Sfc_Emiss_Umd
real(kind=real4), dimension(0:ntype_sfc-1), public,save:: Ch29_Sfc_Emiss_Umd
real(kind=real4), dimension(0:ntype_sfc-1), public,save:: Ch31_Sfc_Emiss_Umd
real(kind=real4), dimension(0:ntype_sfc-1), public,save:: Ch32_Sfc_Emiss_Umd
real(kind=real4), dimension(0:ntype_sfc-1), public,save:: Ch33_Sfc_Emiss_Umd
real(kind=real4), dimension(0:ntype_sfc-1), public,save:: Ch1_Snow_Sfc_Alb_Umd
real(kind=real4), dimension(0:ntype_sfc-1), public,save:: Ch2_Snow_Sfc_Alb_Umd
real(kind=real4), dimension(0:ntype_sfc-1), public,save:: Ch5_Snow_Sfc_Alb_Umd
real(kind=real4), dimension(0:ntype_sfc-1), public,save:: Ch6_Snow_Sfc_Alb_Umd
real(kind=real4), dimension(0:ntype_sfc-1), public,save:: Ch7_Snow_Sfc_Alb_Umd

contains

!-----------------------------------------------------------------
! read into memory the radiative properties for each UMD type
!-----------------------------------------------------------------
subroutine SETUP_UMD_PROPS()

!--- channel 1 surface albedo based on PATMOS-x
!  Ch1_Sfc_Alb_Umd = (/0.0300, & !0

   Ch1_Sfc_Alb_Umd = (/0.0500, & !0
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
     Ch1_Sfc_Alb_Umd = 100.0*Ch1_Sfc_Alb_Umd

     !--- channel 2 surface albedo based on PATMOS-x
     Ch2_Sfc_Alb_Umd = (/0.0163,  &  !0
                         0.1991,  &  !1
                         0.2179,  &  !2
                         0.1991,  &  !3
                         0.2943,  &  !4
                         0.2500,  &  !5
                         0.2500,  &  !6
                         0.2700,  &  !7
                         0.3072,  &  !8
                         0.2942,  &  !9
                         0.2919,  &  !10
                         0.2520,  &  !11
                         0.4019,  &  !12
                         0.2466/)   !13
     Ch2_Sfc_Alb_Umd = 100.0*Ch2_Sfc_Alb_Umd

     !--- channel 3a surface albedo based on PATMOS-x
     Ch6_Sfc_Alb_Umd = (/0.0163, & !0
                         0.2020, & !1
                         0.2210, & !2
                         0.2020, & !3
                         0.2980, & !4
                         0.2500, & !5
                         0.2500, & !6
                         0.2700, & !7
                         0.3400, & !8
                         0.2580, & !9
                         0.2980, & !10
                         0.2520, & !11
                         0.3900, & !12
                         0.2650/)  !13
    Ch6_Sfc_Alb_Umd = 100.0*Ch6_Sfc_Alb_Umd

    !--- channel 3b surface emissivity based on SeeBor et al 2006
    Ch20_Sfc_Emiss_Umd =(/0.985, &   !0
                          0.95,  &   !1
                          0.96,  &   !2
                          0.94,  &   !3
                          0.94,  &   !4
                          0.94,  &   !5
                          0.94,  &   !6
                          0.93,  &   !7
                          0.89,  &   !8
                          0.89,  &   !9
                          0.92,  &   !10
                          0.94,  &   !11
                          0.86,  &   !12
                          0.94/)     !13

    !--- channel 3b surface albedo based on SeeBor et al 2006
    Ch20_Sfc_Alb_Umd = 100.0*(1.0 - Ch20_Sfc_Emiss_Umd)

    !--- channel 28 surface emissivity based on SeeBor et al 2006 -- FAKED
    Ch27_Sfc_Emiss_Umd = (/0.985,  &   !0
                     0.965,  &   !1
                     0.961,  &   !2
                     0.963,  &   !3
                     0.963,  &   !4
                     0.966,  &   !5
                     0.962,  &   !6
                     0.961,  &   !7
                     0.961,  &   !8
                     0.958,  &   !9
                     0.969,  &   !10
                     0.965,  &   !11
                     0.955,  &   !12
                     0.960/)     !13

   !--- channel 28 surface emissivity based on SeeBor et al 2006 -- FAKED
   Ch28_Sfc_Emiss_Umd = (/0.985,  &   !0
                     0.965,  &   !1
                     0.961,  &   !2
                     0.963,  &   !3
                     0.963,  &   !4
                     0.966,  &   !5
                     0.962,  &   !6
                     0.961,  &   !7
                     0.961,  &   !8
                     0.958,  &   !9
                     0.969,  &   !10
                     0.965,  &   !11
                     0.955,  &   !12
                     0.960/)     !13

   !--- channel 29 surface emissivity based on SeeBor et al 2006 -- FAKED
   Ch29_Sfc_Emiss_Umd = (/0.985,  &   !0
                     0.965,  &   !1
                     0.961,  &   !2
                     0.963,  &   !3
                     0.963,  &   !4
                     0.966,  &   !5
                     0.962,  &   !6
                     0.961,  &   !7
                     0.961,  &   !8
                     0.958,  &   !9
                     0.969,  &   !10
                     0.965,  &   !11
                     0.955,  &   !12
                     0.960/)     !13

   !--- channel 31 surface emissivity based on SeeBor et al 2006
   Ch31_Sfc_Emiss_Umd = (/0.985,  &   !0
                     0.965,  &   !1
                     0.961,  &   !2
                     0.963,  &   !3
                     0.963,  &   !4
                     0.966,  &   !5
                     0.962,  &   !6
                     0.961,  &   !7
                     0.961,  &   !8
                     0.958,  &   !9
                     0.969,  &   !10
                     0.965,  &   !11
                     0.955,  &   !12
                     0.960/)     !13

   !--- channel 32 surface emissivity based on SeeBor et al 2006
   Ch32_Sfc_Emiss_Umd = (/0.985,  &   !0
                     0.967,  &   !1
                     0.965,  &   !2
                     0.965,  &   !3
                     0.966,  &   !4
                     0.968,  &   !5
                     0.965,  &   !6
                     0.966,  &   !7
                     0.969,  &   !8
                     0.968,  &   !9
                     0.964,  &   !10
                     0.969,  &   !11
                     0.971,  &   !12
                     0.962/)     !13

   !--- channel 33 surface emissivity based on SeeBor et al 2006 - FAKED
   Ch33_Sfc_Emiss_Umd = (/0.985,  &   !0
                     0.967,  &   !1
                     0.965,  &   !2
                     0.965,  &   !3
                     0.966,  &   !4
                     0.968,  &   !5
                     0.965,  &   !6
                     0.966,  &   !7
                     0.969,  &   !8
                     0.968,  &   !9
                     0.964,  &   !10
                     0.969,  &   !11
                     0.971,  &   !12
                     0.962/)     !13
   
!------------------------------------------------------------------------------
! Snow reflectance values computed as the mode of distributions from 
! MODIS/TERRA for each of the surface types
! This is assuming 0=ocean 1=sea ice, 2-14=Landsurf_type-1, 15 - not used
!------------------------------------------------------------------------------
   Ch1_Snow_Sfc_Alb_Umd = (/0.66,  &  !0
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
   Ch1_Snow_Sfc_Alb_Umd = 100.0*Ch1_Snow_Sfc_Alb_Umd

   Ch2_Snow_Sfc_Alb_Umd = (/0.61,  &  !0
                            0.29 , &  !1   
                            0.60,  &  !2   
                            0.43,  &  !3  
                            0.34,  &  !4  
                            0.31,  &  !5   
                            0.52,  &  !6  
                            0.56,  &  !7  
                            0.60,  &  !8  
                            0.76,  &  !9  
                            0.69,  &  !10 
                            0.71,  &  !11 
                            0.72,  &  !12 
                            0.65/)    !13 
   Ch2_Snow_Sfc_Alb_Umd = 100.0*Ch2_Snow_Sfc_Alb_Umd

   Ch5_Snow_Sfc_Alb_Umd = (/0.12,  &  !0
                            0.18 , &  !1  
                            0.70,  &  !2  
                            0.22,  &  !3  
                            0.27,  &  !4  
                            0.20,  &  !5   
                            0.18,  &  !6  
                            0.18,  &  !7  
                            0.42,  &  !8  
                            0.45,  &  !9  
                            0.41,  &  !10 
                            0.41,  &  !11 
                            0.43,  &  !12 
                            0.52/)    !13 
   Ch5_Snow_Sfc_Alb_Umd = 100.0*Ch5_Snow_Sfc_Alb_Umd


   Ch6_Snow_Sfc_Alb_Umd = (/0.02,  &  !0
                            0.08 , &  !1  
                            0.64,  &  !2   
                            0.15,  &  !3   
                            0.13,  &  !4   
                            0.09,  &  !5   
                            0.09,  &  !6   
                            0.09,  &  !7   
                            0.09,  &  !8   
                            0.10,  &  !9    
                            0.09,  &  !10   
                            0.12,  &  !11   
                            0.09,  &  !12   
                            0.14/)    !13   
   Ch6_Snow_Sfc_Alb_Umd = 100.0*Ch6_Snow_Sfc_Alb_Umd

   Ch7_Snow_Sfc_Alb_Umd = (/0.13,  &  !0
                            0.03 , &  !1   !TBD 
                            0.13,  &  !2    !TBD 
                            0.08,  &  !3    !TBD
                            0.08,  &  !4    !TBD
                            0.08,  &  !5    !TBD
                            0.08,  &  !6    !TBD
                            0.08,  &  !7    !TBD
                            0.08,  &  !8    !TBD
                            0.08,  &  !9    !TBD
                            0.08,  &  !10   !TBD
                            0.08,  &  !11   !TBD
                            0.03,  &  !12   !TBD
                            0.08/)    !13   !TBD
   Ch7_Snow_Sfc_Alb_Umd = 100.0*Ch7_Snow_Sfc_Alb_Umd

end subroutine SETUP_UMD_PROPS

!-----------------------------------------------------------------------
! If no sfc emiss data base read in, base it on the surface type.
! also derive surface reflectance from surface emissivity
! also reset sfc emiss if snow covered
!-----------------------------------------------------------------------
 subroutine GET_PIXEL_SFC_EMISS_FROM_SFC_TYPE(j1,j2)
    integer, intent(in):: j1, j2
    integer:: i, j
    

j_loop:  do j = j1,j1+j2-1

i_loop:    do i = 1,Image%Number_Of_Elements

!--- check for a bad pixel
         if (Bad_Pixel_Mask(i,j) == sym%YES) then
            cycle
         endif
 
!--- based on surface type, assign surface emissivities if seebor not used
!--- if the Sfc_Type is missing, treat pixel as bad
    if (use_seebor == sym%NO) then
      if (Sfc%Sfc_Type(i,j) >= 0 .and. Sfc%Sfc_Type(i,j) < 15) then     !limits depend on Sfc_Type
       if (Sensor%Chan_On_Flag_Default(20)==sym%YES) ch(20)%Sfc_Emiss(i,j) = Ch20_Sfc_Emiss_Umd(Sfc%Sfc_Type(i,j))
       if (Sensor%Chan_On_Flag_Default(21)==sym%YES) ch(21)%Sfc_Emiss(i,j) = Ch20_Sfc_Emiss_Umd(Sfc%Sfc_Type(i,j))
       if (Sensor%Chan_On_Flag_Default(22)==sym%YES) ch(22)%Sfc_Emiss(i,j) = Ch20_Sfc_Emiss_Umd(Sfc%Sfc_Type(i,j))
       if (Sensor%Chan_On_Flag_Default(23)==sym%YES) ch(23)%Sfc_Emiss(i,j) = Ch20_Sfc_Emiss_Umd(Sfc%Sfc_Type(i,j))
       if (Sensor%Chan_On_Flag_Default(24)==sym%YES) ch(24)%Sfc_Emiss(i,j) = Ch20_Sfc_Emiss_Umd(Sfc%Sfc_Type(i,j))
       if (Sensor%Chan_On_Flag_Default(25)==sym%YES) ch(25)%Sfc_Emiss(i,j) = Ch20_Sfc_Emiss_Umd(Sfc%Sfc_Type(i,j))
       if (Sensor%Chan_On_Flag_Default(27)==sym%YES) ch(27)%Sfc_Emiss(i,j) = Ch27_Sfc_Emiss_Umd(Sfc%Sfc_Type(i,j))
       if (Sensor%Chan_On_Flag_Default(28)==sym%YES) ch(28)%Sfc_Emiss(i,j) = Ch28_Sfc_Emiss_Umd(Sfc%Sfc_Type(i,j))
       if (Sensor%Chan_On_Flag_Default(29)==sym%YES) ch(29)%Sfc_Emiss(i,j) = Ch29_Sfc_Emiss_Umd(Sfc%Sfc_Type(i,j))
       if (Sensor%Chan_On_Flag_Default(31)==sym%YES) ch(31)%Sfc_Emiss(i,j) = Ch31_Sfc_Emiss_Umd(Sfc%Sfc_Type(i,j))
       if (Sensor%Chan_On_Flag_Default(32)==sym%YES) ch(32)%Sfc_Emiss(i,j) = Ch32_Sfc_Emiss_Umd(Sfc%Sfc_Type(i,j))
       if (Sensor%Chan_On_Flag_Default(33)==sym%YES) ch(33)%Sfc_Emiss(i,j) = Ch33_Sfc_Emiss_Umd(Sfc%Sfc_Type(i,j))
       if (Sensor%Chan_On_Flag_Default(34)==sym%YES) ch(34)%Sfc_Emiss(i,j) = Ch33_Sfc_Emiss_Umd(Sfc%Sfc_Type(i,j))
       if (Sensor%Chan_On_Flag_Default(35)==sym%YES) ch(35)%Sfc_Emiss(i,j) = Ch33_Sfc_Emiss_Umd(Sfc%Sfc_Type(i,j))
       if (Sensor%Chan_On_Flag_Default(36)==sym%YES) ch(36)%Sfc_Emiss(i,j) = Ch33_Sfc_Emiss_Umd(Sfc%Sfc_Type(i,j))
      else
       if (Sensor%Chan_On_Flag_Default(20)==sym%YES) ch(20)%Sfc_Emiss(i,j) = Missing_Value_Real4
       if (Sensor%Chan_On_Flag_Default(21)==sym%YES) ch(21)%Sfc_Emiss(i,j) = Missing_Value_Real4
       if (Sensor%Chan_On_Flag_Default(27)==sym%YES) ch(27)%Sfc_Emiss(i,j) = Missing_Value_Real4
       if (Sensor%Chan_On_Flag_Default(28)==sym%YES) ch(28)%Sfc_Emiss(i,j) = Missing_Value_Real4
       if (Sensor%Chan_On_Flag_Default(29)==sym%YES) ch(29)%Sfc_Emiss(i,j) = Missing_Value_Real4
       if (Sensor%Chan_On_Flag_Default(31)==sym%YES) ch(31)%Sfc_Emiss(i,j) = Missing_Value_Real4
       if (Sensor%Chan_On_Flag_Default(32)==sym%YES) ch(32)%Sfc_Emiss(i,j) = Missing_Value_Real4
       if (Sensor%Chan_On_Flag_Default(33)==sym%YES) ch(33)%Sfc_Emiss(i,j) = Missing_Value_Real4
       if (Sensor%Chan_On_Flag_Default(34)==sym%YES) ch(34)%Sfc_Emiss(i,j) = Missing_Value_Real4
       if (Sensor%Chan_On_Flag_Default(35)==sym%YES) ch(35)%Sfc_Emiss(i,j) = Missing_Value_Real4
       if (Sensor%Chan_On_Flag_Default(36)==sym%YES) ch(36)%Sfc_Emiss(i,j) = Missing_Value_Real4
       Bad_Pixel_Mask(i,j) = sym%YES
      endif
    endif

!--- for ocean, use this parameterization from Nick Nalli - overwrite seebor
    if (Sfc%Sfc_Type(i,j) == 0) then
!    Sfc_Emiss_Ch31(i,j) =  0.844780 + 0.328921 * coszen(i,j)  -0.182375*(coszen(i,j)**2)
!    Sfc_Emiss_Ch32(i,j) =  0.775019 + 0.474005 * coszen(i,j)  -0.261739*(coszen(i,j)**2)
!--- set to constant until we can verify the above
     if (Sensor%Chan_On_Flag_Default(20)==sym%YES) ch(20)%Sfc_Emiss(i,j) = 0.985
     if (Sensor%Chan_On_Flag_Default(21)==sym%YES) ch(21)%Sfc_Emiss(i,j) = 0.985
     if (Sensor%Chan_On_Flag_Default(22)==sym%YES) ch(22)%Sfc_Emiss(i,j) = 0.985
     if (Sensor%Chan_On_Flag_Default(23)==sym%YES) ch(23)%Sfc_Emiss(i,j) = 0.985
     if (Sensor%Chan_On_Flag_Default(24)==sym%YES) ch(24)%Sfc_Emiss(i,j) = 0.985
     if (Sensor%Chan_On_Flag_Default(25)==sym%YES) ch(25)%Sfc_Emiss(i,j) = 0.985
     if (Sensor%Chan_On_Flag_Default(27)==sym%YES) ch(27)%Sfc_Emiss(i,j) = 0.985
     if (Sensor%Chan_On_Flag_Default(28)==sym%YES) ch(28)%Sfc_Emiss(i,j) = 0.985
     if (Sensor%Chan_On_Flag_Default(29)==sym%YES) ch(29)%Sfc_Emiss(i,j) = 0.985
     if (Sensor%Chan_On_Flag_Default(30)==sym%YES) ch(30)%Sfc_Emiss(i,j) = 0.985
     if (Sensor%Chan_On_Flag_Default(31)==sym%YES) ch(31)%Sfc_Emiss(i,j) = 0.985
     if (Sensor%Chan_On_Flag_Default(32)==sym%YES) ch(32)%Sfc_Emiss(i,j) = 0.985
     if (Sensor%Chan_On_Flag_Default(33)==sym%YES) ch(33)%Sfc_Emiss(i,j) = 0.985
     if (Sensor%Chan_On_Flag_Default(34)==sym%YES) ch(34)%Sfc_Emiss(i,j) = 0.985
     if (Sensor%Chan_On_Flag_Default(35)==sym%YES) ch(35)%Sfc_Emiss(i,j) = 0.985
     if (Sensor%Chan_On_Flag_Default(36)==sym%YES) ch(36)%Sfc_Emiss(i,j) = 0.985
    endif

    !--- if (snow_mask
    if (Sfc%Snow(i,j) /= sym%NO_SNOW) then

      if ((Sfc%Sfc_Type(i,j) /= sym%EVERGREEN_NEEDLE_SFC) .and. &
          (Sfc%Sfc_Type(i,j) /= sym%EVERGREEN_BROAD_SFC)) then

      if (Sensor%Chan_On_Flag_Default(20)==sym%YES) ch(20)%Sfc_Emiss(i,j) = 0.984
      if (Sensor%Chan_On_Flag_Default(27)==sym%YES) ch(27)%Sfc_Emiss(i,j) = 0.979
      if (Sensor%Chan_On_Flag_Default(28)==sym%YES) ch(28)%Sfc_Emiss(i,j) = 0.979
      if (Sensor%Chan_On_Flag_Default(29)==sym%YES) ch(29)%Sfc_Emiss(i,j) = 0.979
      if (Sensor%Chan_On_Flag_Default(31)==sym%YES) ch(31)%Sfc_Emiss(i,j) = 0.979
      if (Sensor%Chan_On_Flag_Default(32)==sym%YES) ch(32)%Sfc_Emiss(i,j) = 0.977
      if (Sensor%Chan_On_Flag_Default(33)==sym%YES) ch(33)%Sfc_Emiss(i,j) = 0.977
      if (Sensor%Chan_On_Flag_Default(34)==sym%YES) ch(34)%Sfc_Emiss(i,j) = 0.977
      if (Sensor%Chan_On_Flag_Default(35)==sym%YES) ch(35)%Sfc_Emiss(i,j) = 0.977
      if (Sensor%Chan_On_Flag_Default(36)==sym%YES) ch(36)%Sfc_Emiss(i,j) = 0.977

     endif
    endif

  end do i_loop
 end do j_loop

 end subroutine GET_PIXEL_SFC_EMISS_FROM_SFC_TYPE
!----------------------------------------------------------------------
! based on the coast and land flags, derive binary land and coast
! masks (yes/no)
!
! Note, coast mask is dependent on sensor resolution
!----------------------------------------------------------------------
 subroutine COMPUTE_BINARY_LAND_COAST_MASKS(j1,j2)
    integer, intent(in):: j1, j2
    integer:: i, j
    
j_loop:  do j = j1,j1+j2-1

i_loop:    do i = 1,Image%Number_Of_Elements

   !--- check for a bad pixel
   if (Bad_Pixel_Mask(i,j) == sym%YES) then
      cycle
   endif

   !--- binary land mask
   Sfc%Land_Mask(i,j) = sym%NO

   !--- if land mask read in, use it
   if (Read_Land_Mask == sym%YES) then
     if (Sfc%Land(i,j) == sym%LAND) then
       Sfc%Land_Mask(i,j) = sym%YES
     endif
   !--- if land mask not read in, base off of surface type
   else   
     if (Sfc%Sfc_Type(i,j) /= sym%WATER_SFC) then
       Sfc%Land_Mask(i,j) = sym%YES
     endif
   endif

   !--- binary coast mask
   if (Read_Coast_Mask == sym%YES) then
    Sfc%Coast_Mask(i,j) = sym%NO

     !-- for gac data
     if ((Sensor%Spatial_Resolution_Meters <= 1000) .and. (Sfc%Coast(i,j) /= sym%NO_COAST)) then
      if (Sfc%Coast(i,j) <= sym%COAST_10KM) then
         Sfc%Coast_Mask(i,j) = sym%YES
      endif
     endif

     !-- for lac,hrpt or frac data
     if ((Sensor%Spatial_Resolution_Meters > 1000) .and. (Sfc%Coast(i,j) /= sym%NO_COAST)) then
      if (Sfc%Coast(i,j) <= sym%COAST_5KM) then
         Sfc%Coast_Mask(i,j) = sym%YES
      endif
     endif

   endif

  end do i_loop
 end do j_loop



!------------- compute coast mask if coast data not read in
 if (Read_Coast_Mask == sym%NO) then 

   call COMPUTE_COAST_MASK_FROM_LAND_MASK(j1,j2)

 endif

 end subroutine COMPUTE_BINARY_LAND_COAST_MASKS


!------------------------------------------------------------------------
! compute the data_mask array - this holds information pertaining to the
! type of pixel this is and how the cloud mask treats it
!------------------------------------------------------------------------------
 subroutine COMPUTE_COAST_MASK_FROM_LAND_MASK(jmin,jmax)
  integer, intent(in):: jmin,jmax
  integer:: i, j, i1, i2, j1,j2,nland,nmax,n

  n = 2    !determine size of box, 2 = 5x5
  do j = jmin,jmax

    j1 = max(jmin,j-n)
    j2 = min(j+n,jmax)

    do i = 1, Image%Number_Of_Elements

     i1 = max(1,i-n)
     i2 = min(i+n,int(Image%Number_Of_Elements, kind=int4))

     !--- check for bad pixels
     if (Bad_Pixel_Mask(i,j) /= 0) then
       cycle
     endif

     !--- compute number of land pixels in box
     Nland = sum(Sfc%Land_Mask(i1:i2,j1:j2))

     !--- compute maximum number of land pixels in box
     Nmax = (j2-j1+1)*(i2-i1+1)

     !--- compute coast mask
     Sfc%Coast_Mask(i,j) = sym%YES
     if ((Nland == Nmax).or.(Nland == 0)) then
      Sfc%Coast_Mask(i,j) = sym%NO
     endif

   enddo
  enddo

end subroutine COMPUTE_COAST_MASK_FROM_LAND_MASK

end module SURFACE_PROPERTIES
