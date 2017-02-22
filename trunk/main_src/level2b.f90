!$Id$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: level2b.f90 (src)
!       LEVEL2B_ROUTINES (program)
!
! PURPOSE: Routines for creating, writing and closing  L2b output files
!
! DESCRIPTION: 
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
! REVISON HISTORY:
!    Creation Date May 2009
!--------------------------------------------------------------------------------------
module LEVEL2B_ROUTINES
 use CONSTANTS
 use HDF
 use LEVEL2_ROUTINES
 use SCALING_PARAMETERS
 use FILE_UTILITY

 implicit none

 public:: DEFINE_SDS_RANK1
 public:: DEFINE_SDS_RANK2
 public:: DEFINE_SDS_RANK3
 public:: READ_SDS     
 public:: UNSCALE_SDS  
 public:: SCALE_SDS     
 public:: WRITE_SDS     
 public:: COPY_GLOBAL_ATTRIBUTES    
 public:: REGRID
 public:: SUBSET_LEVEL2B
 public:: INIT_RANDOM_SEED
 public:: FILE_SEARCH
 public:: COMPUTE_WMO_ID_KNOWING_SENSOR_NAME

 private:: WRITE_FLAG_ATTRIBUTES
 private:: WRITE_SCALING_ATTRIBUTES
 private:: READ_SCALING_ATTRIBUTES

 interface READ_SDS
     module procedure  &
         READ_SDS_RANK1,  &
         READ_SDS_RANK2,  &
         READ_SDS_RANK3
 end interface

 interface WRITE_SDS
     module procedure      &
         WRITE_SDS_RANK1,  &
         WRITE_SDS_RANK2,  &
         WRITE_SDS_RANK3
 end interface

 interface SCALE_SDS
     module procedure      &
         SCALE_SDS_RANK1,  &
         SCALE_SDS_RANK2,  &
         SCALE_SDS_RANK3
 end interface

 interface UNSCALE_SDS
     module procedure      &
         UNSCALE_SDS_RANK2,  &
         UNSCALE_SDS_RANK3
 end interface

!interface DEFINE_SDS
!    module procedure  &
!        DEFINE_SDS_RANK1,  &
!        DEFINE_SDS_RANK2,  &
!        DEFINE_SDS_RANK3
!end interface

 type, public :: Sds_Struct
  character(len=100):: Variable_Name
  integer(kind=int4):: Data_Type
  integer(kind=int4):: Rank
  integer(kind=int1):: Scaling_Type
  character(len=100):: Units
  character(len=300):: Long_Name
  character(len=100):: Standard_Name
  real(kind=real4):: Unscaled_Min
  real(kind=real4):: Unscaled_Max
  real(kind=real4):: Unscaled_Missing
  real(kind=real4):: Scale_Factor
  real(kind=real4):: Add_Offset
  real(kind=real4), dimension(2):: Actual_Range
  real(kind=real4):: Fill_Value
  integer(kind=int4), dimension(2):: Valid_Range
  integer(kind=int4):: Status_Flag
  integer(kind=int4):: Id_Input
  integer(kind=int4):: Id_Output
  integer(kind=int4):: Compression_Type
  integer(kind=int4):: Num_Attrs
 end type Sds_Struct

 character(len=30), parameter, private :: MOD_PROMPT = " PATMOSx_LEVEL2B_ROUTINES: "
 character(len=18), private, parameter:: Coordinates_String = "longitude latitude"

 real(kind=real4), dimension(:,:), allocatable, public, save:: Lat_Input
 real(kind=real4), dimension(:,:), allocatable, public, save:: Lon_Input
 real(kind=real4), dimension(:,:), allocatable, public, save:: Gap_Pixel_Mask_Input
 real(kind=real4), dimension(:,:), allocatable, public, save:: Lat_Output
 real(kind=real4), dimension(:,:), allocatable, public, save:: Lon_Output
 integer(kind=int4), dimension(:,:), allocatable, public, save:: Ielem_Output
 integer(kind=int4), dimension(:,:), allocatable, public, save:: Iline_Output

 !--- no compression
 !integer, parameter, private:: Comp_Type = 0
 !integer, dimension(2), parameter, private:: Comp_Prm = /(0,0)/
 !--- gzip compression
 integer, parameter, private:: Comp_Type = 4
 integer, dimension(4), parameter, private:: Comp_Prm = (/6,0,0,0/)
 !--- szip compression
 !integer, parameter, private:: Comp_Type = 5
 !integer, dimension(2), parameter, private:: Comp_Prm = /(32,0)/


 contains

!====================================================================
! subroutine Name: REGRID
!
! Function: 
!       routines to remap data to a regular lat,lon grid
!
! Description:
!   This subroutine is a routine to remap data from lat/lon to a 
!   static grid.
!
! Inputs:
!       Lon_Output = vector of the longitudes of the grid
!       Lat_Output = vector of the latitudes of the grid
!       Ielem_Output = array of element indices for each grid point
!       Iline_Output = array of line indices for each grid point
!       Random_Flag = set 1 to select a random point in the grid, set
!                     set 0 to select nearest point to grid lattice
!====================================================================
subroutine REGRID( &
        Lat_South, &
        dlat,      &
        Nlat,      &
        Lon_West,  &
        dlon,      &
        Nlon,      &
        Random_Flag)

real(kind=real4), intent(in):: Lat_South
real(kind=real4), intent(in):: dlat
integer(kind=int4), intent(in):: Nlat
real(kind=real4), intent(in):: Lon_West
real(kind=real4), intent(in):: dlon
integer(kind=int4), intent(in):: Nlon
integer(kind=int4), intent(in):: Random_Flag

real(kind=real4):: big_number
real(kind=real4):: Max_Allowable_Distance
real(kind=real4), dimension(:,:),allocatable:: Min_Dist_Grid
integer:: Ielem
integer:: Iline
integer:: Ilon
integer:: Ilat
integer:: nx
integer:: ny
real:: Lat_North
real:: Lon_East
real:: xdist
real:: xdist_01
real:: xdist_02
real:: xdist_10
real:: xdist_20
real:: lat_width
real:: lon_width
integer:: Ilon_1, Ilon_2
integer:: Ilon_min, Ilon_max
real:: lon_1, lon_2
integer:: Ilat_1, Ilat_2
integer:: Ilat_min, Ilat_max
real:: lat_1, lat_2
integer:: i1,i2,j1,j2
integer:: dateline_flag
real:: Lon_East_abs
real:: Lon_Input_abs
real:: xrand

big_number = 999999.9

!--- construct boundaries
Lon_East_abs = Lon_West + dlon * (Nlon-1)
Lon_East = Lon_East_abs
if (Lon_East > 180.0) Lon_East = Lon_East - 360.0

dateline_flag = sym%NO
if (Lon_West > 0 .and. Lon_East < 0.0) dateline_flag = sym%YES

Lat_North = Lat_South + dlat * (Nlat-1)

nx = size(Lat_Input,1)
ny = size(Lat_Input,2)

Ielem_Output = Missing_Value_Int4
Iline_Output = Missing_Value_Int4

Max_Allowable_Distance = 1.0*sqrt(dlat**2 + dlon**2)

allocate(Min_Dist_Grid(Nlon,Nlat))
Min_Dist_Grid = big_number

do Ielem = 1, nx
  do Iline = 1, ny
       
    if ((Lat_Input(Ielem,Iline) >= Lat_South) .and. (Lat_Input(Ielem,Iline) <= Lat_North)) then

      Lon_Input_abs = Lon_Input(Ielem,Iline)
      if (Lon_Input_abs < 0.0) Lon_Input_abs = Lon_Input_abs + 360.0

      if ( ( (dateline_flag == sym%NO) .and. &                !dateline check
             (Lon_Input(Ielem,Iline) >= Lon_West) .and. &
             (Lon_Input(Ielem,Iline) <= Lon_East)) .or. & 
           ( (dateline_flag == sym%YES) .and.  &
             (Lon_Input_abs >= Lon_West) .and. &
             (Lon_Input_abs <= Lon_East_abs))) then 

        i1 = max(1,Ielem-1)
        i2 = min(nx,Ielem+1)
        j1 = max(1,Iline-1)
        j2 = min(ny,Iline+1)

        lon_width = abs(Lon_Input(i1,Iline) - Lon_Input(i2,Iline))
        lon_1 = Lon_Input(Ielem,Iline) + lon_width/2.0
        lon_2 = Lon_Input(Ielem,Iline) - lon_width/2.0

        !-- lon width constraints prevents large values in presence of scanline jumps
        lat_width = min(lon_width,abs(Lat_Input(Ielem,j1) - Lat_Input(Ielem,j2)))
        lat_1 = Lat_Input(Ielem,Iline) + lat_width/2.0
        lat_2 = Lat_Input(Ielem,Iline) - lat_width/2.0

        call POSITION_TO_INDICES (Ilat, Ilon, Lat_Input(Ielem,Iline), Lon_Input(Ielem,Iline),  &
            Lat_South, Lat_North, dlat, Nlat, Lon_West, Lon_East,dlon, Nlon, xdist)

        call POSITION_TO_INDICES (Ilat, Ilon_1, Lat_Input(Ielem,Iline), lon_1,  &
            Lat_South, Lat_North, dlat, Nlat, Lon_West, Lon_East,dlon, Nlon, xdist_10)

        call POSITION_TO_INDICES (Ilat, Ilon_2, Lat_Input(Ielem,Iline), lon_2,  &
            Lat_South, Lat_North, dlat, Nlat, Lon_West, Lon_East,dlon, Nlon, xdist_20)

        call POSITION_TO_INDICES (Ilat_1, Ilon, lat_1, Lon_Input(Ielem,Iline),  &
            Lat_South, Lat_North, dlat, Nlat, Lon_West, Lon_East,dlon, Nlon, xdist_01)

        call POSITION_TO_INDICES (Ilat_2, Ilon, lat_2, Lon_Input(Ielem,Iline),  &
            Lat_South, Lat_North, dlat, Nlat, Lon_West, Lon_East,dlon, Nlon, xdist_02)

        Ilon_min = min(Ilon_1,Ilon_2)
        Ilon_max = max(Ilon_1,Ilon_2)
        Ilat_min = min(Ilat_1,Ilat_2)
        Ilat_max = max(Ilat_1,Ilat_2)

        !update closest grid point

        if ((Ilon >= 1) .and. (Ilat >= 1)) then    !check for valid regrid result

!         Lat_South_This_Gridbox = Lat_South + (Ilat-1)*Dlat
!         Lat_North_This_Gridbox = Lat_South_This_Grid_Box + Dlat
!         Lon_West_This_Gridbox = Lon_West + (Ilon-1)*Dlon
!         Lon_East_This_Gridbox = Lon_West + Dlon

          !--- make sure this pixel falls within the particular gridbox
!         if ( (Lat_Input(Ielem,Iline) >= Lat_South_This_Gridbox) .and. &
!              (Lat_Input(Ielem,Iline) < Lat_North_This_Gridbox) .and. &
!              (Lon_Input(Ielem,Iline) >= Lon_West_This_Gridbox) .and. &
!              (Lon_Input(Ielem,Iline) < Lon_East_This_Gridbox)) then

          if (xdist < Max_Allowable_Distance) then  !check for pixels that fall within grid

          !----- if random, replace true xdist with a random number
          if (Random_Flag == sym%YES) then
              call random_number(xrand)
              xdist = xrand * Max_Allowable_Distance
          endif

          if (xdist < Min_Dist_Grid(Ilon,Ilat)) then  !check for pixels that are closer 

            if (Gap_Pixel_Mask_Input(ielem,iline) /= sym%YES) then
              Ielem_Output(Ilon,Ilat)  = Ielem
              Iline_Output(Ilon,Ilat)  = Iline
              Min_Dist_Grid(Ilon,Ilat) = xdist
            endif

            !update surrounding grid points that fall within the pixel's footprint but
            !have not already been filled in by a closest pixel
            if ((Ilon_min > 1) .and. (Ilat_min > 1)) then
              do Ilon_1 = Ilon_min, Ilon_max
                do Ilat_1 = Ilat_min, Ilat_max
                  !-- do not reset grid point closest to pixel
                  if ((Ilon_1 == Ilon) .and. (Ilat_1 == Ilat)) cycle
                  !-- fill in values for all grid points within pixel footprint
                  !-- ignore those already set previously
                  if (Min_Dist_Grid(Ilon_1,Ilat_1) == big_number .and.  &
                      Gap_Pixel_Mask_Input(Ielem,Iline) /= sym%YES) then
                    Ielem_Output(Ilon_1,Ilat_1)  = Ielem
                    Iline_Output(Ilon_1,Ilat_1)  = Iline
                  endif
                end do
              end do
            endif

          end if   !end check on min distance
          end if   !end check on max distance
        end if     ! valid regrid

      end if      !dateline check

    end if      !latitude check longitude check

  end do
end do

deallocate(Min_Dist_Grid)
!--- make grid lat and lons
  do Ilat = 1, Nlat
    do Ilon = 1, Nlon
       Lat_Output(Ilon,Ilat) = Lat_South + (Ilat-1)*dlat
       Lon_Output(Ilon,Ilat) = Lon_West + (Ilon-1)*dlon
    end do
  end do

  where(Lon_Output > 180.0)
     Lon_Output = Lon_Output - 360.0
  end where

  end subroutine REGRID

!====================================================================
! subroutine Name: indices_to_position
!
! Function: 
!       Converts grid indicies to lat/lon position
!
! Description:
!   This subroutine is a routine to convert the grid index to a lat/lon
!
!====================================================================
subroutine indices_to_position(Ilat,Ilon, lat, lon,  &
                               Lat_South,   &
                               dlat, Nlat, Lon_West, dlon, Nlon)
      integer(kind=int4), intent(in):: Ilat
      integer(kind=int4), intent(in):: Ilon
      integer(kind=int4), intent(in):: Lat_South
      integer(kind=int4), intent(in):: dlat
      integer(kind=int4), intent(in):: Nlat
      integer(kind=int4), intent(in):: Lon_West
      integer(kind=int4), intent(in):: dlon
      integer(kind=int4), intent(in):: Nlon
      real(kind=real4), intent(out):: lat
      real(kind=real4), intent(out):: lon
   
      !--- initialize
      lat = -999.0
      lon = -999.0

      if ((Ilat >= 1) .and. (Ilat <= Nlat) .and. (Ilon >= 1) .and. (Ilon <= Nlon)) then
       lat = Lat_South + (Ilat-1)*dlat
       lon = Lon_West + (Ilon-1)*dlon
      endif

      if (lon > 180.0) then 
           lon = lon - 360.0
      endif
end subroutine indices_to_position

!====================================================================
! subroutine Name: POSITION_TO_INDICES
!
! Function: 
!       Converts lat/lon position to grid indicies 
!
! Description:
!   This subroutine is a routine to convert lat/lon position to grid index
!
!====================================================================

subroutine POSITION_TO_INDICES(Ilat,Ilon, lat, lon,  &
                               Lat_South, Lat_North, &
                               dlat, Nlat, Lon_West, Lon_East, &
                               dlon, Nlon, distance)

       real(kind=real4), intent(in):: lat
       real(kind=real4), intent(in):: lon
       real(kind=real4), intent(in):: Lat_South
       real(kind=real4), intent(in):: Lat_North
       real(kind=real4), intent(in):: dlat
       integer(kind=int4), intent(in):: Nlat
       real(kind=real4), intent(in):: Lon_West
       real(kind=real4), intent(in):: Lon_East
       real(kind=real4), intent(in):: dlon
       integer(kind=int4), intent(in):: Nlon
       integer(kind=int4), intent(out):: Ilat
       integer(kind=int4), intent(out):: Ilon
       real(kind=real4), intent(out):: distance
       real(kind=real4):: lon_abs
       real(kind=real4):: Lon_East_abs
       real(kind=real4):: lat_abs

       Ilat = -999
       Ilon = -999

      Lon_East_abs = Lon_East
      if (Lon_East_abs < 0.0) Lon_East_abs = Lon_East_abs + 360.0
      lon_abs = lon
      if (lon_abs < 0.0) lon_abs = lon_abs + 360.0

      if ((lat >= Lat_South) .and. (lat <= Lat_North)) then
         Ilat = nint((lat - Lat_South) / dlat) + 1
         Ilat = min(Nlat,max(1,Ilat))
      endif

     !-- -no date line
     if (Lon_West  < Lon_East) then
      if ((lon >= Lon_West) .and. (lon <= Lon_East)) then
        Ilon = nint((lon - Lon_West) / dlon) + 1
        Ilon = min(Nlon,max(1,Ilon))
      end if
     else 
      if (lon < 0.0) lon_abs = lon + 360.0
      if ((lon_abs <= Lon_East_abs) .and. (lon_abs >= Lon_West)) then
        Ilon = nint((lon_abs - Lon_West) / dlon) + 1
        Ilon = min(Nlon,max(1,Ilon))
      endif
     end if

     distance = 999.9
     if ((Ilon >= 1) .and. (Ilat >= 1)) then
             lat_abs = Lat_South + (Ilat-1)*dlat
             Lon_East_abs = Lon_West + (Ilon-1)*dlon
             lon_abs = lon
             if (Lon_West > 0.0 .and. lon < 0.0) lon_abs = lon + 360.0
             distance = (lat-lat_abs)**2 + (lon_abs-Lon_East_abs)**2
     end if

end subroutine POSITION_TO_INDICES


!---------------------------------------------------------------
! subroutine Name: DEFINE_SDS
!
! Function:
!       Define the SDS's of the Level 2b file
!---------------------------------------------------------------
subroutine DEFINE_SDS_RANK1(Sd_Id,            &
                            Sds_Dims,         &
                            Sds_Chunk_Size,      &
                            Sds)

    integer, intent(in):: Sd_Id
    integer, dimension(1), intent(in):: Sds_Dims
    integer, dimension(1),intent(in):: Sds_Chunk_Size
    type(Sds_Struct), intent(inout):: Sds
    character(len=100):: Sds_Name_Temp

    integer:: sfcreate
    integer:: sfsnatt
    integer:: sfscatt
    integer:: sfschnk
    integer:: Istatus_Sum

   

    Istatus_Sum = 0

     Sds_Name_Temp = trim(Sds%Variable_Name)

     Sds%Id_Output = sfcreate(Sd_Id,trim(Sds_Name_Temp),Sds%Data_Type,Sds%rank,Sds_Dims)
     Istatus_Sum = sfsnatt(Sds%Id_Output, "SCALED", DFNT_INT8, 1, Sds%Scaling_Type) + Istatus_Sum
     if (Sds%Units /= "none") then
       Istatus_Sum = sfscatt(Sds%Id_Output, "units", DFNT_CHAR8, len_trim(Sds%Units), trim(Sds%Units)) + Istatus_Sum
     endif
     if (Sds%Standard_Name /= "not specified") then
       Istatus_Sum = sfscatt(Sds%Id_Output, "standard_name", DFNT_CHAR8, len_trim(Sds%Standard_Name),  &
                           trim(Sds%Standard_Name)) + Istatus_Sum
     endif
     Istatus_Sum = sfscatt(Sds%Id_Output, "long_name", DFNT_CHAR8, len_trim(Sds%Long_Name),  &
                           trim(Sds%Long_Name)) + Istatus_Sum

!    Dim_Id = sfdimid(Sds%Id_Output, 0)
!    Istatus_Sum = sfsdmname(Dim_Id,"longitude index") + Istatus_Sum

!    Dim_Id = sfdimid(Sds%Id_Output, 1)
!    Istatus_Sum = sfsdmname(Dim_Id,"latitude index") + Istatus_Sum

     Istatus_Sum = sfschnk(Sds%Id_Output,Sds_Chunk_Size,Comp_Type,Comp_Prm) + Istatus_Sum

     call WRITE_SCALING_ATTRIBUTES(Sds,Istatus_Sum)

     Sds%Status_Flag = Istatus_Sum

  end subroutine DEFINE_SDS_RANK1
!---------------------------------------------------------------
! subroutine Name: DEFINE_SDS_RANK2
!
! Function:
!       Define the SDS's of the Level 2b file
!---------------------------------------------------------------
subroutine DEFINE_SDS_RANK2(Sd_Id,            &
                            Sds_Dims,         &
                            Sds_Chunk_Size,   &
                            Sds)

    integer, intent(in):: Sd_Id
    type(Sds_Struct), intent(inout):: Sds
    integer, dimension(2), intent(in):: Sds_Dims
    integer, dimension(2),intent(in):: Sds_Chunk_Size
    character(len=100):: Sds_Name_Temp

    integer:: Dim_Id

    integer:: sfcreate
    integer:: sfsnatt
    integer:: sfscatt
    integer:: sfdimid
    integer:: sfsdmname
    integer:: sfschnk
    integer:: Istatus_Sum

   
    Istatus_Sum = 0

    Sds%Rank = 2 

    Sds_Name_Temp = trim(Sds%Variable_Name)

    Sds%Id_Output = sfcreate(Sd_Id,trim(Sds_Name_Temp),Sds%Data_Type,Sds%Rank,Sds_Dims)
    Istatus_Sum = sfsnatt(Sds%Id_Output, "SCALED", DFNT_INT8, 1, Sds%Scaling_Type) + Istatus_Sum
    if (Sds%Units /= "none") then
       Istatus_Sum = sfscatt(Sds%Id_Output, "units", DFNT_CHAR8,len_trim(Sds%Units), trim(Sds%Units)) + Istatus_Sum
    endif
    if (Sds%Standard_Name /= "not specified") then
      Istatus_Sum = sfscatt(Sds%Id_Output, "standard_name", DFNT_CHAR8, len_trim(Sds%Standard_Name),  &
                           trim(Sds%Standard_Name)) + Istatus_Sum
    endif
    Istatus_Sum = sfscatt(Sds%Id_Output, "long_name", DFNT_CHAR8, len_trim(Sds%Long_Name),  &
                           trim(Sds%Long_Name)) + Istatus_Sum
    Istatus_Sum = sfscatt(Sds%Id_Output, "coordinates", DFNT_CHAR8, len_trim(Coordinates_String),  &
                           trim(Coordinates_String)) + Istatus_Sum

    Dim_Id = sfdimid(Sds%Id_Output, 0)
    Istatus_Sum = sfsdmname(Dim_Id,"longitude") + Istatus_Sum

    Dim_Id = sfdimid(Sds%Id_Output, 1)
    Istatus_Sum = sfsdmname(Dim_Id,"latitude") + Istatus_Sum

    Istatus_Sum = sfschnk(Sds%Id_Output,Sds_Chunk_Size,Comp_Type,Comp_Prm) + Istatus_Sum

    !--- set attributes that communicate scaling
    call WRITE_SCALING_ATTRIBUTES(Sds,Istatus_Sum)

    !--- determine if flag attributes are used, and copy them over if so.
    call WRITE_FLAG_ATTRIBUTES(Sds%Id_Output,Sds_Name_Temp,Istatus_Sum)
  
    Sds%Status_Flag = Istatus_Sum

  end subroutine DEFINE_SDS_RANK2
!---------------------------------------------------------------
! subroutine Name: DEFINE_SDS_RANK3
!
! Function:
!       Define the SDS's of the Level 2b file
!---------------------------------------------------------------
subroutine DEFINE_SDS_RANK3(Sd_Id,            &
                            Sds_Dims,         &
                            Sds_Chunk_Size,   &
                            Sds)

    integer, intent(in):: Sd_Id
    type(Sds_Struct), intent(inout):: Sds
    integer, dimension(3), intent(in):: Sds_Dims
    integer, dimension(3),intent(in):: Sds_Chunk_Size
    character(len=100):: Sds_Name_Temp

    integer:: Dim_Id

    integer:: sfcreate
    integer:: sfsnatt
    integer:: sfscatt
    integer:: sfdimid
    integer:: sfsdmname
    integer:: sfschnk
    integer:: Istatus_Sum

    Istatus_Sum = 0

    Sds%Rank = 3  !This due to added time dimesion for NCDC

    Sds_Name_Temp = trim(Sds%Variable_Name)

    Sds%Id_Output = sfcreate(Sd_Id,trim(Sds_Name_Temp),Sds%Data_Type,Sds%Rank,Sds_Dims)
    Istatus_Sum = sfsnatt(Sds%Id_Output, "SCALED", DFNT_INT8, 1, Sds%Scaling_Type) + Istatus_Sum
    if (Sds%Units /= "none") then
       Istatus_Sum = sfscatt(Sds%Id_Output, "units", DFNT_CHAR8,len_trim(Sds%Units), trim(Sds%Units)) + Istatus_Sum
    endif
    if (Sds%Standard_Name /= "not specified") then
      Istatus_Sum = sfscatt(Sds%Id_Output, "standard_name", DFNT_CHAR8, len_trim(Sds%Standard_Name),  &
                           trim(Sds%Standard_Name)) + Istatus_Sum
    endif
    Istatus_Sum = sfscatt(Sds%Id_Output, "long_name", DFNT_CHAR8, len_trim(Sds%Long_Name),  &
                           trim(Sds%Long_Name)) + Istatus_Sum
    Istatus_Sum = sfscatt(Sds%Id_Output, "coordinates", DFNT_CHAR8, len_trim(Coordinates_String),  &
                           trim(Coordinates_String)) + Istatus_Sum

    Dim_Id = sfdimid(Sds%Id_Output, 0)
    Istatus_Sum = sfsdmname(Dim_Id,"time") + Istatus_Sum

    Dim_Id = sfdimid(Sds%Id_Output, 1)
    Istatus_Sum = sfsdmname(Dim_Id,"longitude") + Istatus_Sum

    Dim_Id = sfdimid(Sds%Id_Output, 2)
    Istatus_Sum = sfsdmname(Dim_Id,"latitude") + Istatus_Sum

    Istatus_Sum = sfschnk(Sds%Id_Output,Sds_Chunk_Size,Comp_Type,Comp_Prm) + Istatus_Sum

    !--- set attributes that communicate scaling
    call WRITE_SCALING_ATTRIBUTES(Sds,Istatus_Sum)

    !--- determine if flag attributes are used, and copy them over if so.
    call WRITE_FLAG_ATTRIBUTES(Sds%Id_Output,Sds_Name_Temp,Istatus_Sum)
  
    Sds%Status_Flag = Istatus_Sum

  end subroutine DEFINE_SDS_RANK3
!------------------------------------------------------------------------------------------------
! subroutine Name: READ_SDS_RANK1
!
! Function:
!     Reads 1dimensional SDSs from Level 2b files. 
!-----------------------------------------------------------------------------------------------
  subroutine READ_SDS_RANK1(Sd_Id,     &
                            Scaled_Sds_Data,   &
                            Sds)

    integer, intent(in):: Sd_Id
    type(Sds_Struct), intent(inout):: Sds
    real(kind=real4), dimension(:), intent(out):: Scaled_Sds_Data
    integer, dimension(1):: Sds_Dims
    integer, dimension(1):: Sds_Start
    integer, dimension(1):: Sds_Stride
    integer, dimension(1):: Sds_Edges
    integer(kind=int1), dimension(:), allocatable:: Temp_I1
    integer(kind=int2), dimension(:), allocatable:: Temp_I2
    integer(kind=int4), dimension(:), allocatable:: Temp_I4
    real(kind=real4), dimension(:), allocatable:: Temp_R4
    integer:: Num_Attrs
    integer:: sfn2index
    integer:: sfselect
    integer:: sfginfo
    integer:: sfrdata
    integer:: sfendacc
    integer:: Istatus
    character(100):: Sds_Name_Temp

    Istatus = 0

    !--- open Sds for reading
    Sds%Id_Input = sfselect(Sd_Id, sfn2index(Sd_Id,trim(Sds%Variable_Name)))
    if ( Sds%Id_Input < 0 ) then
       Sds%Data_Type = -999
       return
    endif
    Istatus = sfginfo(Sds%Id_Input, Sds_Name_Temp, Sds%rank, Sds_Dims, Sds%Data_Type, Num_Attrs) + Istatus
    Sds_Start = (/ 0 /)
    Sds_Stride = (/ 1 /)
    Sds_Edges = Sds_Dims

    !--- read sds scaling attributes
    call READ_SCALING_ATTRIBUTES(Sds,Istatus)

    !--- allocate arrays for holding data, read data and store in output array
    if (Sds%Data_Type == DFNT_INT8) then
       allocate(Temp_I1(Sds_Dims(1)))
       Istatus = sfrdata(Sds%Id_Input, Sds_Start, Sds_Stride, Sds_Edges, Temp_I1) + Istatus
       Scaled_Sds_Data = real(Temp_I1)
    elseif (Sds%Data_Type == DFNT_INT16) then
       allocate(Temp_I2(Sds_Dims(1)))
       Istatus = sfrdata(Sds%Id_Input, Sds_Start, Sds_Stride, Sds_Edges, Temp_I2) + Istatus
       Scaled_Sds_Data = real(Temp_I2)
    elseif (Sds%Data_Type == DFNT_INT32) then
       allocate(Temp_I4(Sds_Dims(1)))
       Istatus = sfrdata(Sds%Id_Input, Sds_Start, Sds_Stride, Sds_Edges, Temp_I4) + Istatus
       Scaled_Sds_Data = real(Temp_I4)
    elseif (Sds%Data_Type == DFNT_FLOAT32) then
       allocate(Temp_R4(Sds_Dims(1)))
       Istatus = sfrdata(Sds%Id_Input, Sds_Start, Sds_Stride, Sds_Edges, Temp_R4) + Istatus
       Scaled_Sds_Data = Temp_R4
    else
       print *, "Possibly fatal error at location 1:"
       print *, "data type was ", Sds%Data_Type
       print *, MOD_PROMPT, "attempt to read unsupported data type, stopping "
       Scaled_Sds_Data = Missing_Value_Real4
       Sds%Data_Type = -999
       return
    endif

    !---deallocate temp arrays
    if (allocated(Temp_I1)) deallocate(Temp_I1)
    if (allocated(Temp_I2)) deallocate(Temp_I2)
    if (allocated(Temp_I4)) deallocate(Temp_I4)
    if (allocated(Temp_R4)) deallocate(Temp_R4)

    !--- close Sds
    Istatus = sfendacc(Sds%Id_Input) + Istatus

   end subroutine READ_SDS_RANK1


!------------------------------------------------------------------------------------------------
! subroutine Name: READ_SDS_RANK2
!
! Function:
!     Reads 2-dimensional SDSs from Level 2b files. 
!-----------------------------------------------------------------------------------------------
  subroutine READ_SDS_RANK2(Sd_Id,     &
                            Scaled_Sds_Data,   &
                            Sds)

    integer, intent(in):: Sd_Id
    type(Sds_Struct), intent(inout):: Sds
    real(kind=real4), dimension(:,:), intent(out):: Scaled_Sds_Data
    integer, dimension(2):: Sds_Dims
    integer, dimension(2):: Sds_Start
    integer, dimension(2):: Sds_Stride
    integer, dimension(2):: Sds_Edges
    integer(kind=int1), dimension(:,:), allocatable:: Temp_I1
    integer(kind=int2), dimension(:,:), allocatable:: Temp_I2
    integer(kind=int4), dimension(:,:), allocatable:: Temp_I4
    real(kind=real4), dimension(:,:), allocatable:: Temp_R4
    integer:: Num_Attrs
    integer:: sfn2index
    integer:: sfselect
    
    
   
    integer:: sfginfo
    integer:: sfrdata
    integer:: sfendacc
    integer:: Istatus
    character(100):: Sds_Name_Temp

    Istatus = 0

    !--- open Sds for reading
    Sds%Id_Input = sfselect(Sd_Id, sfn2index(Sd_Id,trim(Sds%Variable_Name)))
    if ( Sds%Id_Input < 0 ) then
       Sds%Data_Type = -999
       return
    endif
    Istatus = sfginfo(Sds%Id_Input, Sds_Name_Temp, Sds%rank, Sds_Dims, Sds%Data_Type, Num_Attrs) + Istatus
    Sds_Start = (/ 0, 0 /)
    Sds_Stride = (/ 1, 1 /)
    Sds_Edges = Sds_Dims

    !--- read sds scaling attributes
    call READ_SCALING_ATTRIBUTES(Sds,Istatus)

    !--- allocate arrays for holding data, read data and store in output array
    if (Sds%Data_Type == DFNT_INT8) then
       allocate(Temp_I1(Sds_Dims(1),Sds_Dims(2)))
       Istatus = sfrdata(Sds%Id_Input, Sds_Start, Sds_Stride, Sds_Edges, Temp_I1) + Istatus
       Scaled_Sds_Data = real(Temp_I1)
    elseif (Sds%Data_Type == DFNT_INT16) then
       allocate(Temp_I2(Sds_Dims(1), Sds_Dims(2)))
       Istatus = sfrdata(Sds%Id_Input, Sds_Start, Sds_Stride, Sds_Edges, Temp_I2) + Istatus
       Scaled_Sds_Data = real(Temp_I2)
    elseif (Sds%Data_Type == DFNT_INT32) then
       allocate(Temp_I4(Sds_Dims(1),Sds_Dims(2)))
       Istatus = sfrdata(Sds%Id_Input, Sds_Start, Sds_Stride, Sds_Edges, Temp_I4) + Istatus
       Scaled_Sds_Data = real(Temp_I4)
    elseif (Sds%Data_Type == DFNT_FLOAT32) then
       allocate(Temp_R4(Sds_Dims(1),Sds_Dims(2)))
       Istatus = sfrdata(Sds%Id_Input, Sds_Start, Sds_Stride, Sds_Edges, Temp_R4) + Istatus
       Scaled_Sds_Data = Temp_R4
    else
       print *, "Possibly fatal error at location 2:"
       print *, "data type was ", Sds%Data_Type
       print *, "sds complete value is "
       print *, Sds
       print *, MOD_PROMPT, "attempt to read unsupported data type, stopping"
       Scaled_Sds_Data = Missing_Value_Real4 
       Sds%Data_Type = -999
       return
    endif

    !---deallocate temp arrays
    if (allocated(Temp_I1)) deallocate(Temp_I1)
    if (allocated(Temp_I2)) deallocate(Temp_I2)
    if (allocated(Temp_I4)) deallocate(Temp_I4)
    if (allocated(Temp_R4)) deallocate(Temp_R4)

    !--- close Sds
    Istatus = sfendacc(Sds%Id_Input) + Istatus

   end subroutine READ_SDS_RANK2

!------------------------------------------------------------------------------------------------
! subroutine Name: READ_SDS_RANK3
!
! Function:
!     Reads 3-dimensional SDSs from Level 2b files and reshape to 2 dimensions
!     (remove time coordinate). 
!-----------------------------------------------------------------------------------------------
  subroutine READ_SDS_RANK3(Sd_Id,     &
                            Scaled_Sds_Data,   &
                            Sds)

    integer, intent(in):: Sd_Id
    type(Sds_Struct), intent(inout):: Sds
    real(kind=real4), dimension(:,:,:), intent(out):: Scaled_Sds_Data
    integer, dimension(3):: Sds_Dims
    integer, dimension(3):: Sds_Start
    integer, dimension(3):: Sds_Stride
    integer, dimension(3):: Sds_Edges
    integer(kind=int1), dimension(:,:,:), allocatable:: Temp_I1
    integer(kind=int2), dimension(:,:,:), allocatable:: Temp_I2
    integer(kind=int4), dimension(:,:,:), allocatable:: Temp_I4
    real(kind=real4), dimension(:,:,:), allocatable:: Temp_R4
    integer:: Num_Attrs
    integer:: sfn2index
    integer:: sfselect
    
    integer:: sfginfo
    integer:: sfrdata
    integer:: sfendacc
    integer:: Istatus
    character(100):: Sds_Name_Temp

    Istatus = 0

    !--- open Sds for reading
    Sds%Id_Input = sfselect(Sd_Id, sfn2index(Sd_Id,trim(Sds%Variable_Name)))
    if ( Sds%Id_Input < 0 ) then
       Sds%Data_Type = -999
       return
    endif
    Istatus = sfginfo(Sds%Id_Input, Sds_Name_Temp, Sds%rank, Sds_Dims, Sds%Data_Type, Num_Attrs) + Istatus
    Sds_Start = (/ 0, 0, 0 /)
    Sds_Stride = (/ 1, 1, 1 /)
    Sds_Edges = Sds_Dims

    !--- read sds scaling attributes
    call READ_SCALING_ATTRIBUTES(Sds,Istatus)

    !--- allocate arrays for holding data, read data and store in output array
    if (Sds%Data_Type == DFNT_INT8) then
       allocate(Temp_I1(Sds_Dims(1),Sds_Dims(2),Sds_Dims(3)))
       Istatus = sfrdata(Sds%Id_Input, Sds_Start, Sds_Stride, Sds_Edges, Temp_I1) + Istatus
       Scaled_Sds_Data = real(Temp_I1)
    elseif (Sds%Data_Type == DFNT_INT16) then
       allocate(Temp_I2(Sds_Dims(1), Sds_Dims(2),Sds_Dims(3)))
       Istatus = sfrdata(Sds%Id_Input, Sds_Start, Sds_Stride, Sds_Edges, Temp_I2) + Istatus
       Scaled_Sds_Data = real(Temp_I2)
    elseif (Sds%Data_Type == DFNT_INT32) then
       allocate(Temp_I4(Sds_Dims(1),Sds_Dims(2),Sds_Dims(3)))
       Istatus = sfrdata(Sds%Id_Input, Sds_Start, Sds_Stride, Sds_Edges, Temp_I4) + Istatus
       Scaled_Sds_Data = real(Temp_I4)
    elseif (Sds%Data_Type == DFNT_FLOAT32) then
       allocate(Temp_R4(Sds_Dims(1),Sds_Dims(2),Sds_Dims(3)))
       Istatus = sfrdata(Sds%Id_Input, Sds_Start, Sds_Stride, Sds_Edges, Temp_R4) + Istatus
       Scaled_Sds_Data = real(Temp_R4)
    else
       print *, "Possibly fatal error at location 2:"
       print *, "data type was ", Sds%Data_Type
       print *, "sds complete value is "
       print *, Sds
       print *, MOD_PROMPT, "attempt to read unsupported data type, stopping"
       Scaled_Sds_Data = Missing_Value_Real4
       Sds%Data_Type = -999
       return
    endif

    !---deallocate temp arrays
    if (allocated(Temp_I1)) deallocate(Temp_I1)
    if (allocated(Temp_I2)) deallocate(Temp_I2)
    if (allocated(Temp_I4)) deallocate(Temp_I4)
    if (allocated(Temp_R4)) deallocate(Temp_R4)

    !--- close Sds
    Istatus = sfendacc(Sds%Id_Input) + Istatus

   end subroutine READ_SDS_RANK3
!------------------------------------------------------------------------------------------------
! subroutine Name: UNSCALE_SDS_RANK1
!
! Function:
!     Unscales scaled SDS data from level2b files
!-----------------------------------------------------------------------------------------------
   subroutine UNSCALE_SDS_RANK1(Sds, Scaled_Data, Unscaled_Data)
      type(Sds_Struct), intent(in):: Sds
      real(kind=real4), dimension(:), intent(in):: Scaled_Data
      real(kind=real4), dimension(:), intent(out):: Unscaled_Data

      !--- unscale Sds

      if (Sds%Scaling_Type /= sym%NO_SCALING) then

        !---- linear
        if (Sds%Scaling_Type == sym%LINEAR_SCALING) then
           Unscaled_Data = Sds%Scale_Factor * real(Scaled_Data) + Sds%Add_Offset
        endif

        where (Scaled_Data == Sds%Fill_Value)
         Unscaled_Data = Missing_Value_Real4
        endwhere

      else

        Unscaled_Data = Scaled_Data

      endif

   end subroutine UNSCALE_SDS_RANK1
!------------------------------------------------------------------------------------------------
! subroutine Name: UNSCALE_SDS_RANK2
!
! Function:
!     Unscales scaled SDS data from level2b files
!-----------------------------------------------------------------------------------------------
   subroutine UNSCALE_SDS_RANK2(Sds, Scaled_Data, Unscaled_Data)
      type(Sds_Struct), intent(in):: Sds
      real(kind=real4), dimension(:,:), intent(in):: Scaled_Data
      real(kind=real4), dimension(:,:), intent(out):: Unscaled_Data

      !--- unscale Sds

      if (Sds%Scaling_Type /= sym%NO_SCALING) then

        !---- linear
        if (Sds%Scaling_Type == sym%LINEAR_SCALING) then
           Unscaled_Data = Sds%Scale_Factor * real(Scaled_Data) + Sds%Add_Offset
        endif

        where (Scaled_Data == Sds%Fill_Value)
         Unscaled_Data = Missing_Value_Real4
        endwhere

      else

        Unscaled_Data = Scaled_Data

      endif

   end subroutine UNSCALE_SDS_RANK2

!------------------------------------------------------------------------------------------------
! subroutine Name: UNSCALE_SDS_RANK3
!
! Function:
!     Unscales scaled SDS data from level2b files
!-----------------------------------------------------------------------------------------------
   subroutine UNSCALE_SDS_RANK3(Sds, Scaled_Data, Unscaled_Data)
      type(Sds_Struct), intent(in):: Sds
      real(kind=real4), dimension(:,:,:), intent(in):: Scaled_Data
      real(kind=real4), dimension(:,:,:), intent(out):: Unscaled_Data

      !--- unscale Sds
      if (Sds%Scaling_Type /= sym%NO_SCALING) then

        !---- linear
        if (Sds%Scaling_Type == sym%LINEAR_SCALING) then
           Unscaled_Data = Sds%Scale_Factor * real(Scaled_Data) + Sds%Add_Offset
        endif

        !--- handle missing
        where (Scaled_Data == Sds%Fill_Value)
         Unscaled_Data = Sds%Unscaled_Missing
        endwhere
        
      else
        Unscaled_Data = Scaled_Data
      endif

   end subroutine UNSCALE_SDS_RANK3
!------------------------------------------------------------------------------------------------
! subroutine Name: SCALE_SDS_RANK1
!
! Function:
!     Scales unscaled SDS data from level2b files
!-----------------------------------------------------------------------------------------------
   subroutine SCALE_SDS_RANK1(Sds, Unscaled_Data, Scaled_Data)
      type(Sds_Struct), intent(in):: Sds
      real(kind=real4), dimension(:), intent(in):: Unscaled_Data
      real(kind=real4), dimension(:), intent(out):: Scaled_Data

      !--- scale Sds
      if (Sds%Scaling_Type /= sym%NO_SCALING) then

        !---- linear
        if (Sds%Scaling_Type == sym%LINEAR_SCALING) then
           Scaled_Data = (Unscaled_Data - Sds%Add_Offset)/(Sds%Scale_Factor)
           where (Scaled_Data < Sds%Valid_Range(1))
             Scaled_Data = Sds%Valid_Range(1)
           endwhere
           where (Scaled_Data > Sds%Valid_Range(2))
             Scaled_Data = Sds%Valid_Range(2)
           endwhere
        endif

        !--- set scaled missing values
        where (Unscaled_Data == Sds%Unscaled_Missing)
         Scaled_Data = Sds%Fill_Value
        endwhere

     else
       Scaled_Data = Unscaled_Data
     endif

   end subroutine SCALE_SDS_RANK1
!------------------------------------------------------------------------------------------------
! subroutine Name: SCALE_SDS_RANK2
!
! Function:
!     Scales unscaled SDS data from level2b files
!-----------------------------------------------------------------------------------------------
   subroutine SCALE_SDS_RANK2(Sds, Unscaled_Data, Scaled_Data)
      type(Sds_Struct), intent(in):: Sds
      real(kind=real4), dimension(:,:), intent(in):: Unscaled_Data
      real(kind=real4), dimension(:,:), intent(out):: Scaled_Data

      !--- scale Sds
      if (Sds%Scaling_Type /= sym%NO_SCALING) then

        !---- linear
        if (Sds%Scaling_Type == sym%LINEAR_SCALING) then
           Scaled_Data = (Unscaled_Data - Sds%Add_Offset)/(Sds%Scale_Factor)
           where (Scaled_Data < Sds%Valid_Range(1))
             Scaled_Data = Sds%Valid_Range(1)
           endwhere
           where (Scaled_Data > Sds%Valid_Range(2))
             Scaled_Data = Sds%Valid_Range(2)
           endwhere
        endif

        !--- set scaled missing values
        where (Unscaled_Data == Sds%Unscaled_Missing)
         Scaled_Data = Sds%Fill_Value
        endwhere

     else
       Scaled_Data = Unscaled_Data
     endif

   end subroutine SCALE_SDS_RANK2

!------------------------------------------------------------------------------------------------
! subroutine Name: SCALE_SDS_RANK3
!
! Function:
!     Scales unscaled SDS data from level2b files
!-----------------------------------------------------------------------------------------------
   subroutine SCALE_SDS_RANK3(Sds, Unscaled_Data, Scaled_Data)
      type(Sds_Struct), intent(in):: Sds
      real(kind=real4), dimension(:,:,:), intent(in):: Unscaled_Data
      real(kind=real4), dimension(:,:,:), intent(out):: Scaled_Data

      !--- scale Sds
      if (Sds%Scaling_Type /= sym%NO_SCALING) then

        !---- linear
        if (Sds%Scaling_Type == sym%LINEAR_SCALING) then
           Scaled_Data = (Unscaled_Data - Sds%Add_Offset)/(Sds%Scale_Factor)
           where (Scaled_Data < Sds%Valid_Range(1))
             Scaled_Data = Sds%Valid_Range(1)
           endwhere
           where (Scaled_Data > Sds%Valid_Range(2))
             Scaled_Data = Sds%Valid_Range(2)
           endwhere
        endif

        !--- set scaled missing values
        where (Unscaled_Data == Sds%Unscaled_Missing)
         Scaled_Data = Sds%Fill_Value
        endwhere

     else
       Scaled_Data = Unscaled_Data
     endif

   end subroutine SCALE_SDS_RANK3


!-----------------------------------------------------------------------------
! subroutine Name: WRITE_SDS_RANK1
!
! Function:
!     Write 1 dimensional data to a HDF file
!-----------------------------------------------------------------------------
  subroutine WRITE_SDS_RANK1(Scaled_Sds_Data,   &
                             Sds)

    type(Sds_Struct), intent(in):: Sds
    real(kind=real4), dimension(:), intent(in):: Scaled_Sds_Data
    integer, dimension(1):: Sds_Dims
    integer, dimension(1):: Sds_Start
    integer, dimension(1):: Sds_Stride
    integer, dimension(1):: Sds_Edges
    integer(kind=int1), dimension(:), allocatable:: Temp_I1
    integer(kind=int2), dimension(:), allocatable:: Temp_I2
    integer(kind=int4), dimension(:), allocatable:: Temp_I4
    real(kind=real4), dimension(:), allocatable:: Temp_R4
    integer(kind=int4):: Istatus
    integer:: sfwdata

    Sds_Dims(1) = size(Scaled_Sds_Data,1)
    Sds_Start = (/ 0 /)
    Sds_Stride = (/ 1 /)
    Sds_Edges = Sds_Dims


    Istatus = 0

    Sds_Dims(1) = size(Scaled_Sds_Data,1)

    if (Sds%Data_Type == DFNT_INT8) then
       allocate(Temp_I1(Sds_Dims(1)))
       Temp_I1 = int(Scaled_Sds_Data,kind=int1)
       Istatus = sfwdata(Sds%Id_Output, Sds_Start, Sds_Stride, Sds_Edges, Temp_I1) + Istatus
    elseif (Sds%Data_Type == DFNT_INT16) then
       allocate(Temp_I2(Sds_Dims(1)))
       Temp_I2 = int(Scaled_Sds_Data,kind=int2)
       Istatus = sfwdata(Sds%Id_Output, Sds_Start, Sds_Stride, Sds_Edges, Temp_I2) + Istatus
    elseif (Sds%Data_Type == DFNT_INT32) then
       allocate(Temp_I4(Sds_Dims(1)))
       Temp_I4 = int(Scaled_Sds_Data,kind=int4)
       Istatus = sfwdata(Sds%Id_Output, Sds_Start, Sds_Stride, Sds_Edges, Temp_I4) + Istatus
    elseif (Sds%Data_Type == DFNT_FLOAT32) then
       allocate(Temp_R4(Sds_Dims(1)))
       Temp_R4 = real(Scaled_Sds_Data,kind=real4)
       Istatus = sfwdata(Sds%Id_Output, Sds_Start, Sds_Stride, Sds_Edges, Temp_R4) + Istatus
    endif


   !---deallocate temp arrays
   if (allocated(Temp_I1)) deallocate(Temp_I1)
   if (allocated(Temp_I2)) deallocate(Temp_I2)
   if (allocated(Temp_I4)) deallocate(Temp_I4)
   if (allocated(Temp_R4)) deallocate(Temp_R4)

   end subroutine WRITE_SDS_RANK1

!------------------------------------------------------------------------------------------------
! subroutine Name: WRITE_SDS_RANK2
!
! Function:
!     Write 2 dimensional data to a HDF file
!-----------------------------------------------------------------------------------------------
    subroutine WRITE_SDS_RANK2(Scaled_Sds_Data,Sds)

     type(Sds_Struct), intent(in):: Sds
     real(kind=real4), dimension(:,:), intent(in):: Scaled_Sds_Data
     integer, dimension(2):: Sds_Dims
     integer, dimension(2):: Sds_Start
     integer, dimension(2):: Sds_Stride
     integer, dimension(2):: Sds_Edges
     integer(kind=int1), dimension(:,:), allocatable:: Temp_I1
     integer(kind=int2), dimension(:,:), allocatable:: Temp_I2
     integer(kind=int4), dimension(:,:), allocatable:: Temp_I4
     real(kind=real4), dimension(:,:), allocatable:: Temp_R4
     integer(kind=int4):: Istatus
     integer:: sfwdata

     Sds_Dims(1) = size(Scaled_Sds_Data,1)
     Sds_Dims(2) = size(Scaled_Sds_Data,2)
     Sds_Start = (/ 0, 0 /)
     Sds_Stride = (/ 1, 1 /)
     Sds_Edges = Sds_Dims

     Istatus = 0

     if (Sds%Data_Type == DFNT_INT8) then
        allocate(Temp_I1(Sds_Dims(1),Sds_Dims(2)))
        Temp_I1 = int(Scaled_Sds_Data,kind=int1)
        Istatus = sfwdata(Sds%Id_Output, Sds_Start, Sds_Stride, Sds_Edges, Temp_I1) + Istatus
     elseif (Sds%Data_Type == DFNT_INT16) then
        allocate(Temp_I2(Sds_Dims(1),Sds_Dims(2)))
        Temp_I2 = int(Scaled_Sds_Data,kind=int2)
        Istatus = sfwdata(Sds%Id_Output, Sds_Start, Sds_Stride, Sds_Edges, Temp_I2) + Istatus
     elseif (Sds%Data_Type == DFNT_INT32) then
        allocate(Temp_I4(Sds_Dims(1),Sds_Dims(2)))
        Temp_I4 = int(Scaled_Sds_Data,kind=int4)
        Istatus = sfwdata(Sds%Id_Output, Sds_Start, Sds_Stride, Sds_Edges, Temp_I4) + Istatus
     elseif (Sds%Data_Type == DFNT_FLOAT32) then
        allocate(Temp_R4(Sds_Dims(1),Sds_Dims(2)))
        Temp_R4 = real(Scaled_Sds_Data,kind=real4)
        Istatus = sfwdata(Sds%Id_Output, Sds_Start, Sds_Stride, Sds_Edges, Temp_R4) + Istatus
     endif

   !---deallocate temp arrays
   if (allocated(Temp_I1)) deallocate(Temp_I1)
   if (allocated(Temp_I2)) deallocate(Temp_I2)
   if (allocated(Temp_I4)) deallocate(Temp_I4)
   if (allocated(Temp_R4)) deallocate(Temp_R4)

end subroutine WRITE_SDS_RANK2

!------------------------------------------------------------------------------------------------
! subroutine Name: WRITE_SDS_RANK3
!
! Function:
!     Write 3 dimensional data to a HDF file
!-----------------------------------------------------------------------------------------------
    subroutine WRITE_SDS_RANK3(Scaled_Sds_Data,Sds)

     type(Sds_Struct), intent(in):: Sds
     real(kind=real4), dimension(:,:,:), intent(in):: Scaled_Sds_Data
     integer, dimension(3):: Sds_Dims
     integer, dimension(3):: Sds_Start
     integer, dimension(3):: Sds_Stride
     integer, dimension(3):: Sds_Edges
     integer(kind=int1), dimension(:,:,:), allocatable:: Temp_I1
     integer(kind=int2), dimension(:,:,:), allocatable:: Temp_I2
     integer(kind=int4), dimension(:,:,:), allocatable:: Temp_I4
     real(kind=real4), dimension(:,:,:), allocatable:: Temp_R4
     integer(kind=int4):: Istatus
     integer:: sfwdata

     Sds_Dims(1) = size(Scaled_Sds_Data,1)
     Sds_Dims(2) = size(Scaled_Sds_Data,2)
     Sds_Dims(3) = size(Scaled_Sds_Data,3)
     Sds_Start = (/ 0, 0, 0 /)
     Sds_Stride = (/ 1, 1, 1 /)
     Sds_Edges = Sds_Dims

     Istatus = 0

     if (Sds%Data_Type == DFNT_INT8) then
        allocate(Temp_I1(Sds_Dims(1),Sds_Dims(2),Sds_Dims(3)))
        Temp_I1 = int(Scaled_Sds_Data,kind=int1)
        Istatus = sfwdata(Sds%Id_Output, Sds_Start, Sds_Stride, Sds_Edges, Temp_I1) + Istatus
     elseif (Sds%Data_Type == DFNT_INT16) then
        allocate(Temp_I2(Sds_Dims(1),Sds_Dims(2),Sds_Dims(3)))
        Temp_I2 = int(Scaled_Sds_Data,kind=int2)
        Istatus = sfwdata(Sds%Id_Output, Sds_Start, Sds_Stride, Sds_Edges, Temp_I2) + Istatus
     elseif (Sds%Data_Type == DFNT_INT32) then
        allocate(Temp_I4(Sds_Dims(1),Sds_Dims(2),Sds_Dims(3)))
        Temp_I4 = int(Scaled_Sds_Data,kind=int4)
        Istatus = sfwdata(Sds%Id_Output, Sds_Start, Sds_Stride, Sds_Edges, Temp_I4) + Istatus
     elseif (Sds%Data_Type == DFNT_FLOAT32) then
        allocate(Temp_R4(Sds_Dims(1),Sds_Dims(2),Sds_Dims(3)))
        Temp_R4 = real(Scaled_Sds_Data,kind=real4)
        Istatus = sfwdata(Sds%Id_Output, Sds_Start, Sds_Stride, Sds_Edges, Temp_R4) + Istatus
     endif

   !---deallocate temp arrays
   if (allocated(Temp_I1)) deallocate(Temp_I1)
   if (allocated(Temp_I2)) deallocate(Temp_I2)
   if (allocated(Temp_I4)) deallocate(Temp_I4)
   if (allocated(Temp_R4)) deallocate(Temp_R4)

end subroutine WRITE_SDS_RANK3

!------------------------------------------------------------------------------------------------
! subroutine Name: COPY_GLOBAL_ATTRIBUTES
!
! Function:
!    Copies the global attributes from one file and puts them in a different file
!-----------------------------------------------------------------------------------------------

subroutine COPY_GLOBAL_ATTRIBUTES(Sd_Id_Input,Sd_Id_Output)
   
   integer(kind=int4):: Sd_Id_Input
   integer(kind=int4):: Sd_Id_Output
   integer:: Istatus
   integer(kind=int4):: num_global_attrs
   integer(kind=int4):: Num_Sds_Input
   integer(kind=int4):: iattr
   character(len=100):: Attr_Name
   integer(kind=int4):: Data_Type
   integer(kind=int4):: count
   integer(kind=int4):: Attr_Buffer_i4
   integer(kind=int2):: Attr_Buffer_i2
   integer(kind=int1):: Attr_Buffer_i1
   real(kind=real4):: Attr_Buffer_r4
   character(len=500):: Attr_Buffer_char

   integer:: sffinfo
   integer:: sfgainfo
   integer:: sfrnatt
   integer:: sfsnatt
   
   
   Istatus = sffinfo(Sd_Id_Input,Num_Sds_Input,num_global_attrs)

   do iattr = 0, num_global_attrs-1

       Istatus = sfgainfo(Sd_Id_Input,iattr,Attr_Name,Data_Type,count)

       !--- skip certain level2 attributes in level2b
       if (Attr_Name == "FILENAME" .or. &
           Attr_Name == "L1B" .or. &
           Attr_Name == "NUMBER_OF_SCANS_LEVEL1B" .or. &
           Attr_Name == "NUMBER_OF_SCANS_LEVEL2" .or. &
           Attr_Name == "START_YEAR" .or.  &
           Attr_Name == "END_YEAR" .or.  &
           Attr_Name == "START_DAY" .or.  &
           Attr_Name == "END_DAY" .or.  &
           Attr_Name == "START_TIME" .or.  &
           Attr_Name == "END_TIME") then
           cycle
        endif

       if (Data_Type == DFNT_INT8) then
           Istatus = sfrnatt(Sd_Id_Input,iattr,Attr_Buffer_i1)
           Istatus = sfsnatt(Sd_Id_Output,trim(Attr_Name),Data_Type,count,Attr_Buffer_i1)
       endif

       if (Data_Type == DFNT_INT16) then
           Istatus = sfrnatt(Sd_Id_Input,iattr,Attr_Buffer_i2)
           Istatus = sfsnatt(Sd_Id_Output,trim(Attr_Name),Data_Type,count,Attr_Buffer_i2)
       endif

       if (Data_Type == DFNT_INT32) then
           Istatus = sfrnatt(Sd_Id_Input,iattr,Attr_Buffer_i4)
           Istatus = sfsnatt(Sd_Id_Output,trim(Attr_Name),Data_Type,count,Attr_Buffer_i4)
       endif

       if (Data_Type == DFNT_FLOAT32) then
           Istatus = sfrnatt(Sd_Id_Input,iattr,Attr_Buffer_r4)
           Istatus = sfsnatt(Sd_Id_Output,trim(Attr_Name),Data_Type,count,Attr_Buffer_r4)
       endif

       if (Data_Type == DFNT_CHAR) then
           Istatus = sfrnatt(Sd_Id_Input,iattr,Attr_Buffer_char)
           Istatus = sfsnatt(Sd_Id_Output,trim(Attr_Name),Data_Type,count,trim(Attr_Buffer_char))
       endif
   enddo

   end subroutine 

!------------------------------------------------------------------------------------------------
! subroutine Name: SUBSET_LEVEL2b
!
! Function:
!    Reads in a level-2b file and then write out a subset of it
!
! Notes:
!    This routine assumes there are only 1d and 2d arrays
!    This routine assumes the only 1d arrays are latitude and longitude
!   
!-----------------------------------------------------------------------------------------------
  subroutine SUBSET_LEVEL2b(File_Input, &
                            Path_Input, & 
                            File_Output, &
                            Path_Output, & 
                            Lon_West_Out_Temp, &
                            Lon_East_Out_Temp, &
                            Dlon_Out, &
                            Lat_South_Out_Temp, &
                            Lat_North_Out_Temp, &
                            DLat_Out, &
                            Data_Description_String, &
                            Sds_Output_Names)

  character(len=*), intent(in):: File_Input
  character(len=*), intent(in):: Path_Input
  character(len=*), intent(in):: File_Output
  character(len=*), intent(in):: Path_Output
  character(len=*), intent(in):: Data_Description_String
  character(len=*), dimension(:), intent(in):: Sds_Output_Names
  real(kind=real4), intent(in):: Lon_West_Out_Temp
  real(kind=real4), intent(in):: Lon_East_Out_Temp
  real(kind=real4), intent(in):: Dlon_Out
  real(kind=real4), intent(in):: Lat_South_Out_Temp
  real(kind=real4), intent(in):: Lat_North_Out_Temp
  real(kind=real4), intent(in):: DLat_Out

  real(kind=real4) :: Lon_West_Out
  real(kind=real4) :: Lon_East_Out
  real(kind=real4) :: Lat_South_Out
  real(kind=real4) :: Lat_North_Out

  integer(kind=int4):: Sd_Id_Input
  integer(kind=int4):: Sd_Id_Output

  type(Sds_Struct) :: Sds_Input
  type(Sds_Struct) :: Sds_Output

  integer(kind=int4) :: Nlon_Out
  integer(kind=int4) :: NLat_Out
  real(kind=real4):: Lon_West_In
  real(kind=real4):: Lon_East_In
  real(kind=real4):: Dlon_In
  integer(kind=int4):: Nlon_In
  real(kind=real4):: Lat_South_In
  real(kind=real4):: Lat_North_In
  real(kind=real4):: Dlat_In
  integer(kind=int4):: Nlat_In

  integer(kind=int4):: Lon_Stride
  integer(kind=int4):: Lat_Stride
  integer(kind=int4):: Lon_Start
  integer(kind=int4):: Lat_Start

  real(kind=real4), dimension(:), allocatable:: Scaled_Sds_Longitude_Input
  real(kind=real4), dimension(:), allocatable:: Scaled_Sds_Latitude_Input
  real(kind=real4), dimension(:,:), allocatable:: Scaled_Sds_Data_Input
  real(kind=real4), dimension(:), allocatable:: Scaled_Sds_Longitude_Output
  real(kind=real4), dimension(:), allocatable:: Scaled_Sds_Latitude_Output
  real(kind=real4), dimension(:,:), allocatable:: Scaled_Sds_Data_Output

  integer:: sfstart, sfend, sfscatt, sfsnatt, sfrnatt,  &
            sfginfo, sfselect, sffattr, sffinfo
  integer:: Isds, Isds_Inner, Num_Sds_Input, Num_Global_Attrs, Num_Sds_Output
  integer:: Istatus
  integer:: Ifound
  integer(kind=int4), dimension(2):: Sds_Dims

  Istatus = 0
  Num_Sds_Output = size(Sds_Output_Names)

  !---- prevent round-off
  Lon_West_Out = nint(Lon_West_Out_Temp*100)/100.0
  Lon_East_Out = nint(Lon_East_Out_Temp*100)/100.0
  Lat_North_Out = nint(Lat_North_Out_Temp*100)/100.0
  Lat_South_Out = nint(Lat_South_Out_Temp*100)/100.0

  !-----------------------------------------------------------------
  ! determine size of output grid and shift if necessary to align
  ! with the input grid
  !-----------------------------------------------------------------

  !------------------------------------------------------------------
  ! open input file for reading
  !------------------------------------------------------------------
  Sd_Id_Input = sfstart(trim(Path_Input)//trim(File_Input), DFACC_READ)


  !------------------------------------------------------------------
  ! open output file for writing 
  !------------------------------------------------------------------
  Sd_Id_Output = sfstart(trim(Path_Output)//trim(File_Output), DFACC_CREATE)


  !--- determine number of sds's in these files (assume the same)
  Istatus = sffinfo(Sd_Id_Input,Num_Sds_Input,Num_Global_Attrs) + Istatus

  !------------------------------------------------------------------
  ! copy attributes from input to output file 
  !------------------------------------------------------------------
  call COPY_GLOBAL_ATTRIBUTES(Sd_Id_Input,Sd_Id_Output)

  !------------------------------------------------------------------
  !  write data description attribute
  !------------------------------------------------------------------
  Istatus = sfscatt(Sd_Id_Output,"DATA_DESCRIPTION", DFNT_CHAR8, &
                    len_trim(Data_Description_String),trim(Data_Description_String)) + Istatus

  !------------------------------------------------------------------
  !  determine size of input fields
  !------------------------------------------------------------------
  Istatus = sfrnatt(Sd_Id_Input,sffattr(Sd_Id_Input,"NUMBER_OF_LONGITUDES"),Nlon_In) + Istatus
  Istatus = sfrnatt(Sd_Id_Input,sffattr(Sd_Id_Input,"NUMBER_OF_LATITUDES"),Nlat_In) + Istatus
  Istatus = sfrnatt(Sd_Id_Input,sffattr(Sd_Id_Input,"LONGITUDE_SPACING"),Dlon_In) + Istatus
  Istatus = sfrnatt(Sd_Id_Input,sffattr(Sd_Id_Input,"LATITUDE_SPACING"),Dlat_In) + Istatus
  Istatus = sfrnatt(Sd_Id_Input,sffattr(Sd_Id_Input,"WESTERN_MOST_LONGITUDE"),Lon_West_In) + Istatus
  Istatus = sfrnatt(Sd_Id_Input,sffattr(Sd_Id_Input,"EASTERN_MOST_LONGITUDE"),Lon_East_In) + Istatus
  Istatus = sfrnatt(Sd_Id_Input,sffattr(Sd_Id_Input,"NORTHERN_MOST_LATITUDE"),Lat_North_In) + Istatus
  Istatus = sfrnatt(Sd_Id_Input,sffattr(Sd_Id_Input,"SOUTHERN_MOST_LATITUDE"),Lat_South_In) + Istatus

  !---- prevent round-off
  Lon_West_In = nint(Lon_West_In*100)/100.0
  Lon_East_In = nint(Lon_East_In*100)/100.0
  Lat_North_In = nint(Lat_North_In*100)/100.0
  Lat_South_In = nint(Lat_South_In*100)/100.0


  !----------------------------------------------------------------------------------
  !  Constrain
  !----------------------------------------------------------------------------------
  Lat_South_Out = max(Lat_South_Out,Lat_South_In)
  Lat_North_Out = min(Lat_North_Out,Lat_North_In)
  Lon_West_Out = min(Lon_West_Out,Lon_West_In)
  Lon_East_Out = max(Lon_East_Out,Lon_East_In)

  !--- allocate input arrays
  allocate(Scaled_Sds_Data_Input(Nlon_In,Nlat_In))
  allocate(Scaled_Sds_Latitude_Input(Nlat_In))
  allocate(Scaled_Sds_Longitude_Input(Nlon_In))

  !------------------------------------------------------------------
  ! determine stride and length of output data
  !------------------------------------------------------------------
  Lon_Stride = nint(Dlon_Out / Dlon_In)
  Lat_Stride = nint(DLat_Out / Dlat_In)
  Lon_Start = (Lon_West_Out - Lon_West_In) / Dlon_In + 1
  Lat_Start = (Lat_South_Out - Lat_South_In) / Dlat_In + 1

  Nlon_Out = int((Lon_East_Out - Lon_West_Out) / Dlon_Out)  + 1
  NLat_Out = int((Lat_North_Out - Lat_South_Out) / DLat_Out) + 1

  !-----------------------------------------------------------------
  ! allocate arrays for output
  !-----------------------------------------------------------------
  allocate(Scaled_Sds_Data_Output(Nlon_Out,NLat_Out))
  allocate(Scaled_Sds_Latitude_Output(NLat_Out))
  allocate(Scaled_Sds_Longitude_Output(Nlon_Out))

  !------------------------------------------------------------------
  ! replace attributes that differ from the input and output files 
  !------------------------------------------------------------------
  Istatus = sfsnatt(Sd_Id_Output,"NUMBER_OF_LONGITUDES",DFNT_INT32,1,Nlon_Out)
  Istatus = sfsnatt(Sd_Id_Output,"NUMBER_OF_LATITUDES",DFNT_INT32,1,NLat_Out)
  Istatus = sfsnatt(Sd_Id_Output,"LONGITUDE_SPACING",DFNT_FLOAT32,1,Dlon_Out)
  Istatus = sfsnatt(Sd_Id_Output,"LATITUDE_SPACING",DFNT_FLOAT32,1,DLat_Out)
  Istatus = sfsnatt(Sd_Id_Output,"WESTERN_MOST_LONGITUDE",DFNT_FLOAT32,1,Lon_West_Out)
  Istatus = sfsnatt(Sd_Id_Output,"EASTERN_MOST_LONGITUDE",DFNT_FLOAT32,1,Lon_East_Out)
  Istatus = sfsnatt(Sd_Id_Output,"NORTHERN_MOST_LATITUDE",DFNT_FLOAT32,1,Lat_North_Out)
  Istatus = sfsnatt(Sd_Id_Output,"SOUTHERN_MOST_LATITUDE",DFNT_FLOAT32,1,Lat_South_Out)

  Sds_Loop: do Isds = 0, Num_Sds_Input-1

  Sds_Input%Id_Input = sfselect(Sd_Id_Input,Isds) !note, Sds index starts at zero

  Istatus = sfginfo(Sds_Input%Id_Input,  &
                    Sds_Input%Variable_Name, &
                    Sds_Input%Rank, &
                    Sds_Dims, &
                    Sds_Input%Data_Type, &
                    Sds_Input%Num_Attrs) + Istatus

      !--- determine if this sds name is included in list of sds in output
      Ifound = 0
      Sds_Loop_Inner: do Isds_Inner = 1, Num_Sds_Output
         if (trim(Sds_Input%Variable_Name) == trim(Sds_Output_Names(Isds_Inner))) then
           Ifound = 1
           exit
         endif
      enddo Sds_Loop_Inner   

      if (Ifound == 0) cycle

       !--- read sds from input file
       if (Sds_Input%Variable_Name == "longitude") then
           call READ_SDS(Sd_Id_Input,Scaled_Sds_Longitude_Input,Sds_Input)
           Sds_Output = Sds_Input
           Sds_Output%Id_Output = Sd_Id_Output
           call DEFINE_SDS_RANK1(Sd_Id_Output, (/Nlon_Out/), (/Nlon_Out/), Sds_Output)
           Scaled_Sds_Longitude_Input = nint(Scaled_Sds_Longitude_Input*100)/100.0
       else if (Sds_Input%Variable_Name == "latitude") then
           call READ_SDS(Sd_Id_Input,Scaled_Sds_Latitude_Input,Sds_Input)
           Sds_Output = Sds_Input
           Sds_Output%Id_Output = Sd_Id_Output
           call DEFINE_SDS_RANK1(Sd_Id_Output, (/NLat_Out/), (/NLat_Out/), Sds_Output)
           Scaled_Sds_Latitude_Input = nint(Scaled_Sds_Latitude_Input*100)/100.0
       else 
           call READ_SDS(Sd_Id_Input,Scaled_Sds_Data_Input,Sds_Input)
           Sds_Output = Sds_Input
           Sds_Output%Id_Output = Sd_Id_Output
           call DEFINE_SDS_RANK2(Sd_Id_Output, (/1,Nlon_Out,NLat_Out/), (/1,Nlon_Out,NLat_Out/), Sds_Output)
       endif
  
  !--- crop the data to match the desired output size and write to the file
  if (Sds_Input%Variable_Name == "longitude") then

     Scaled_Sds_Longitude_Output = Scaled_Sds_Longitude_Input(Lon_Start:Lon_Start+Nlon_Out-1:Lon_Stride) 
     call WRITE_SDS(Scaled_Sds_Longitude_Output,Sds_Output)

  elseif (Sds_Input%Variable_Name == "latitude") then

     Scaled_Sds_Latitude_Output = Scaled_Sds_Latitude_Input(Lat_Start:Lat_Start+NLat_Out-1:Lat_Stride) 
     call WRITE_SDS(Scaled_Sds_Latitude_Output,Sds_Output)

  else

     Scaled_Sds_Data_Output = Scaled_Sds_Data_Input(Lon_Start:Lon_Start+Nlon_Out-1:Lon_Stride,Lat_Start:Lat_Start+NLat_Out-1:Lat_Stride) 
     call WRITE_SDS(Scaled_Sds_Data_Output,Sds_Output)

  endif

  end do Sds_Loop

  !--- deallocate memory
  if (allocated(Scaled_Sds_Data_Input)) deallocate(Scaled_Sds_Data_Input)
  if (allocated(Scaled_Sds_Latitude_Input)) deallocate(Scaled_Sds_Latitude_Input)
  if (allocated(Scaled_Sds_Longitude_Input)) deallocate(Scaled_Sds_Longitude_Input)
  if (allocated(Scaled_Sds_Data_Output)) deallocate(Scaled_Sds_Data_Output)
  if (allocated(Scaled_Sds_Latitude_Output)) deallocate(Scaled_Sds_Latitude_Output)
  if (allocated(Scaled_Sds_Longitude_Output)) deallocate(Scaled_Sds_Longitude_Output)

  !----- close input hdf file
  Istatus = sfend(Sd_Id_Input) + Istatus

  !----- close output hdf file
  Istatus = sfend(Sd_Id_Output) + Istatus

  end subroutine SUBSET_LEVEL2b
!----------------------------------------------------------------------
! An example of how to reset a random seed based on system time
!
! taken from:
! http://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html
!----------------------------------------------------------------------
  subroutine init_random_seed()
            integer :: i, n, clock
            integer, DIMENSION(:), ALLOCATABLE :: seed
          
            call RANDOM_SEED(size = n)
            allocate(seed(n))
          
            call SYSTEM_CLOCK(COUNT=clock)
          
            seed = clock + 37 * (/ (i - 1, i = 1, n) /)
            call RANDOM_SEED(PUT = seed)
          
            deallocate(seed)
  end subroutine init_random_seed
!--------------------------------------------------------------------
! File Search - similar to IDL
!--------------------------------------------------------------------
  subroutine FILE_SEARCH(Path,Search_String,Num_Files,Files)
        character(len=*), intent(in):: Path
        character(len=*), intent(in):: Search_String
        character(len=1020), dimension(:), allocatable, intent(out):: Files
        integer(kind=int4), intent(out):: Num_Files
        integer(kind=int4):: Lun
        integer(kind=int4):: Idx

        Lun = GET_LUN()

        !-- first determine number of files
!       call system('ls -1 -p '//trim(Path)//' | grep level2b.hdf | wc -l > subset_file_list')
        call system('ls -1 -p '//trim(Path)//' | grep '//trim(Search_String)//' | wc -l > subset_file_list')

        !-- next, list the files
!       call system('ls -1 -p '//trim(Path)//' | grep level2b.hdf >> subset_file_list')
        call system('ls -1 -p '//trim(Path)//' | grep '//trim(Search_String)//' >> subset_file_list')

        open (unit=Lun, file = 'subset_file_list', status="old",action="read")

        read(unit=Lun,fmt=*) Num_Files

        if (Num_Files == 0) then
           print *, "FILE_SEARCH found no files with ",trim(Search_String)," in  ", trim(Path)
           return
        endif

        allocate(Files(Num_Files))

        do Idx = 1, Num_Files
          read(unit=Lun,fmt=*) Files(Idx)
        enddo 

        close(unit=Lun) 
        
  end subroutine
!---------------------------------------------------------------------------------
!
!---------------------------------------------------------------------------------
subroutine COMPUTE_WMO_ID_KNOWING_SENSOR_NAME(Sensor_Name,WMO_Id)

character(len=*), intent(in):: Sensor_Name
integer(kind=int2), intent(out):: WMO_Id

select case (trim(Sensor_Name))


case("zen")
  WMO_Id = 0
case("metop-02")
  WMO_Id = 4
case("metop-01")
  WMO_Id = 3
case("metop-03")
  WMO_Id = 5
case("met8")
  WMO_Id = 55
case("met9")
  WMO_Id = 56
case("met10")
  WMO_Id = 57
case("met11")
  WMO_Id = 70
case("mtsat-1r")
  WMO_Id = 171 
case("mtsat-2")
  WMO_Id = 172 
case("noaa-06")
  WMO_Id = 706
case("noaa-07")
  WMO_Id = 707
case("noaa-05")
  WMO_Id = 708
case("noaa-08")
  WMO_Id = 200
case("noaa-09")
  WMO_Id = 201
case("noaa-10")
  WMO_Id = 202
case("noaa-11")
  WMO_Id = 203
case("noaa-12")
  WMO_Id = 204
case("noaa-14")
  WMO_Id = 205
case("noaa-15")
  WMO_Id = 206
case("noaa-16")
  WMO_Id = 207
case("noaa-17")
  WMO_Id = 208
case("noaa-18")
  WMO_Id = 209
case("noaa-19")
  WMO_Id = 223
case("npp")
  WMO_Id = 224
case("iff_npp")
  WMO_Id = 224
case("npp-nasa")
  WMO_Id = 224
case("goes-08")
  WMO_Id = 252
case("goes-09")
  WMO_Id = 253
case("goes-10")
  WMO_Id = 254
case("goes-11")
  WMO_Id = 255
case("goes-12")
  WMO_Id = 256
case("goes-13")
  WMO_Id = 257
case("goes-14")
  WMO_Id = 258
case("goes-15")
  WMO_Id = 259
case("terra-modis")
  WMO_Id = 783
case("aqua-modis")
  WMO_Id = 784
case("coms")
  WMO_Id = 810

case default

print *, "Unknown Sensor Name, Can not assign Sc_Id"
stop

end select

end subroutine COMPUTE_WMO_ID_KNOWING_SENSOR_NAME

!--------------------------------------------------------------------
!  Write the flag_values and flag_meanings attributes for
!  relevant variables
!--------------------------------------------------------------------
subroutine WRITE_FLAG_ATTRIBUTES(Sds_Id,Sds_Name,IStatus_Sum)
   integer, intent(in):: Sds_Id
   integer, intent(inout):: Istatus_Sum
   character(len=*), intent(in):: Sds_Name
   character(len=500):: Temp_Name
   integer(kind=int1), dimension(2):: yes_no_flag_values
   integer(kind=int1), dimension(4):: qf_flag_values

   integer:: sfsnatt
   integer:: sfscatt

   yes_no_flag_values = int((/0,1/),kind=int1)
   qf_flag_values = int((/0,1,2,3/),kind=int1)

     if (Sds_Name == "acha_quality") then
       Temp_Name = "Processed "// &
                   "valid_Tc_retrieval "// &
                   "valid_ec_retrieval "// &
                   "valid_beta_retrieval "// &
                   "degraded_Tc_retrieval "// &
                   "degraded_ec_retrieval "// &
                   "degraded_beta_retrieval "
       Istatus_Sum = sfsnatt(Sds_Id, "flag_masks", DFNT_INT8, 7, int((/1,2,4,8,16,32,64/),kind=int1)) + Istatus_Sum
       Istatus_Sum = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, len_trim(Temp_Name), trim(Temp_Name)) + Istatus_Sum
     end if

     if (Sds_Name == "acha_info") then
       Temp_Name = "Cloud_Height_Attempted "// &
                   "Bias_Correction_Employed "// &
                   "Ice_Cloud_Retrieval "// &
                   "Local_Radiative_Center_Processing_Used "// &
                   "Multi_Layer_Retrieval "// &
                   "Lower_Cloud_Interpolation_Used "// &
                   "Boundary_Layer_Inversion_Assumed "
       Istatus_Sum = sfsnatt(Sds_Id, "flag_masks", DFNT_INT8, 7, int((/1,2,4,8,16,32,64/),kind=int1)) + Istatus_Sum
       Istatus_Sum = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, len_trim(Temp_Name), trim(Temp_Name)) + Istatus_Sum
     end if

     if (Sds_Name == "dcomp_quality") then
       Temp_Name = "Processed "// &
                   "valid_COD_retrieval "// &
                   "valid_REF_retrieval "// &
                   "degraded_COD_retrieval "// &
                   "degraded_REF_retrieval "// &
                   "convergency "// &
                   "glint "
       Istatus_Sum = sfsnatt(Sds_Id, "flag_masks", DFNT_INT8, 7, int((/1,2,4,8,16,32,64/),kind=int1)) + Istatus_Sum
       Istatus_Sum = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, len_trim(Temp_Name), trim(Temp_Name)) + Istatus_Sum
     end if

     if (Sds_Name == "dcomp_info") then
       Temp_Name = "land_or_sea_mask "// &
                   "day_or_night_mask "// &
                   "twilight_@65-82_solar_zenith "// &
                   "snow "// &
                   "sea_ice "// &
                   "liquid_or_ice_phase "// &
                   "thick_cloud "// &
                   "thin_cloud "
       Istatus_Sum = sfsnatt(Sds_Id, "flag_masks", DFNT_INT16, 8, int((/1,2,4,8,16,32,64,128/),kind=int2)) + Istatus_Sum
       Istatus_Sum = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, len_trim(Temp_Name), trim(Temp_Name)) + Istatus_Sum
     end if

     if (Sds_Name == "cld_opd_dcomp_qf") then
       Istatus_Sum = sfsnatt(Sds_Id, "flag_values", DFNT_INT8, 4, qf_flag_values) + Istatus_Sum
       Istatus_Sum = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, 49, "not_attempted "// &
                                " failed "// &
                                " low_quality "// &
                                " high_quality ") + Istatus_Sum
     end if

     if (Sds_Name == "cld_reff_dcomp_qf") then
       Istatus_Sum = sfsnatt(Sds_Id, "flag_values", DFNT_INT8, 4, qf_flag_values) + Istatus_Sum
       Istatus_Sum = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, 49, "not_attempted "// &
                                " failed "// &
                                " low_quality "// &
                                " high_quality ") + Istatus_Sum
     end if

     if (Sds_Name == "cld_temp_acha_qf") then
       Istatus_Sum = sfsnatt(Sds_Id, "flag_values", DFNT_INT8, 4, qf_flag_values) + Istatus_Sum
       Istatus_Sum = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, 49, "not_attempted "// &
                                " failed "// &
                                " low_quality "// &
                                " high_quality ") + Istatus_Sum
     end if

     if (Sds_Name == "cld_emiss_acha_qf") then
       Istatus_Sum = sfsnatt(Sds_Id, "flag_values", DFNT_INT8, 4, qf_flag_values) + Istatus_Sum
       Istatus_Sum = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, 49, "not_attempted "// &
                                " failed "// &
                                " low_quality "// &
                                " high_quality ") + Istatus_Sum
     end if

     if (Sds_Name == "cld_beta_acha_qf") then
       Istatus_Sum = sfsnatt(Sds_Id, "flag_values", DFNT_INT8, 4, qf_flag_values) + Istatus_Sum
       Istatus_Sum = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, 49, "not_attempted "// &
                                " failed "// &
                                " low_quality "// &
                                " high_quality ") + Istatus_Sum
     end if

     if (Sds_Name == "cloud_mask") then
       Istatus_Sum = sfsnatt(Sds_Id, "flag_values", DFNT_INT8, 4, qf_flag_values) + Istatus_Sum
       Istatus_Sum = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, 47, "clear "// &
                                " probably_clear "// &
                                " probably_cloudy "// &
                                " cloudy ") + Istatus_Sum
     end if

     if (Sds_Name == "cloud_type") then
       Temp_Name = "clear "// &
                   "probably_clear "// &
                   "fog "// &
                   "water "// &
                   "supercooled_water "// &
                   "mixed "// &
                   "opaque_ice "// &
                   "cirrus "// &
                   "overlapping "// &
                   "overshooting "// &
                   "unknown "// &
                   "dust "// &
                   "smoke "
       Istatus_Sum = sfsnatt(Sds_Id, "flag_values", DFNT_INT8, 13, int((/0,1,2,3,4,5,6,7,8,9,10,11,12/),kind=int1)) + Istatus_Sum
       Istatus_Sum = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, len_trim(Temp_Name), trim(Temp_Name)) + Istatus_Sum
     end if

     if (Sds_Name == "surface_type") then
       Temp_Name = "water "// &
                   "evergreen_needle "// &
                   "evergreen_broad "// &
                   "deciduous_needle "// &
                   "deciduous_broad "// &
                   "mixed_forest "// &
                   "woodlands "// &
                   "wooded_grass "// &
                   "closed_shrubs "// &
                   "open_shrubs "// &
                   "grasses "// &
                   "croplands "// &
                   "bare "// &
                   "urban "
       Istatus_Sum = sfsnatt(Sds_Id, "flag_values", DFNT_INT8, 14, int((/0,1,2,3,4,5,6,7,8,9,10,11,12,13/),kind=int1)) + Istatus_Sum
       Istatus_Sum = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, len_trim(Temp_Name), trim(Temp_Name)) + Istatus_Sum
     end if

     if (Sds_Name == "land_class") then
       Istatus_Sum = sfsnatt(Sds_Id, "flag_values", DFNT_INT8, 8, int((/0,1,2,3,4,5,6,7/),kind=int1)) + Istatus_Sum
       Istatus_Sum = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, 109, "ocean "// &
                                " land "// &
                                " coastline "// &
                                " shallow_inland_water "// &
                                " ephemeral_water "// &
                                " deep_inland_water "// &
                                " moderate_ocean "// &
                                " deep_ocean ") + Istatus_Sum
     end if

     if (Sds_Name == "snow_class") then
       Istatus_Sum = sfsnatt(Sds_Id, "flag_values", DFNT_INT8, 3, int((/1,2,3/),kind=int1)) + Istatus_Sum
       Istatus_Sum = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, 30, "no_snow_or_ice "// &
                                " sea_ice "// &
                                " snow ") + Istatus_Sum
     end if

end subroutine WRITE_FLAG_ATTRIBUTES

!------------------------------------------------------------------------------------------------
!  write the sds attributes that communicate how data is scaled
!------------------------------------------------------------------------------------------------
subroutine WRITE_SCALING_ATTRIBUTES(Sds,Istatus_Sum)
    type(Sds_Struct), intent(in):: Sds
    integer, intent(inout):: Istatus_Sum
    integer:: sfsnatt


     !--- determined if a scaled Sds, if so write needed attributes for scaling
     if (Sds%Scaling_Type > 0) then

      Istatus_Sum = sfsnatt(Sds%Id_Output, "unscaled_missing", DFNT_FLOAT32, 1, Sds%Unscaled_Missing) + Istatus_Sum
      Istatus_Sum = sfsnatt(Sds%Id_Output, "actual_range", DFNT_FLOAT32, 2, Sds%Actual_Range) + Istatus_Sum
      if (Sds%Data_Type == DFNT_INT8) then
         Istatus_Sum = sfsnatt(Sds%Id_Output, "valid_range", Sds%Data_Type, 2, int(Sds%valid_range,kind=int1)) + Istatus_Sum
      elseif (Sds%Data_Type == DFNT_INT16) then
         Istatus_Sum = sfsnatt(Sds%Id_Output, "valid_range", Sds%Data_Type, 2, int(Sds%valid_range,kind=int2)) + Istatus_Sum
      elseif (Sds%Data_Type == DFNT_FLOAT32) then
         Istatus_Sum = sfsnatt(Sds%Id_Output, "valid_range", Sds%Data_Type, 2, Sds%valid_range) + Istatus_Sum
     endif



      !--- write remaining attributes
      if (Sds%Scaling_Type == 1) then
       Istatus_Sum = sfsnatt(Sds%Id_Output, "scale_factor", DFNT_FLOAT32, 1, Sds%Scale_Factor) + Istatus_Sum
       Istatus_Sum = sfsnatt(Sds%Id_Output, "add_offset", DFNT_FLOAT32, 1, Sds%Add_Offset) + Istatus_Sum
!     else  !deprecated
!      Istatus_Sum = sfsnatt(Sds%Id_Output, "RANGE_MIN", DFNT_FLOAT32, 1, Sds%Unscaled_Min) + Istatus_Sum
!      Istatus_Sum = sfsnatt(Sds%Id_Output, "RANGE_MAX", DFNT_FLOAT32, 1, Sds%Unscaled_Max) + Istatus_Sum
!      Istatus_Sum = sfsnatt(Sds%Id_Output, "SCALED_MIN", DFNT_INT32, 1, Sds%Scaled_Min) + Istatus_Sum
!      Istatus_Sum = sfsnatt(Sds%Id_Output, "SCALED_MAX", DFNT_INT32, 1, Sds%Scaled_Max) + Istatus_Sum
!      Istatus_Sum = sfsnatt(Sds%Id_Output, "SCALED_MISSING", DFNT_INT32, 1, Sds%Scaled_Missing) + Istatus_Sum
      endif

     endif

     !--- write out fill value for all but packed variables
     if (Sds%Unscaled_Missing /= No_Attribute_Missing_Value) then
         if (Sds%Data_Type == DFNT_INT8) then
           Istatus_Sum = sfsnatt(Sds%Id_Output, "_FillValue", Sds%Data_Type, 1, int(Sds%Fill_Value,kind=int1)) + Istatus_Sum
         elseif (Sds%Data_Type == DFNT_INT16) then
           Istatus_Sum = sfsnatt(Sds%Id_Output, "_FillValue", Sds%Data_Type, 1, int(Sds%Fill_Value,kind=int2)) + Istatus_Sum
         elseif (Sds%Data_Type == DFNT_FLOAT32) then
           Istatus_Sum = sfsnatt(Sds%Id_Output, "_FillValue", Sds%Data_Type, 1, Sds%Fill_Value) + Istatus_Sum
         endif
     endif

end subroutine WRITE_SCALING_ATTRIBUTES
!------------------------------------------------------------------------------------------------
!  write the sds attributes that communicate how data is scaled
!------------------------------------------------------------------------------------------------
subroutine READ_SCALING_ATTRIBUTES(Sds,Istatus)
    type(Sds_Struct), intent(inout):: Sds
    integer, intent(inout):: Istatus

    integer(kind=int1):: fill_value_i1
    integer(kind=int2):: fill_value_i2
    integer(kind=int4):: fill_value_i4
    real(kind=real4):: fill_value_r4
    integer(kind=int1), dimension(2):: valid_range_i1
    integer(kind=int2), dimension(2):: valid_range_i2
    integer(kind=int4), dimension(2):: valid_range_i4
    real(kind=real4), dimension(2):: valid_range_r4

    integer:: sffattr, sfrnatt, sfrcatt

    !--- read Sds attributes
    Istatus = sfrnatt(Sds%Id_Input, sffattr(Sds%Id_Input,"SCALED"), Sds%Scaling_Type)
    Istatus = sfrcatt(Sds%Id_Input, sffattr(Sds%Id_Input,"units"), Sds%Units)
    if (Istatus /= 0) Sds%Units = "none"
    Istatus = sfrcatt(Sds%Id_Input, sffattr(Sds%Id_Input,"standard_name"), Sds%Standard_Name)
    if (Istatus /= 0) Sds%Standard_Name = "none"
    Istatus = sfrcatt(Sds%Id_Input, sffattr(Sds%Id_Input,"long_name"), Sds%Long_Name)
    if (Istatus /= 0) Sds%Long_Name = "none"


    if (Sds%Scaling_Type /= sym%NO_SCALING) then
       Istatus = sfrnatt(Sds%Id_Input, sffattr(Sds%Id_Input,"unscaled_missing"), Sds%Unscaled_Missing)
       Istatus = sfrnatt(Sds%Id_Input, sffattr(Sds%Id_Input,"actual_range"), Sds%Actual_Range)
    endif

    if (Sds%Scaling_Type /= sym%NO_SCALING) then
       if (Sds%Data_Type == DFNT_INT8) then
         Istatus = sfrnatt(Sds%Id_Input, sffattr(Sds%Id_Input,"valid_range"), valid_range_i1)
         Sds%Valid_Range = int(valid_range_i1,kind=int4)
       elseif (Sds%Data_Type == DFNT_INT16) then
         Istatus = sfrnatt(Sds%Id_Input, sffattr(Sds%Id_Input,"valid_range"), valid_range_i2)
         Sds%Valid_Range = int(valid_range_i2,kind=int4)
       elseif (Sds%Data_Type == DFNT_INT32) then
         Istatus = sfrnatt(Sds%Id_Input, sffattr(Sds%Id_Input,"valid_range"), valid_range_i4)
         Sds%Valid_Range = int(valid_range_i4,kind=int4)
       elseif (Sds%Data_Type == DFNT_FLOAT32) then
         Istatus = sfrnatt(Sds%Id_Input, sffattr(Sds%Id_Input,"valid_range"), valid_range_r4)
         Sds%Valid_Range = int(valid_range_r4,kind=real4)
       endif
    endif

    !--- fill value (if absent, assume this is packed variable and set
    !--- Unscaled_Missing to -888
    if (Sds%Data_Type == DFNT_INT8) then
         Istatus = sfrnatt(Sds%Id_Input, sffattr(Sds%Id_Input,"_FillValue"), fill_value_i1)
         if (Istatus == 0) then
           Sds%Fill_Value = int(fill_value_i1,kind=int4)
         else
           Sds%Fill_Value = Missing_Value_int4
           Sds%Unscaled_Missing = No_Attribute_Missing_Value
         endif
    elseif (Sds%Data_Type == DFNT_INT16) then
         Istatus = sfrnatt(Sds%Id_Input, sffattr(Sds%Id_Input,"_FillValue"), fill_value_i2)
         if (Istatus == 0) then
            Sds%Fill_Value = int(fill_value_i2,kind=int4)
         else
           Sds%Fill_Value = Missing_Value_int4
           Sds%Unscaled_Missing = No_Attribute_Missing_Value
         endif
    elseif (Sds%Data_Type == DFNT_INT32) then
         Istatus = sfrnatt(Sds%Id_Input, sffattr(Sds%Id_Input,"_FillValue"), fill_value_i4)
         if (Istatus == 0) then
            Sds%Fill_Value = int(fill_value_i4,kind=int4)
         else
           Sds%Fill_Value = Missing_Value_int4
           Sds%Unscaled_Missing = No_Attribute_Missing_Value
         endif
    elseif (Sds%Data_Type == DFNT_FLOAT32) then
         Istatus = sfrnatt(Sds%Id_Input, sffattr(Sds%Id_Input,"_FillValue"), fill_value_r4)
         if (Istatus == 0) then
            Sds%Fill_Value = int(fill_value_r4,kind=int4)
         else
           Sds%Fill_Value = Missing_Value_int4
           Sds%Unscaled_Missing = No_Attribute_Missing_Value
         endif
    endif

    !-- if scaled, read attributes that allow unscaling
    if (Sds%Scaling_Type > 0) then
      if (Sds%Scaling_Type == sym%LINEAR_SCALING) then
       Istatus = sfrnatt(Sds%Id_Input, sffattr(Sds%Id_Input,"scale_factor"), Sds%Scale_Factor) + Istatus
       Istatus = sfrnatt(Sds%Id_Input, sffattr(Sds%Id_Input,"add_offset"), Sds%Add_Offset) + Istatus
!     else
!      Istatus = sfrnatt(Sds%Id_Input, sffattr(Sds%Id_Input,"SCALED_MISSING"), Sds%Scaled_Missing) + Istatus
!      Istatus = sfrnatt(Sds%Id_Input, sffattr(Sds%Id_Input,"SCALED_MIN"), Sds%Scaled_Min) + Istatus
!      Istatus = sfrnatt(Sds%Id_Input, sffattr(Sds%Id_Input,"SCALED_MAX"), Sds%Scaled_Max) + Istatus
!      Istatus = sfrnatt(Sds%Id_Input, sffattr(Sds%Id_Input,"RANGE_MIN"), Sds%Unscaled_Min) + Istatus
!      Istatus = sfrnatt(Sds%Id_Input, sffattr(Sds%Id_Input,"RANGE_MAX"), Sds%Unscaled_Max) + Istatus
      endif
    endif

end subroutine READ_SCALING_ATTRIBUTES


end module LEVEL2B_ROUTINES
