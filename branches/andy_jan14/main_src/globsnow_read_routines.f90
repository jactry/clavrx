! $Id$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: globsnow_read_routines.f90 (src)
!       GLOBSNOW_READ_ROUTINES (program)
!
! PURPOSE: Routines for opening, reading and closing the GlobSnow snow coverage files
!
! DESCRIPTION:
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
! Includes: GET_GlobSnow_FILENAME
!           READ_GlobSnow_ANALYSIS_MAP
!           GET_PIXEL_GlobSnow_ANALYSIS 
!
!--------------------------------------------------------------------------------------
module GLOBSNOW_READ_ROUTINES
   use HDF
   use CONSTANTS
   use NUMERICAL_ROUTINES
   use FILE_UTILITY
   implicit none
   private
   public:: GET_PIXEL_GLOBSNOW_ANALYSIS
   public:: GET_GLOBSNOW_FILENAME
   public:: READ_GLOBSNOW_ANALYSIS_MAP
   integer, parameter, public:: Num_Lon_GlobSnow = 721
   integer, parameter, public:: Num_Lat_GlobSnow = 721 
   real, parameter, private:: GLOBSNOW_SWE_THRESHOLD = 0.001
   real, parameter, private:: Min_Lat_GlobSnow = 35.0
   real, parameter, private:: Max_Lat_GlobSnow = 85.0
   real(kind=real4), dimension(Num_Lon_GlobSnow,Num_Lat_GlobSnow), public, save:: GlobSnow_Map

   integer(kind=int4), parameter, private :: MAX_GLOBSNOW_LATENCY = 3 !including current day

contains

!----------------------------------------------------------------------------------------
!
! Name: GET_GlobSnow_FILENAME
!
! Description: Function to find the GlobSnow map filename.
!
!  Input: year_in - year being processed
!         day_of_year - day being processed
!         GlobSnow_path - path to globsnow files
!
!  Output: GlobSnow_filename - GlobSnow path and filename
!
!----------------------------------------------------------------------------------------

 FUNCTION GET_GlobSnow_FILENAME(year_in,day_of_year,GlobSnow_path) result(GlobSnow_filename)
   CHARACTER(*), intent(in) :: GlobSnow_path
   integer(kind=int2), intent(in):: year_in
   integer(kind=int2), intent(in):: day_of_year
   integer:: GlobSnow_option
   CHARACTER(len=256) :: GlobSnow_filename,GlobSnow_full
   CHARACTER(len=256) :: GlobSnow_filename_tmp
   integer(kind=int4) :: iday, year, month, day, jday, ileap
   CHARACTER (len=2)   :: day_string, month_string
   CHARACTER (len=4)   :: year_string

   GlobSnow_filename = "no_file"

   GlobSnow_option = 1

   if (GlobSnow_option == 1) then

     do iday=0, MAX_GLOBSNOW_LATENCY - 1
       jday = day_of_year - iday
       year = year_in 
       ileap = leap_year_fct(year)
       if (jday < 1) then
         year = year - 1
         ileap = leap_year_fct(year)
         jday = (365 + ileap) + jday
       endif 
       month = compute_month(jday, ileap)
       day = compute_day(jday, ileap)
       write (year_string,fmt="(I4)") year
       write (month_string, '(I2.2)') month
       write (day_string,   '(I2.2)') day

       GlobSnow_filename_tmp= trim(year_string)//"/" &
                 //"GlobSnow_SWE_L3A_"//year_string//month_string//day_string//"_v1.0.hdf"
       GlobSnow_full = trim(GlobSnow_path)//trim(GlobSnow_filename_tmp)
       
       if (file_exists(trim(GlobSnow_full)) .eqv. .true.) then
         GlobSnow_filename = GlobSnow_filename_tmp
         print *, EXE_PROMPT, "Found ", trim(GlobSnow_filename_tmp)
         exit
       endif
     end do
   endif
   return

 END FUNCTION GET_GlobSnow_FILENAME

!----------------------------------------------------------------------------------------
! 
! Name: READ_GlobSnow_ANALYSIS_MAP
!
! Description: Opens and reads in globsnow hdf file and converts snow water equivalent
!              into a snow mask using a threshold of 0.001 mm
!
!  Input: GlobSnow_name - globsnow path and filename
!
!  Output: GlobSnow_Map - snow cover mask map
!
!----------------------------------------------------------------------------------------
 subroutine READ_GLOBSNOW_ANALYSIS_MAP(GlobSnow_name)
   character(len=*), intent(in):: GlobSnow_name
   character(len=128):: file_temp
   character(len=100) :: sds_name
   integer(kind=int4) :: sd_id
   integer(kind=int4) :: astatus, istatus, sds_id, sds_rank, sds_type, sds_nattr
   integer(kind=int4), dimension(2) :: sds_dims
   integer(kind=int4):: sfend, sfstart, sfselect, sfrdata, sfendacc, &
                        sfginfo, sfn2index
   real(kind=real4), dimension(:,:), allocatable :: buffer_map
   
   
   integer, dimension(2):: start, stride, edge 
     

   !--- uncompress
 !  system_string = "gunzip -c "//trim(GlobSnow_name)// &
 !       " > temp_GlobSnow_file"
  ! call system(system_string)
 
  ! file_temp = "temp_GlobSnow_file"
  
   file_temp = GlobSnow_name

!-------------------------------------------------------------------------------------
! read in data
!-------------------------------------------------------------------------------------
   istatus = 0
   start(1)=0
   start(2)=0
   stride(1)=1
   stride(2)=1
   edge(1) = Num_Lon_GlobSnow
   edge(2) = Num_Lat_GlobSnow
   sds_name = "swe"

   sd_id = sfstart(GlobSnow_name, DFACC_RDONLY)

   sds_id = sfselect(sd_id, sfn2index(sd_id,trim(sds_name)))    ! 0 is the index for SWE

   istatus = sfginfo(sds_id, sds_name, sds_rank, sds_dims, sds_type, sds_nattr) + istatus
   if (istatus /= 0) then
     print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
     stop
   else
!    print "(a,'Successful read from ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
   endif
 
!   if (edge(1) < 0) edge(1) = sds_dims(1)
!   if (edge(2) < 0) edge(2) = sds_dims(2)
 
!   edge(1) = min((sds_dims(1) - start(1)),edge(1))
!   edge(2) = min((sds_dims(2) - start(2)),edge(2))
 
   allocate(buffer_map(edge(1),edge(2)),stat=astatus)
   if (astatus /= 0) then
     print "(a,'Not enough memory to allocate 2d buffer.')",EXE_PROMPT
     stop
   endif
 
   !--- read from this sds
   istatus = sfrdata(sds_id, start, stride, edge, buffer_map) + istatus

   !--- end access to this sds
   istatus = sfendacc(sds_id) + istatus

   !-- close hdf file
   istatus = sfend(sd_id) + istatus

   if (istatus /= 0) then
     print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
     stop
   else
!    print "(a,'Successful Read of GlobSnow from ',a)",EXE_PROMPT,trim(GlobSnow_name)
   endif
   
   GlobSnow_Map = buffer_map

 end subroutine READ_GLOBSNOW_ANALYSIS_MAP

!-----------------------------------------------------------------------------------------
!
! Name: GET_PIXEL_GLOBSNOW_ANALYSIS
!
! Description: Locate GlobSnow map value for each AVHRR pixel. The GlobSnow data is on an
!              Equal-Area Scalable Earth Grid (EASE-Grid) and provides the whole Northern 
!              Hemisphere (lambert's equal-area azimuthal projection) in a single data 
!              field. The nominal resolution of a single pixel is 25 km x 25 km
!
!  Input: Latitude - pixel Latitudes of the segment
!         Longitude - pixel Longitudes of the segment
!         Land_Class - land type classification
!         Invalid_Pixel_Mask - yes/no mask for bad pixels
!
!  Output: GlobSnow snow mask value
!
!-----------------------------------------------------------------------------------------
 subroutine GET_PIXEL_GLOBSNOW_ANALYSIS(Latitude,Longitude,Land_Class,Invalid_Pixel_Mask,Snow_Out)
  real(kind=real4), dimension(:,:), intent(in) :: Latitude
  real(kind=real4), dimension(:,:), intent(in) :: Longitude
  integer(kind=int1), dimension(:,:), intent(in)::Land_Class
  integer(kind=int1), dimension(:,:), intent(in)::Invalid_Pixel_Mask
  integer(kind=int1), dimension(:,:), intent(out)::Snow_Out
  integer:: ielem
  integer:: iline
  real:: xLon, xLat  
! radius of the earth (km), authalic sphere based on International datum
  real,parameter:: Re_km = 6371.228
! nominal cell size in kilometers
  real,parameter ::Cell_km = 25.067525
! scale factor for standard paralles at +/-30.00 degrees
  real,parameter ::COS_PHI1 = .866025403
  integer:: col, row
  real Rg, phi, lam, rho
  integer::c0,r0
  integer ::  nx
  integer ::  ny
  integer ::  nx_map
  integer ::  ny_map
  
  nx = size(Snow_Out,1)
  ny = size(Snow_Out,2)
  nx_map = Num_Lon_GlobSnow
  ny_map = Num_Lat_GlobSnow 

  !-- compute derived constants here outside of loop
  Rg = Re_km/Cell_km
  c0 = (Num_Lon_GlobSnow-1)/2.0
  r0 = (Num_Lat_GlobSnow-1)/2.0
  
  !---------------------------------------------------------------------
  ! loop of pixels in this segment
  !---------------------------------------------------------------------

  line_loop:  do iline = 1 , ny
    element_loop:    do ielem = 1, nx

      !--- skip bad pixels
      if (Invalid_Pixel_Mask(ielem,iline) == sym%YES) then
         cycle
      endif

      if (Longitude(ielem,iline) < 0.0) then
         xLon = Longitude(ielem,iline) + 360.0
      else
         xLon = Longitude(ielem,iline)
      endif
      xLat = Latitude(ielem,iline)

      if ((xLat < Min_Lat_GlobSnow) .or. (xLat > Max_Lat_GlobSnow)) then   ! GlobSnow data is only available from 35N-85N
        cycle
      endif

      !----------------------------------------------------------------------------------------
      ! convert lat/lon coordinates to EASE-grid row/column coordinates
      !
      !  Input : latitude - local variable Latitude 
      !          longitude - local variable Longitude
      !
      !  Output: col - column coordinate for EASE-grid
      !          row - row coordinate for EASE-grid
      !----------------------------------------------------------------------------------------
      phi = dtor * xLat
      lam = dtor * xLon

      rho = 2 * Rg * sin(pi/4.0 - phi/2.0)

      col = min(nx_map,max(1,nint(c0 + rho * sin(lam))))
      row = min(ny_map,max(1,nint(r0 + rho * cos(lam))))

      !---------------------------------------------------------------------------------------
      ! default to no snow
      !---------------------------------------------------------------------------------------
      Snow_Out(ielem,iline) = sym%NO_SNOW

      !---------------------------------------------------------------------------------------
      ! convert from snow water equivalent (mm) to snow mask (> 0.001)
      !---------------------------------------------------------------------------------------
      if (GlobSnow_Map(col,row) >= GLOBSNOW_SWE_THRESHOLD) then
        Snow_Out(ielem,iline) = sym%SNOW
      endif

      !---------------------------------------------------------------------------------------
      ! convert to sea-ice over water
      !---------------------------------------------------------------------------------------
      if ((Snow_Out(ielem,iline) == sym%SNOW) .and. (Land_Class(ielem,iline) /= sym%LAND)) then
        Snow_Out(ielem,iline) = sym%SEA_ICE
      endif

    end do element_loop
  end do line_loop

 end subroutine GET_PIXEL_GLOBSNOW_ANALYSIS

end module GLOBSNOW_READ_ROUTINES
