! $Id$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: ncep_reanalysis.f90 (src)
!       NCEP_REANALYSIS (program)
!
! PURPOSE: This module houses all of the routines necessary to interface with
!          NCEP Reanalysis Data
!
! DESCRIPTION: 
!            Note this is hardcoded for the current 2.5x2.5 degree data.  It checks
!            to make sure that this is the case, if not, it reports this and stops
!
!            this restriction comes from the mapping of the T62 gaussian fields to
!            to the 2.5x2.5 fields.  This step uses a nearest neighbor approach and 
!            should be revisted or this step should be moved outside of CLAVR-x.
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
! File I/O: Two netcdf files are opened and closed used HDF routines
!
! Public Routines
!  READ_NCEP_REANALYSIS_DATA - main routine to read in the fields from the
!                             NCEP Reanalysis
! 
! Private Routines
!  READ_DATA_1D - routines to read in one-dimensional fields
!  READ_DATA_2D - routines to read in two-dimensional fields
!  READ_DATA_3D - routines to read in three-dimensional fields
!  LOCATE_OUR_TIME - given a time, locate the correct NCEP Reanalysis fields
!--------------------------------------------------------------------------------------
module NCEP_REANALYSIS

use CONSTANTS
use NWP_COMMON
use SORT_MODULE
use NUMERICAL_TOOLS_MOD

implicit none
include 'hdf.f90'
private:: READ_DATA_1D, READ_DATA_2D, READ_DATA_3D, LOCATE_OUR_TIME
public:: READ_NCEP_REANALYSIS_DATA

real (kind=real4), parameter, private :: missing_ncep = -9.999E+20

integer, parameter, private:: nlon_ncep = 144, nlat_ncep = 73

integer, dimension(nlon_ncep), parameter, private:: &
  lon_index_gauss = (/ &
   1,   2,   4,   5,   6,   8,   9,  10,  12,  13, &
  14,  16,  17,  18,  20,  21,  22,  24,  25,  26, &
  28,  29,  30,  32,  33,  34,  36,  37,  38,  40, &
  41,  42,  44,  45,  46,  48,  49,  50,  52,  53, &
  54,  56,  57,  58,  60,  61,  62,  64,  65,  66, &
  68,  69,  70,  72,  73,  74,  76,  77,  78,  80, &
  81,  82,  84,  85,  86,  88,  89,  90,  92,  93, &
  94,  96,  97,  98, 100, 101, 102, 104, 105, 106, &
 108, 109, 110, 112, 113, 114, 116, 117, 118, 120, &
 121, 122, 124, 125, 126, 128, 129, 130, 132, 133, &
 134, 136, 137, 138, 140, 141, 142, 144, 145, 146, &
 148, 149, 150, 152, 153, 154, 156, 157, 158, 160, &
 161, 162, 164, 165, 166, 168, 169, 170, 172, 173, &
 174, 176, 177, 178, 180, 181, 182, 184, 185, 186, &
 188, 189, 190, 192 /)

integer, dimension(nlat_ncep), parameter, private::&
lat_index_gauss =  (/ &
   1,   2,   3,   4,   5,   7,   8,   9,  11,  12, &
  13,  15,  16,  17,  19,  20,  21,  23,  24,  25, &
  26,  28,  29,  30,  32,  33,  34,  36,  37,  38, &
  40,  41,  42,  44,  45,  46,  47,  49,  50,  51, &
  53,  54,  55,  57,  58,  59,  61,  62,  63,  65, &
  66,  67,  69,  70,  71,  72,  74,  75,  76,  78, &
  79,  80,  82,  83,  84,  86,  87,  88,  90,  91, &
  92,  93,  94 /)

contains

subroutine READ_NCEP_REANALYSIS_DATA(start_year, start_jday, start_itime, &
                                     end_jday, end_itime, data_dir)

! Input/output arguments                                     
character(len=*), intent(in) :: data_dir
integer(kind=int2), intent(in) :: start_year, start_jday, end_jday
integer(kind=int4), intent(in) :: start_itime, end_itime

! Other vars
integer :: n, nlevels
character(len=4) :: year_string
real(kind=real4), dimension(:,:,:), allocatable :: ncep_airtemp_profile
real(kind=real4), dimension(:,:,:), allocatable :: ncep_geopothgt
real(kind=real4), dimension(:,:,:), allocatable :: ncep_relhum_profile
real(kind=real4), dimension(:,:), allocatable :: ncep_psfc
real(kind=real4), dimension(:,:), allocatable :: ncep_airsfc
real(kind=real4), dimension(:,:), allocatable :: ncep_tpw
real(kind=real4), dimension(:,:), allocatable :: ncep_rhsfc
real(kind=real4), dimension(:,:), allocatable :: ncep_zsfc
real(kind=real4), dimension(:,:), allocatable :: ncep_skint
real(kind=real4), dimension(:,:), allocatable :: ncep_skint_before
real(kind=real4), dimension(:,:), allocatable :: ncep_skint_after
real(kind=real4), dimension(:,:), allocatable :: ncep_weasd
real(kind=real4), dimension(:,:), allocatable :: ncep_u_wnd_10m
real(kind=real4), dimension(:,:), allocatable :: ncep_v_wnd_10m
real(kind=real4), dimension(:,:), allocatable :: ncep_t_trop
real(kind=real4), dimension(:,:), allocatable :: ncep_p_trop
real(kind=real4), dimension(:,:), allocatable :: ncep_land
real(kind=real4), dimension(:), allocatable :: ncep_airtemp_level
real(kind=real4), dimension(:), allocatable :: ncep_relhum_level
integer :: levels(50) = -1
integer :: i
integer :: j
integer :: ii
integer :: jj
integer(kind=int4) :: bot
integer(kind=int4) :: top
integer(kind=int4) :: level
real (kind=real4), allocatable, dimension(:) :: lat
real (kind=real4), allocatable, dimension(:) :: lon
real, parameter :: R_wv = 461.5
real, parameter :: g = 9.8

real::level1b_file_start_time
real::level1b_file_end_time
real::level1b_file_mean_time
integer:: mean_doy
real:: mean_hours
real:: x_interp


!---- make year string for filename
write (year_string, '(I4.4)') start_year


!--- determine ncep field indices that bound the mean orbit time
!--- year starts with bot=1. Read routines account for zero index
level1b_file_start_time = start_jday + start_itime/86400000.0_real4 
level1b_file_end_time = end_jday + end_itime/86400000.0_real4 
level1b_file_mean_time = 0.5*(level1b_file_start_time + level1b_file_end_time)
mean_doy = int(level1b_file_mean_time)
mean_hours = (level1b_file_mean_time - mean_doy)*24.0
bot = (mean_doy-1)*4 + int(mean_hours/6.0) + 1
top = bot + 1
x_interp = (mean_hours - 6.0*int(mean_hours/6.0))/6.0
x_interp = max(0.0,min(1.0,x_interp))

!-- these are public variables used for temporal interpolation
nwp_start_hour = 6*int(mean_hours/6.0)
nwp_end_hour = nwp_start_hour + 6

!print *, "level1b times = ", level1b_file_start_time, level1b_file_end_time, level1b_file_mean_time, mean_doy
!print *, "indices = ", bot, top
!print *, "x_interp = ", x_interp
!print *, "start, end, mean hours = ", nwp_start_hour, nwp_end_hour, mean_hours

!
! Read 1D data
!

call READ_DATA_1D(trim(data_dir)//"air."//year_string//".nc", &
   "lat", lat)
call READ_DATA_1D(trim(data_dir)//"air."//year_string//".nc", &
   "lon", lon)
call READ_DATA_1D(trim(data_dir)//"air."//year_string//".nc", &
   "level", ncep_airtemp_level)

!
! Reading 3D data
!

! air temperature
call READ_DATA_3D(trim(data_dir)//"air."//year_string//".nc", &
   "air", ncep_airtemp_profile, bot, top, x_interp)

! geopotential height
call READ_DATA_3D(trim(data_dir)//"hgt."//year_string//".nc", &
   "hgt", ncep_geopothgt, bot, top, x_interp)

! relative humidity
call READ_DATA_3D(trim(data_dir)//"rhum."//year_string//".nc", &
   "rhum", ncep_relhum_profile, bot, top, x_interp)

!
! Reading 2D data
!

! surface air temperature
call READ_DATA_2D( &
   trim(data_dir)//"air.sig995."//year_string//".nc", &
   "air", ncep_airsfc, bot, top, x_interp)

! surface skin temperature - note this is a gaussian grid
call READ_DATA_2D( &
    trim(data_dir)//"skt.sfc.gauss."//year_string//".nc", &
    "skt", ncep_skint, bot, top, x_interp)
call READ_DATA_2D( &
   trim(data_dir)//"skt.sfc.gauss."//year_string//".nc", &
   "skt", ncep_skint_before, bot, top, 0.0)
call READ_DATA_2D( &
   trim(data_dir)//"skt.sfc.gauss."//year_string//".nc", &
   "skt", ncep_skint_after, bot, top, 1.0)

! water equiv snow depth - note this is a gaussian grid
call READ_DATA_2D( &
   trim(data_dir)//"weasd.sfc.gauss."//year_string//".nc", &
   "weasd", ncep_weasd, bot, top, x_interp)


! u wind at 10m - note this is a gaussian grid
call READ_DATA_2D( &
   trim(data_dir)//"uwnd.10m.gauss."//year_string//".nc", &
   "uwnd", ncep_u_wnd_10m, bot, top, x_interp)

! v wind at 10m - note this is a gaussian grid
call READ_DATA_2D( &
   trim(data_dir)//"vwnd.10m.gauss."//year_string//".nc", &
   "vwnd", ncep_v_wnd_10m, bot, top, x_interp)


! precipitable water
call READ_DATA_2D( &
   trim(data_dir)//"pr_wtr.eatm."//year_string//".nc", &
   "pr_wtr", ncep_tpw, bot, top, x_interp)

! surface pressure
call READ_DATA_2D( &
   trim(data_dir)//"pres.sfc."//year_string//".nc", "pres", &
   ncep_psfc, bot, top, x_interp)

! surface rh
call READ_DATA_2D( &
   trim(data_dir)//"rhum.sig995."//year_string//".nc", "rhum", &
   ncep_rhsfc, bot, top, x_interp)

! Surface geopotential height
call READ_UNI_DATA_2D(trim(data_dir)//"hgt.sfc.nc", "hgt", &
   ncep_zsfc)

! Surface NCEP land mask
call READ_UNI_DATA_2D(trim(data_dir)//"land.nc", "land", ncep_land)

! tropopause temperature

call READ_DATA_2D( &
   trim(data_dir)//"air.tropp."//year_string//".nc", "air", &
   ncep_t_trop, bot, top, x_interp)

! tropopause pressure
call READ_DATA_2D( &
   trim(data_dir)//"pres.tropp."//year_string//".nc", "pres", &
   ncep_p_trop, bot, top, x_interp)

!----- convert surface pressure from Pa to hPa
ncep_psfc = ncep_psfc / 100.0
ncep_p_trop = ncep_p_trop / 100.0

!------
missing_nwp = missing_ncep
nlat_nwp = size(lat, 1)
nlon_nwp = size(lon, 1)

!--- check to see x,y dimensions are as expected
 if (nlat_nwp /= nlat_ncep) then
     print *, "latitude dimensions differ from expectations for NCEP reanalysis data, stopping"
     stop
 endif
 if (nlon_nwp /= nlon_ncep) then
     print *, "longitude dimensions differ from expectations for NCEP reanalysis data, stopping"
     stop
 endif

!--- make other dimensions
npoints = nlon_nwp*nlat_nwp
nlevels = size(ncep_airtemp_level, 1)
nlevels_nwp = nlevels
lat1_nwp = lat(1)
lon1_nwp = lon(1)
dlat_nwp = lat(2) - lat(1)
dlon_nwp = lon(2) - lon(1)

call CREATE_NWP_ARRAYS()

levels(1:nlevels) = int(ncep_airtemp_level)
call sort(levels(1:nlevels))
do n = 1, nlevels
   P_std_nwp(n) = real(levels(n))
enddo

!--- there is no ozone profile in ncep reanalysis
ozone_prof_nwp = 0.0

!------------------------------------------------------------------
! reorder NCEP arrays
!-----------------------------------------------------------------
!--- reorder
do level = 1, nlevels
  z_prof_nwp(level,:,:) = ncep_geopothgt(:,:,nlevels+1-level)
  t_prof_nwp(level,:,:) = ncep_airtemp_profile(:,:,nlevels+1-level)
enddo

n = size(ncep_relhum_profile,3)
!rh_prof_nwp(:,:,nlevels-n+1:nlevels) = ncep_relhum_profile(:,:,n:1:-1)
do level = nlevels-n+1,nlevels
 rh_prof_nwp(level,:,:) = ncep_relhum_profile(:,:,nlevels-level+1)
enddo

!------------------------------------------------------------
! fill in 2d sfc arrays
!------------------------------------------------------------
psfc_nwp = ncep_psfc
tmpair_nwp = ncep_airsfc
tpw_nwp = ncep_tpw / 10.0         !now in cm
zsfc_nwp = 0.0    ! this is missing
rhsfc_nwp = ncep_rhsfc
land_nwp = int(ncep_land,kind=int1)
t_trop_nwp = ncep_t_trop
p_trop_nwp = ncep_p_trop

!--- Convert those things on the T62 Gaussian grid to the 2.5x2.5 grid
  do j = 1, nlat_nwp
    do i = 1, nlon_nwp
      tmpsfc_nwp_before(i,j) = ncep_skint_before(lon_index_gauss(i),lat_index_gauss(j))
      tmpsfc_nwp_after(i,j) = ncep_skint_after(lon_index_gauss(i),lat_index_gauss(j))
      tmpsfc_nwp(i,j) = ncep_skint(lon_index_gauss(i),lat_index_gauss(j))
      weasd_nwp(i,j) = ncep_weasd(lon_index_gauss(i),lat_index_gauss(j))
      u_wnd_10m_nwp(i,j) = ncep_u_wnd_10m(lon_index_gauss(i),lat_index_gauss(j))
      v_wnd_10m_nwp(i,j) = ncep_v_wnd_10m(lon_index_gauss(i),lat_index_gauss(j))
    enddo
  enddo



!--- store 500 mb heights
hght500_nwp = z_prof_nwp(12,:,:)  


!---- compute wind speed 
wnd_spd_10m_nwp = missing_value_real4
where((u_wnd_10m_nwp /= missing_ncep) .and. (v_wnd_10m_nwp /= missing_ncep)) 
 wnd_spd_10m_nwp = sqrt(u_wnd_10m_nwp**2 + v_wnd_10m_nwp**2)
endwhere

!---- compute wind direction 
wnd_dir_10m_nwp = missing_value_real4
where((wnd_spd_10m_nwp >= 0.0) .and. (v_wnd_10m_nwp /= missing_ncep))
 wnd_dir_10m_nwp = acos(-1.0*v_wnd_10m_nwp/wnd_spd_10m_nwp) / dtor     !in degrees  (0 to 180)
endwhere
where((u_wnd_10m_nwp >= 0.0 ) .and. (u_wnd_10m_nwp /= missing_ncep))
 wnd_dir_10m_nwp = 360.0 - wnd_dir_10m_nwp     
endwhere

!---- compute the 3x3 uniformity of tmpair
do i = 1, nlon_nwp
 ii = max(2,min(nlon_nwp-1,i))
 do j = 1, nlat_nwp
  jj = max(2,min(nlat_nwp-1,j))
  tmpair_uni_nwp(i,j) = (maxval(tmpair_nwp(ii-1:ii+1,jj-1:jj+1)) - &
                        minval(tmpair_nwp(ii-1:ii+1,jj-1:jj+1)) ) / 3.0
  tmpsfc_uni_nwp(i,j) = (maxval(tmpsfc_nwp(ii-1:ii+1,jj-1:jj+1)) - &
                        minval(tmpsfc_nwp(ii-1:ii+1,jj-1:jj+1)) ) / 3.0
 enddo
enddo

! convert z_prof_nwp to km
where (z_prof_nwp /= missing_nwp) z_prof_nwp= z_prof_nwp/1000.0_real4


!--- deallocate temp arrays
if (allocated(ncep_airtemp_level)) deallocate(ncep_airtemp_level)
if (allocated(ncep_relhum_level)) deallocate(ncep_relhum_level)
if (allocated(ncep_psfc)) deallocate(ncep_psfc)
if (allocated(ncep_t_trop)) deallocate(ncep_t_trop)
if (allocated(ncep_p_trop)) deallocate(ncep_p_trop)
if (allocated(ncep_airsfc)) deallocate(ncep_airsfc)
if (allocated(ncep_rhsfc)) deallocate(ncep_rhsfc)
if (allocated(ncep_tpw)) deallocate(ncep_tpw)
if (allocated(ncep_airtemp_profile)) deallocate(ncep_airtemp_profile)
if (allocated(ncep_geopothgt)) deallocate(ncep_geopothgt)
if (allocated(ncep_relhum_profile)) deallocate(ncep_relhum_profile)
if (allocated(lat)) deallocate(lat)
if (allocated(lon)) deallocate(lon)
if (allocated(ncep_skint)) deallocate(ncep_skint)
if (allocated(ncep_skint_before)) deallocate(ncep_skint_before)
if (allocated(ncep_skint_after)) deallocate(ncep_skint_after)
if (allocated(ncep_weasd)) deallocate(ncep_weasd)
if (allocated(ncep_u_wnd_10m)) deallocate(ncep_u_wnd_10m)
if (allocated(ncep_v_wnd_10m)) deallocate(ncep_v_wnd_10m)
if (allocated(ncep_zsfc)) deallocate(ncep_zsfc)
if (allocated(ncep_land)) deallocate(ncep_land)

end subroutine READ_NCEP_REANALYSIS_DATA


!
! Private subroutines
!

subroutine LOCATE_OUR_TIME(ncep_file, our_time, bot, top, time_before, &
   time_after)

! Input/output args
character(len=*), intent(in) :: ncep_file
real(kind=real8), intent(in) :: our_time
integer(kind=int4), intent(out) :: bot, top
real(kind=real8), intent(out) :: time_before, time_after

! Other vars.
character(len=4) :: var_name="time"
real(kind=real8), allocatable, dimension(:) :: time
integer:: istatus, sd_id, sds_idx, sds_id, sds_rank, sds_type, &
          sds_nattr
integer, dimension(1):: sds_dims, start, edges, stride
character(len=128):: sds_name
integer(kind=int4) :: mid

! HDF API declaration
integer:: sfstart, sfend, sfselect, sfn2index, sfginfo, sfrdata, sfendacc



!print "(a)", "reading 1D variable '"//trim(var_name)//"' from "//trim(ncep_file)
          
sd_id = sfstart(ncep_file, DFACC_READ)
if (sd_id == FAIL) then
   print *, "failed open on "//trim(ncep_file)//" file"
   stop 261
endif

start = 0
stride = 1
istatus = 0

sds_idx = sfn2index(sd_id, var_name)
if (sds_idx == FAIL) then
   print *, "SDS '"//trim(var_name)//"' not found"
   stop 262
endif

sds_id = sfselect(sd_id, sds_idx)
if (sds_id == FAIL) then
   print *, "sfselect() on '"//trim(var_name)//"' failed"
   stop 263
endif

istatus = sfginfo(sds_id, sds_name, sds_rank, sds_dims, sds_type, sds_nattr)
if (istatus == FAIL) then
   print *, "sfginfo() on '"//trim(var_name)//"' failed"
   stop 264
endif

if (sds_rank /= 1) then
   print *, "Rank of SDS '"//trim(var_name)//"' not 1"
   stop 265
endif

allocate( time(sds_dims(1)) )

edges(1) = sds_dims(1)
istatus = sfrdata(sds_id, start, stride, edges, time)
if (istatus == FAIL) then
   print *, "failed reading data from SDS '"//trim(var_name)//"'"
   stop 266
endif

istatus = sfendacc(sds_id)
if (istatus == FAIL) then
   print *, "sfendacc() on '"//trim(var_name)//"' failed"
   stop 267
endif

istatus = sfend(sd_id)
if (istatus == FAIL) then
   print "(a)", "sfend() failed on "//trim(ncep_file)
   stop 268
endif

bot = lbound(time, 1)
top = ubound(time, 1)
if ( our_time < time(bot) ) then
   if ( (time(bot)-our_time) > 6.0 ) then
      print *, "orbit's time more than 6 hours before data time range"
      stop "298-before"
   else
      print *, "orbit's time just outside of data time range; extrapolation"
      top = bot+1
   endif
else if ( our_time > time(top) ) then
   if ( (our_time-time(top)) > 6.0 ) then
      print *, "orbit's time more than 6 hours after data time range"
      stop "298-after"
   else
      print *, "orbit's time just outside of data time range; extrapolation"
      bot = top-1
   endif
else
   do
      if (top-bot == 1) exit
      mid = bot+(top-bot)/2  ! <-- this must be integer division
      if (our_time < time(mid)) then
         top = mid
      else
         bot = mid
      endif
   end do
endif

time_before = time(bot)
time_after = time(top)

deallocate(time)

end subroutine LOCATE_OUR_TIME

subroutine READ_DATA_1D(ncep_file, var_name, var)

! Input/output args
character(len=*), intent(in) :: ncep_file, var_name
real(kind=real4), intent(out), allocatable, dimension(:) :: var

integer:: istatus, sd_id, sds_idx, sds_id, sds_rank, sds_type, &
          sds_nattr
integer, dimension(1):: sds_dims, start, edges, stride
character(len=128):: sds_name

! HDF API declaration
integer:: sfstart, sfend, sfselect, sfn2index, sfginfo, sfrdata, sfendacc

          

!print "(a)", "reading 1D variable '"//trim(var_name)//"' from "//trim(ncep_file)
          
sd_id = sfstart(ncep_file, DFACC_READ)
if (sd_id == FAIL) then
   print *, "failed open on "//trim(ncep_file)//" file"
   stop 251
endif

start = 0
stride = 1
istatus = 0

sds_idx = sfn2index(sd_id, var_name)
if (sds_idx == FAIL) then
   print *, "SDS '"//trim(var_name)//"' not found"
   stop 252
endif

sds_id = sfselect(sd_id, sds_idx)
if (sds_id == FAIL) then
   print *, "sfselect() on '"//trim(var_name)//"' failed"
   stop 253
endif

istatus = sfginfo(sds_id, sds_name, sds_rank, sds_dims, sds_type, sds_nattr)
if (istatus == FAIL) then
   print *, "sfginfo() on '"//trim(var_name)//"' failed"
   stop 254
endif

if (sds_rank /= 1) then
   print *, "Rank of SDS '"//trim(var_name)//"' not 1"
   stop 255
endif

allocate( var(sds_dims(1)) )

edges(1) = sds_dims(1)
istatus = sfrdata(sds_id, start, stride, edges, var)
if (istatus == FAIL) then
   print *, "failed reading data from SDS '"//trim(var_name)//"'"
   stop 256
endif

istatus = sfendacc(sds_id)
if (istatus == FAIL) then
   print *, "sfendacc() on '"//trim(var_name)//"' failed"
   stop 257
endif

istatus = sfend(sd_id)
if (istatus == FAIL) then
   print "(a)", "sfend() failed on "//trim(ncep_file)
   stop 258
endif

end subroutine READ_DATA_1D

!-----------------------------------------------------------------------
! This reads in a time-invariant 2d array - such as surface altitude
!-----------------------------------------------------------------------
subroutine READ_UNI_DATA_2D(ncep_file, var_name, var)

! Input/output arguments                                     
character(len=*), intent(in) :: ncep_file, var_name
real(kind=real4), intent(out), allocatable, dimension(:,:) :: var

! Other variables
integer:: istatus, sd_id, sds_idx, sds_id, sds_rank, sds_type, &
          sds_nattr, attr_idx
integer, dimension(3):: sds_dims, start, edges, stride
character(len=128):: sds_name
integer(kind=int2), allocatable, dimension(:,:) :: temp_var
real(kind=real4), dimension(1) :: add_offset, scale_factor
integer(kind=int2), dimension(1) :: missing_value

integer:: sfstart, sfend, sfselect, sfn2index, sfginfo, sfrdata, sfendacc, &
          sffattr, sfrnatt


!print "(a)", "reading 2D variable '"//trim(var_name)//"' from "//trim(ncep_file)

sd_id = sfstart(ncep_file, DFACC_READ)
if (sd_id == FAIL) then
   print *, "failed open on "//trim(ncep_file)//" file"
   stop 230
endif

start = 0
stride = 1
istatus = 0

!
! Read in the data
!
istatus = 0
sds_idx = sfn2index(sd_id, var_name)
if (sds_idx == FAIL) then
   print *, "SDS '"//trim(var_name)//"' not found"
   stop 238
endif

sds_id = sfselect(sd_id, sds_idx)
if (sds_id == FAIL) then
   print *, "sfselect() on '"//trim(var_name)//"' failed"
   stop 239
endif

istatus = sfginfo(sds_id, sds_name, sds_rank, sds_dims, sds_type, sds_nattr)
if (istatus == FAIL) then
   print *, "sfginfo() on '"//trim(var_name)//"' failed"
   stop 240
endif

if (sds_rank /= 3) then
   print *, "Rank of SDS '"//trim(var_name)//"' not 3"
   stop 241
endif

allocate( temp_var(sds_dims(1), sds_dims(2)) )
allocate( var(sds_dims(1), sds_dims(2)) )

edges = sds_dims
edges(3) = 1
start = 0
istatus = sfrdata(sds_id, start, stride, edges, temp_var)
if (istatus == FAIL) then
   print *, "failed reading data from SDS '"//trim(var_name)//"'"
   stop 242
endif

! Read in now SDS's scaling attributes: add_offset, scale_factor, and
! missing_value.
istatus = 0
attr_idx = sffattr(sds_id, "add_offset")
if (attr_idx == FAIL) then
   print *, "sffattr() failed on 'add_offset'"
   stop 243
endif
istatus = sfrnatt(sds_id, attr_idx, add_offset)
if (istatus == FAIL) then
   print *, "sfrnatt() failed on 'add_offset'"
   stop 244
endif

attr_idx = sffattr(sds_id, "scale_factor")
if (attr_idx == FAIL) then
   print *, "sffattr() failed on 'scale_factor'"
   stop 245
endif
istatus = sfrnatt(sds_id, attr_idx, scale_factor)
if (istatus == FAIL) then
   print *, "sfrnatt() failed on 'scale_factor'"
   stop 246
endif

attr_idx = sffattr(sds_id, "missing_value")
if (attr_idx == FAIL) then
   print *, "sffattr() failed on 'missing_value'"
   stop 247
endif
istatus = sfrnatt(sds_id, attr_idx, missing_value)
if (istatus == FAIL) then
   print *, "sfrnatt() failed on 'missing_value'"
   stop 248
endif

! End access to the SDS and the file
istatus = 0
istatus = sfendacc(sds_id)
if (istatus == FAIL) then
   print *, "sfendacc() on '"//trim(var_name)//"' failed"
   stop 250
endif
istatus = sfend(sd_id)
if (istatus == FAIL) then
   print "(a)", "sfend() failed on "//trim(ncep_file)
   stop "250-sfend"
endif

var = missing_ncep
where (temp_var /= missing_value(1))
   var = scale_factor(1)*temp_var+add_offset(1)
end where
deallocate(temp_var)

!print *, "minval("//trim(var_name)//") = ", minval(var)
!print *, "maxval("//trim(var_name)//") = ", maxval(var)

end subroutine READ_UNI_DATA_2D

!-----------------------------------------------------------------------
! Read two-dimensional fields
!
!  ncep_file = name of the file
!  var_name = name of the variable
!  bot = index time in the file of the time_before
!  top = index time in the file of the time_after
!  time_before = the time of the before field
!  time_after = the time of the after field
!-----------------------------------------------------------------------
subroutine READ_DATA_2D(ncep_file, var_name, var, bot, top, x_interp)

! Input/output arguments                                     
character(len=*), intent(in) :: ncep_file, var_name
integer(kind=int4), intent(in) :: bot, top
real(kind=real4), intent(in) :: x_interp
real(kind=real4), intent(out), allocatable, dimension(:,:) :: var

! Other variables
integer:: istatus, sd_id, sds_idx, sds_id, sds_rank, sds_type, &
          sds_nattr, attr_idx
integer, dimension(3):: sds_dims, start, edges, stride
character(len=128):: sds_name
integer(kind=int2), allocatable, dimension(:,:) :: temp_var
real(kind=real4), allocatable, dimension(:,:) :: before, after    !originally coded as integer - bug?
real(kind=real4), dimension(1) :: add_offset, scale_factor
integer(kind=int2), dimension(1) :: missing_value

integer:: sfstart, sfend, sfselect, sfn2index, sfginfo, sfrdata, sfendacc, &
          sffattr, sfrnatt

!print "(a)", "reading 2D variable '"//trim(var_name)//"' from "//trim(ncep_file)

sd_id = sfstart(ncep_file, DFACC_READ)
if (sd_id == FAIL) then
   print *, "failed open on "//trim(ncep_file)//" file"
   stop 230
endif

start = 0
stride = 1
istatus = 0

!
! Read in the slices of data defined by 'bot' and 'top' indices.
!

! First the "before" slice.
istatus = 0
sds_idx = sfn2index(sd_id, var_name)
if (sds_idx == FAIL) then
   print *, "SDS '"//trim(var_name)//"' not found"
   stop 238
endif

sds_id = sfselect(sd_id, sds_idx)
if (sds_id == FAIL) then
   print *, "sfselect() on '"//trim(var_name)//"' failed"
   stop 239
endif

istatus = sfginfo(sds_id, sds_name, sds_rank, sds_dims, sds_type, sds_nattr)
if (istatus == FAIL) then
   print *, "sfginfo() on '"//trim(var_name)//"' failed"
   stop 240
endif

if (sds_rank /= 3) then
   print *, "Rank of SDS '"//trim(var_name)//"' not 3"
   stop 241
endif

allocate( temp_var(sds_dims(1), sds_dims(2)) )
allocate( before(sds_dims(1), sds_dims(2)) )

edges = sds_dims
edges(3) = 1
start = 0
start(3) = bot-1
istatus = sfrdata(sds_id, start, stride, edges, temp_var)
if (istatus == FAIL) then
   print *, "failed reading data from SDS '"//trim(var_name)//"'"
   stop 242
endif

before = temp_var
deallocate(temp_var)

! Read in now SDS's scaling attributes: add_offset, scale_factor, and
! missing_value.
istatus = 0
attr_idx = sffattr(sds_id, "add_offset")
if (attr_idx == FAIL) then
   print *, "sffattr() failed on 'add_offset'"
   stop 243
endif
istatus = sfrnatt(sds_id, attr_idx, add_offset)
if (istatus == FAIL) then
   print *, "sfrnatt() failed on 'add_offset'"
   stop 244
endif

attr_idx = sffattr(sds_id, "scale_factor")
if (attr_idx == FAIL) then
   print *, "sffattr() failed on 'scale_factor'"
   stop 245
endif
istatus = sfrnatt(sds_id, attr_idx, scale_factor)
if (istatus == FAIL) then
   print *, "sfrnatt() failed on 'scale_factor'"
   stop 246
endif

attr_idx = sffattr(sds_id, "missing_value")
if (attr_idx == FAIL) then
   print *, "sffattr() failed on 'missing_value'"
   stop 247
endif
istatus = sfrnatt(sds_id, attr_idx, missing_value)
if (istatus == FAIL) then
   print *, "sfrnatt() failed on 'missing_value'"
   stop 248
endif

! Now the "after" slice.
istatus = 0

allocate( temp_var(sds_dims(1), sds_dims(2)) )
allocate( after(sds_dims(1), sds_dims(2)) )

edges = sds_dims
edges(3) = 1
start = 0
start(3) = top-1
istatus = sfrdata(sds_id, start, stride, edges, temp_var)
if (istatus == FAIL) then
   print *, "failed reading data from SDS '"//trim(var_name)//"'"
   stop 249
endif

! End access to the SDS
istatus = 0
istatus = sfendacc(sds_id)
if (istatus == FAIL) then
   print *, "sfendacc() on '"//trim(var_name)//"' failed"
   stop 250
endif

after = temp_var
deallocate(temp_var)

! Unscale the "before" and "after" data.
where (before /= missing_value(1)) before = scale_factor(1)*before+add_offset(1)
where (after /= missing_value(1)) after = scale_factor(1)*after+add_offset(1)

! Interpolate between the two slices using a simple straight line formula.
allocate( var(sds_dims(1), sds_dims(2)) )
var = missing_ncep
where ( (before /= missing_value(1)) .and. (after /= missing_value(1)) )
   var = (1.0 - x_interp) * before + x_interp * after
end where

deallocate(before, after)

istatus = sfend(sd_id)
if (istatus == FAIL) then
   print "(a)", "sfend() failed on "//trim(ncep_file)
   stop "250-sfend"
endif

end subroutine READ_DATA_2D


subroutine READ_DATA_3D(ncep_file, var_name, var, bot, top, x_interp)

! Input/output arguments                                     
character(len=*), intent(in) :: ncep_file, var_name
integer(kind=int4), intent(in) :: bot, top
real(kind=real4), intent(in) :: x_interp
real(kind=real4), intent(out), allocatable, dimension(:,:,:) :: var

! Other variables
integer:: istatus, sd_id, sds_idx, sds_id, sds_rank, sds_type, &
          sds_nattr, attr_idx
integer, dimension(4):: sds_dims, start, edges, stride
character(len=128):: sds_name
integer(kind=int2), allocatable, dimension(:,:,:) :: temp_var
integer(kind=real4), allocatable, dimension(:,:,:) :: before, after
real(kind=real4), dimension(1) :: add_offset, scale_factor
integer(kind=int2), dimension(1) :: missing_value

integer:: sfstart, sfend, sfselect, sfn2index, sfginfo, sfrdata, sfendacc, &
          sffattr, sfrnatt

!print "(a)", "reading 3D variable '"//trim(var_name)//"' from "//trim(ncep_file)

sd_id = sfstart(ncep_file, DFACC_READ)
if (sd_id == FAIL) then
   print *, "failed open on "//trim(ncep_file)//" file"
   stop 200
endif

start = 0
stride = 1
istatus = 0

!
! Read in the slices of data defined by 'bot' and 'top' indices.
!

! First the "before" slice.
istatus = 0
sds_idx = sfn2index(sd_id, var_name)
if (sds_idx == FAIL) then
   print *, "SDS '"//trim(var_name)//"' not found"
   stop 208
endif

sds_id = sfselect(sd_id, sds_idx)
if (sds_id == FAIL) then
   print *, "sfselect() on '"//trim(var_name)//"' failed"
   stop 209
endif

istatus = sfginfo(sds_id, sds_name, sds_rank, sds_dims, sds_type, sds_nattr)
if (istatus == FAIL) then
   print *, "sfginfo() on '"//trim(var_name)//"' failed"
   stop 210
endif

if (sds_rank /= 4) then
   print *, "Rank of SDS '"//trim(var_name)//"' not 4"
   stop 211
endif

allocate( temp_var(sds_dims(1), sds_dims(2), sds_dims(3)) )
allocate( before(sds_dims(1), sds_dims(2), sds_dims(3)) )

edges = sds_dims
edges(4) = 1
start = 0
start(4) = bot-1
istatus = sfrdata(sds_id, start, stride, edges, temp_var)
if (istatus == FAIL) then
   print *, "failed reading data from SDS '"//trim(var_name)//"'"
   stop 212
endif

before = temp_var
deallocate(temp_var)

! Read in now SDS's scaling attributes: add_offset, scale_factor, and
! missing_value.
istatus = 0
attr_idx = sffattr(sds_id, "add_offset")
if (attr_idx == FAIL) then
   print *, "sffattr() failed on 'add_offset'"
   stop 213
endif
istatus = sfrnatt(sds_id, attr_idx, add_offset)
if (istatus == FAIL) then
   print *, "sfrnatt() failed on 'add_offset'"
   stop 214
endif

attr_idx = sffattr(sds_id, "scale_factor")
if (attr_idx == FAIL) then
   print *, "sffattr() failed on 'scale_factor'"
   stop 215
endif
istatus = sfrnatt(sds_id, attr_idx, scale_factor)
if (istatus == FAIL) then
   print *, "sfrnatt() failed on 'scale_factor'"
   stop 216
endif

attr_idx = sffattr(sds_id, "missing_value")
if (attr_idx == FAIL) then
   print *, "sffattr() failed on 'missing_value'"
   stop 213
endif
istatus = sfrnatt(sds_id, attr_idx, missing_value)
if (istatus == FAIL) then
   print *, "sfrnatt() failed on 'missing_value'"
   stop 214
endif

! Now the "after" slice.
istatus = 0

allocate( temp_var(sds_dims(1), sds_dims(2), sds_dims(3)) )
allocate( after(sds_dims(1), sds_dims(2), sds_dims(3)) )

edges = sds_dims
edges(4) = 1
start = 0
start(4) = top-1
istatus = sfrdata(sds_id, start, stride, edges, temp_var)
if (istatus == FAIL) then
   print *, "failed reading data from SDS '"//trim(var_name)//"'"
   stop 217
endif

! End access to the SDS
istatus = 0
istatus = sfendacc(sds_id)
if (istatus == FAIL) then
   print *, "sfendacc() on '"//trim(var_name)//"' failed"
   stop 218
endif

after = temp_var
deallocate(temp_var)

! Unscale the "before" and "after" data.
where (before /= missing_value(1)) before = scale_factor(1)*before+add_offset(1)
where (after /= missing_value(1)) after = scale_factor(1)*after+add_offset(1)

! Interpolate between the two slices using a simple straight line formula.
allocate( var(sds_dims(1), sds_dims(2), sds_dims(3)) )
var = missing_ncep
where ( (before /= missing_value(1)) .and. (after /= missing_value(1)) )
   var = (1.0 - x_interp) * before + x_interp * after
end where

deallocate(before, after)

istatus = sfend(sd_id)
if (istatus == FAIL) then
   print "(a)", "sfend() failed on "//trim(ncep_file)
   stop 219
endif

!print *, "minval("//trim(var_name)//") = ", minval(var)
!print *, "maxval("//trim(var_name)//") = ", maxval(var)

end subroutine READ_DATA_3D

end module NCEP_REANALYSIS
