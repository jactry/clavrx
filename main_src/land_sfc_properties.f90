!$Id: land_sfc_properties.f90,v 1.14.2.2 2014/01/26 04:48:34 heidinger Exp $
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: land_sfc_properties.f90 (src)
!       land_sfc_properties (program)
!
! PURPOSE: 
!
! DESCRIPTION: This module taken from GEOCAT and modified for CLAVR-x
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
!--------------------------------------------------------------------------------------

MODULE land_sfc_properties
  use HDF
  use CONSTANTS
  use NUMERICAL_ROUTINES
  use FILE_UTILITY
  
  implicit none
  
  private :: read_hdf_sds
  private :: read_hdf_global_attribute_float64
  private :: read_hdf_sds_dimenions
  public  :: read_land_sfc_hdf
  public  :: get_snow_map_filename
  
   interface read_land_sfc_hdf
        module procedure  &
           read_land_sfc_hdf_i1,  &
           read_land_sfc_hdf_i2
    end interface

   interface read_hdf_sds
        module procedure  &
           read_hdf_sds_i1,  &
           read_hdf_sds_i2
    end interface

  INTEGER(kind=int4), parameter, private :: num_lat_default = 4500
  INTEGER(kind=int4), parameter, private :: num_lon_default = 9000
  REAL(kind=real8), parameter, private :: first_lat_default = -90.0_real8
  REAL(kind=real8), parameter, private :: last_lat_default = 90.0_real8
  REAL(kind=real8), parameter, private :: first_lon_default = -180.0_real8
  REAL(kind=real8), parameter, private :: last_lon_default = 180.0_real8
  REAL(kind=real8), parameter, private :: del_lat_default = 0.04_real8
  REAL(kind=real8), parameter, private :: del_lon_default = 0.04_real8

  INTEGER(kind=int4), parameter, private :: MAX_SNOW_LATENCY = 4 !including current day

  TYPE, public :: land_grid_description
    CHARACTER(len=256) :: sds_name
    INTEGER(kind=int4) :: num_lat
    INTEGER(kind=int4) :: num_lon
    REAL(kind=real8) :: del_lat
    REAL(kind=real8) :: del_lon
    REAL(kind=real8) :: first_lat
    REAL(kind=real8) :: first_lon
  END TYPE land_grid_description
    
  CONTAINS
  
!-------------------------------------------------------------------
! Subroutine to open the land surface file.
!-------------------------------------------------------------------

FUNCTION open_land_sfc_hdf(data_dir, filename, grid_str) result(id)
  CHARACTER(len=*), intent(in) :: data_dir, filename
  TYPE(land_grid_description), optional, intent(inout) :: grid_str
  
  INTEGER(kind=int4) :: id  
  CHARACTER(len=256) :: filename_full
  
  logical :: file_exists
  
  INTEGER :: sfstart
  
  filename_full = trim(data_dir)//trim(filename)
  
  inquire(file = filename_full, exist = file_exists)
  if (.not. file_exists) then
    print "(/,a,'Land surface file, ',a,' does not exist.')",EXE_PROMPT,trim(filename_full)
    stop
  endif
  
  id = sfstart(trim(filename_full), DFACC_READ)
  if (id == FAIL) then
    print "(/,a,'Failed to open, ',a)",EXE_PROMPT,trim(filename_full)
    stop
  endif
  
  if (present(grid_str)) then
    grid_str%del_lat = read_hdf_global_attribute_float64(id, "dlat")
    grid_str%del_lon = read_hdf_global_attribute_float64(id, "dlon")
    grid_str%first_lat = read_hdf_global_attribute_float64(id, "first_lat")
    grid_str%first_lon = read_hdf_global_attribute_float64(id, "first_lon")
  
    if (grid_str%del_lat == missing_value_real8) grid_str%del_lat = del_lat_default
    if (grid_str%del_lon == missing_value_real8) grid_str%del_lon = del_lon_default
    if (grid_str%first_lat == missing_value_real8) grid_str%first_lat = first_lat_default
    if (grid_str%first_lon == missing_value_real8) grid_str%first_lon = first_lon_default
  
    call read_hdf_sds_dimenions(id, grid_str%sds_name, grid_str%num_lat, grid_str%num_lon) 

  endif
  
  return

END FUNCTION open_land_sfc_hdf

!-------------------------------------------------------------------
! Subroutine to close the land surface file.
!-------------------------------------------------------------------

SUBROUTINE close_land_sfc_hdf(id)
  INTEGER(kind=int4), intent(in) :: id
  
  INTEGER(kind=int4) :: istatus  
  INTEGER :: sfend
  
  istatus = sfend(id)
  if (istatus /= 0) then
    print "(/,a,'Error closing land surface hdf file.')",EXE_PROMPT
    stop
  endif
  

END SUBROUTINE close_land_sfc_hdf


!-------------------------------------------------------------------
! Function to find the snow map name.
!-------------------------------------------------------------------

  FUNCTION get_snow_map_filename(year_in,day_of_year,snow_path) result(snow_filename)
   CHARACTER(*), intent(in) :: snow_path
   INTEGER(kind=int2), intent(in):: year_in
   INTEGER(kind=int2), intent(in):: day_of_year
   CHARACTER(len=256) :: snow_filename
   CHARACTER(len=256) :: snow_filename_tmp
   INTEGER(kind=int4) :: iday, year, month, day, jday, ileap
   CHARACTER (len=2)   :: year_string, day_string, month_string

    snow_filename = "no_file"

            do iday=0, MAX_SNOW_LATENCY - 1
              jday = day_of_year - iday
               year = year_in 
               ileap = leap_year_fct(year)
               if (jday < 1) then
                 year = year - 1
                 jday = (365 + ileap) + jday
               endif 
               month = compute_month(jday, ileap)
               day = compute_day(jday, ileap)
			   write (year_string,  '(I2.2)') year - 100*(year/100)
			   write (month_string, '(I2.2)') month
			   write (day_string,   '(I2.2)') day
			   
			   snow_filename_tmp = "snow_map_4km_" //year_string//month_string// &
                       day_string//".hdf"
				   
                if (file_exists(trim(snow_path)//trim(snow_filename_tmp)) .eqv. .true.) then
                   snow_filename = snow_filename_tmp
                   exit
                endif
              end do
      return

 END FUNCTION get_snow_map_filename
!-------------------------------------------------------------------
! Subroutine to read in global attributes.
!-------------------------------------------------------------------

FUNCTION read_hdf_global_attribute_float64(id, attr_name) result(attr)
  INTEGER(kind=int4), intent(in) :: id
  CHARACTER(len=*), intent(in) :: attr_name
  
  REAL(kind=real8), dimension(1) :: buffer  
  REAL(kind=real8) :: attr
  INTEGER(kind=int4) :: istatus, attr_index
  INTEGER :: sffattr, sfrnatt
  
  attr_index = sffattr(id, trim(attr_name))
  istatus = sfrnatt(id, attr_index, buffer)
  if (istatus /= 0) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(attr_name),id
    buffer(1) = missing_value_real8
  endif
  
  attr = buffer(1)
  return

END FUNCTION read_hdf_global_attribute_float64

!-------------------------------------------------------------------
! Subroutine to read in SDS dimensions.
!-------------------------------------------------------------------

SUBROUTINE read_hdf_sds_dimenions(id, sds_name, num_lat, num_lon)
  INTEGER(kind=int4), intent(in) :: id
  CHARACTER(len=*), intent(in) :: sds_name
  INTEGER(kind=int4), intent(out) :: num_lat, num_lon
  
  INTEGER(kind=int4) :: sds_id, istatus, sds_rank, sds_type, sds_nattr
  INTEGER(kind=int4), dimension(2) :: sds_dims
    
  INTEGER :: sfselect, sfn2index, sfginfo, sfendacc
  
  sds_id = sfselect(id, sfn2index(id,trim(sds_name)))
    
  istatus = sfginfo(sds_id, sds_name, sds_rank, sds_dims, sds_type, sds_nattr)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),id    
    stop
  endif
 
  num_lon = sds_dims(1)
  num_lat = sds_dims(2)
  
  istatus = sfendacc(sds_id)
  
END SUBROUTINE read_hdf_sds_dimenions

!-------------------------------------------------------------------
! Subroutine to read a given land surface hdf file.
!-------------------------------------------------------------------
SUBROUTINE read_land_sfc_hdf_i1(id, grid_str, lat, lon, space_mask, land)
  INTEGER(kind=int4), intent(in) :: id
  TYPE(land_grid_description), intent(in) :: grid_str
  REAL(kind=real4), dimension(:,:), intent(in) :: lat, lon
  INTEGER(kind=int1), dimension(:,:), intent(in) :: space_mask
  INTEGER(kind=int1), dimension(:,:), intent(out) :: land
    
  INTEGER :: astatus
  INTEGER :: ilat1, ilat2, ilon1, ilon2, ilat, ilon, ilat_ad, ilon_ad, &
             ilon1_2, ilon2_2
  INTEGER :: temp, nx, ny, i, j
  INTEGER(kind=int1), dimension(:,:), allocatable :: land_grid, land_grid_2
  REAL(kind=real4) :: wlon, elon, slat, nlat
  INTEGER(kind=int1) :: dateline_flg, space_check
  
  INTEGER, dimension(2) :: start_2d, stride_2d, edge_2d, &
                           start_2d_2, stride_2d_2, edge_2d_2
  
  space_check = minval(space_mask)
  if (space_check == 1) then
    land = missing_value_int1
    return
  endif
  
  nx = size(land,1)
  ny = size(land,2)
  
  call find_bounds(lat,lon,wlon,elon,slat,nlat,dateline_flg)

  if (dateline_flg == 0) then
    
    ilat1 = max(0,min(grid_str%num_lat,int(abs(nlat - grid_str%first_lat)/grid_str%del_lat) + 0))
    ilat2 = max(0,min(grid_str%num_lat,int(abs(slat - grid_str%first_lat)/grid_str%del_lat) + 0))
  
    ilon1 = max(0,min(grid_str%num_lon,int(abs(wlon - grid_str%first_lon)/grid_str%del_lon) + 0))
    ilon2 = max(0,min(grid_str%num_lon,int(abs(elon - grid_str%first_lon)/grid_str%del_lon) + 0))
    
    if (ilat1 > ilat2) then
      temp = ilat1
      ilat1 = ilat2
      ilat2 = temp
    endif
  
    if (ilon1 > ilon2) then
      temp = ilon1
      ilon1 = ilon2
      ilon2 = temp
    endif
  
    start_2d = (/ilon1, ilat1/)
    stride_2d = (/1, 1/)
    edge_2d = (/(ilon2-ilon1)+1, (ilat2-ilat1)+1/)
  
    call read_hdf_sds(id, trim(grid_str%sds_name), start_2d, stride_2d, edge_2d, land_grid)
  
    do j = 1, ny
      do i = 1, nx
    
        if (space_mask(i,j) == sym%NO_SPACE) then
                
              ilat = max(1,min(grid_str%num_lat,int(abs(lat(i,j) - grid_str%first_lat)/grid_str%del_lat) + 1))
              ilon = max(1,min(grid_str%num_lon,int(abs(lon(i,j) - grid_str%first_lon)/grid_str%del_lon) + 1))
              ilat_ad = max(1,min((ilat - start_2d(2)) + 1,size(land_grid,2)))
              ilon_ad = max(1,min((ilon - start_2d(1)) + 1,size(land_grid,1)))
              land(i,j) = land_grid(ilon_ad,ilat_ad)
  
        endif
      
      enddo
    enddo
    
    deallocate(land_grid, stat=astatus)
    if (astatus /= 0) then
      print "(a,'Error deallocating land surface grid.')",EXE_PROMPT
      stop
    endif
    
  else
  
    ilat1 = max(0,min(grid_str%num_lat,int(abs(nlat - grid_str%first_lat)/grid_str%del_lat) + 0))
    ilat2 = max(0,min(grid_str%num_lat,int(abs(slat - grid_str%first_lat)/grid_str%del_lat) + 0))
  
    ilon1 = max(0,min(grid_str%num_lon,int(abs(wlon - grid_str%first_lon)/grid_str%del_lon) + 0))
    ilon2 = max(0,min(grid_str%num_lon,int(abs(180.0 - grid_str%first_lon)/grid_str%del_lon) + 0))
    
    ilon1_2 = max(0,min(grid_str%num_lon,int(abs(-180.0 - grid_str%first_lon)/grid_str%del_lon) + 0))
    ilon2_2 = max(0,min(grid_str%num_lon,int(abs((elon-360.0) - grid_str%first_lon)/grid_str%del_lon) + 0))
  
    if (ilat1 > ilat2) then
      temp = ilat1
      ilat1 = ilat2
      ilat2 = temp
    endif
  
    if (ilon1 > ilon2) then
      temp = ilon1
      ilon1 = ilon2
      ilon2 = temp
    endif
    
    if (ilon1_2 > ilon2_2) then
      temp = ilon1_2
      ilon1_2 = ilon2_2
      ilon2_2 = temp
    endif
  
    start_2d = (/ilon1, ilat1/)
    stride_2d = (/1, 1/)
    edge_2d = (/(ilon2-ilon1)+1, (ilat2-ilat1)+1/)
  
    call read_hdf_sds(id, trim(grid_str%sds_name), start_2d, stride_2d, edge_2d, land_grid)
    
    start_2d_2 = (/ilon1_2, ilat1/)
    stride_2d_2 = (/1, 1/)
    edge_2d_2 = (/(ilon2_2-ilon1_2)+1, (ilat2-ilat1)+1/)
  
    call read_hdf_sds(id, trim(grid_str%sds_name), start_2d_2, stride_2d_2, edge_2d_2, land_grid_2)
    
    do j = 1, ny
      do i = 1, nx
    
        if (space_mask(i,j) == sym%NO_SPACE) then

              ilat = max(1,min(grid_str%num_lat,int(abs(lat(i,j) - grid_str%first_lat)/grid_str%del_lat) + 1))
              ilon = max(1,min(grid_str%num_lon,int(abs(lon(i,j) - grid_str%first_lon)/grid_str%del_lon) + 1))

              if (lon(i,j) >= 0.0) then
                ilat_ad = max(1,min((ilat - start_2d(2)) + 1,size(land_grid,2)))
                ilon_ad = max(1,min((ilon - start_2d(1)) + 1,size(land_grid,1)))
                land(i,j) = land_grid(ilon_ad,ilat_ad)
              else
                ilat_ad = max(1,min((ilat - start_2d_2(2)) + 1,size(land_grid_2,2)))
                ilon_ad = max(1,min((ilon - start_2d_2(1)) + 1,size(land_grid_2,1)))
                land(i,j) = land_grid_2(ilon_ad,ilat_ad)
              endif

        endif
      
      enddo
    enddo
  
    deallocate(land_grid, land_grid_2, stat=astatus)
    if (astatus /= 0) then
      print "(a,'Error deallocating land surface grid.')",EXE_PROMPT
      stop
    endif
    
  endif
  
END SUBROUTINE read_land_sfc_hdf_i1

!-------------------------------------------------------------------
! Subroutine to read a given land surface hdf file.
!-------------------------------------------------------------------
SUBROUTINE read_land_sfc_hdf_i2(id, grid_str, lat, lon, space_mask, land)
  INTEGER(kind=int4), intent(in) :: id
  TYPE(land_grid_description), intent(in) :: grid_str
  REAL(kind=real4), dimension(:,:), intent(in) :: lat, lon
  INTEGER(kind=int1), dimension(:,:), intent(in) :: space_mask
  INTEGER(kind=int2), dimension(:,:), intent(out) :: land
    
  INTEGER :: astatus
  INTEGER :: ilat1, ilat2, ilon1, ilon2, ilat, ilon, ilat_ad, ilon_ad, &
             ilon1_2, ilon2_2
  INTEGER :: temp, nx, ny, i, j
  INTEGER(kind=int2), dimension(:,:), allocatable :: land_grid, land_grid_2
  REAL(kind=real4) :: wlon, elon, slat, nlat
  INTEGER(kind=int1) :: dateline_flg, space_check
  
  INTEGER, dimension(2) :: start_2d, stride_2d, edge_2d, &
                           start_2d_2, stride_2d_2, edge_2d_2
  
  space_check = minval(space_mask)
  if (space_check == 1) then
    land = missing_value_int2
    return
  endif
 

  nx = size(land,1)
  ny = size(land,2)

  call find_bounds(lat,lon,wlon,elon,slat,nlat,dateline_flg)

  if (dateline_flg == 0) then
    
    ilat1 = max(0,min(grid_str%num_lat,int(abs(nlat - grid_str%first_lat)/grid_str%del_lat) + 0))
    ilat2 = max(0,min(grid_str%num_lat,int(abs(slat - grid_str%first_lat)/grid_str%del_lat) + 0))
  
    ilon1 = max(0,min(grid_str%num_lon,int(abs(wlon - grid_str%first_lon)/grid_str%del_lon) + 0))
    ilon2 = max(0,min(grid_str%num_lon,int(abs(elon - grid_str%first_lon)/grid_str%del_lon) + 0))
    
    if (ilat1 > ilat2) then
      temp = ilat1
      ilat1 = ilat2
      ilat2 = temp
    endif
  
    if (ilon1 > ilon2) then
      temp = ilon1
      ilon1 = ilon2
      ilon2 = temp
    endif
  
    start_2d = (/ilon1, ilat1/)
    stride_2d = (/1, 1/)
    edge_2d = (/(ilon2-ilon1)+1, (ilat2-ilat1)+1/)

    call read_hdf_sds(id, trim(grid_str%sds_name), start_2d, stride_2d, edge_2d, land_grid)
  
    do j = 1, ny
      do i = 1, nx
    
        if (space_mask(i,j) == sym%NO_SPACE) then
                
              ilat = max(1,min(grid_str%num_lat,int(abs(lat(i,j) - grid_str%first_lat)/grid_str%del_lat) + 1))
              ilon = max(1,min(grid_str%num_lon,int(abs(lon(i,j) - grid_str%first_lon)/grid_str%del_lon) + 1))
              ilat_ad = max(1,min((ilat - start_2d(2)) + 1,size(land_grid,2)))
              ilon_ad = max(1,min((ilon - start_2d(1)) + 1,size(land_grid,1)))
              land(i,j) = land_grid(ilon_ad,ilat_ad)
  
        endif
      
      enddo
    enddo
    
    deallocate(land_grid, stat=astatus)
    if (astatus /= 0) then
      print "(a,'Error deallocating land surface grid.')",EXE_PROMPT
      stop
    endif
    
  else
  
    ilat1 = max(0,min(grid_str%num_lat,int(abs(nlat - grid_str%first_lat)/grid_str%del_lat) + 0))
    ilat2 = max(0,min(grid_str%num_lat,int(abs(slat - grid_str%first_lat)/grid_str%del_lat) + 0))
  
    ilon1 = max(0,min(grid_str%num_lon,int(abs(wlon - grid_str%first_lon)/grid_str%del_lon) + 0))
    ilon2 = max(0,min(grid_str%num_lon,int(abs(180.0 - grid_str%first_lon)/grid_str%del_lon) + 0))
    
    ilon1_2 = max(0,min(grid_str%num_lon,int(abs(-180.0 - grid_str%first_lon)/grid_str%del_lon) + 0))
    ilon2_2 = max(0,min(grid_str%num_lon,int(abs((elon-360.0) - grid_str%first_lon)/grid_str%del_lon) + 0))
  
    if (ilat1 > ilat2) then
      temp = ilat1
      ilat1 = ilat2
      ilat2 = temp
    endif
  
    if (ilon1 > ilon2) then
      temp = ilon1
      ilon1 = ilon2
      ilon2 = temp
    endif
    
    if (ilon1_2 > ilon2_2) then
      temp = ilon1_2
      ilon1_2 = ilon2_2
      ilon2_2 = temp
    endif
  
    start_2d = (/ilon1, ilat1/)
    stride_2d = (/1, 1/)
    edge_2d = (/(ilon2-ilon1)+1, (ilat2-ilat1)+1/)
  
    call read_hdf_sds(id, trim(grid_str%sds_name), start_2d, stride_2d, edge_2d, land_grid)
    
    start_2d_2 = (/ilon1_2, ilat1/)
    stride_2d_2 = (/1, 1/)
    edge_2d_2 = (/(ilon2_2-ilon1_2)+1, (ilat2-ilat1)+1/)
  
    call read_hdf_sds(id, trim(grid_str%sds_name), start_2d_2, stride_2d_2, edge_2d_2, land_grid_2)
    
    do j = 1, ny
      do i = 1, nx
    
        if (space_mask(i,j) == sym%NO_SPACE) then

              ilat = max(1,min(grid_str%num_lat,int(abs(lat(i,j) - grid_str%first_lat)/grid_str%del_lat) + 1))
              ilon = max(1,min(grid_str%num_lon,int(abs(lon(i,j) - grid_str%first_lon)/grid_str%del_lon) + 1))
              if (lon(i,j) >= 0.0) then
                ilat_ad = max(1,min((ilat - start_2d(2)) + 1,size(land_grid,2)))
                ilon_ad = max(1,min((ilon - start_2d(1)) + 1,size(land_grid,1)))
                land(i,j) = land_grid(ilon_ad,ilat_ad)
              else
                ilat_ad = max(1,min((ilat - start_2d_2(2)) + 1,size(land_grid_2,2)))
                ilon_ad = max(1,min((ilon - start_2d_2(1)) + 1,size(land_grid_2,1)))
                land(i,j) = land_grid_2(ilon_ad,ilat_ad)
              endif

        endif
      
      enddo
    enddo
  
    deallocate(land_grid, land_grid_2, stat=astatus)
    if (astatus /= 0) then
      print "(a,'Error deallocating land surface grid.')",EXE_PROMPT
      stop
    endif
    
  endif

END SUBROUTINE read_land_sfc_hdf_i2

!-------------------------------------------------------------------
! Subroutine to read hdf data.
!-------------------------------------------------------------------
  
SUBROUTINE read_hdf_sds_i1(sd_id, sds_name, istart, istride, iedge, buffer)
  INTEGER(kind=int4), intent(in) :: sd_id
  CHARACTER(*), intent(in) :: sds_name
  INTEGER(kind=int4), dimension(2), intent(in) :: istart, istride
  INTEGER(kind=int4), dimension(2), intent(inout) :: iedge
  INTEGER(kind=int1), dimension(:,:), allocatable, intent(out) :: buffer
  INTEGER(kind=int4) :: astatus, istatus, sds_id, sds_rank, sds_type, sds_nattr
  INTEGER(kind=int4), dimension(2) :: sds_dims
    
  integer :: sfselect, sfn2index, sfrdata, sfendacc, sfginfo
  
  sds_id = sfselect(sd_id, sfn2index(sd_id,trim(sds_name)))
    
  istatus = sfginfo(sds_id, sds_name, sds_rank, sds_dims, sds_type, sds_nattr)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    stop
  endif
  
  if (iedge(1) < 0) iedge(1) = sds_dims(1)
  if (iedge(2) < 0) iedge(2) = sds_dims(2)
  
  iedge(1) = min((sds_dims(1) - istart(1)),iedge(1))
  iedge(2) = min((sds_dims(2) - istart(2)),iedge(2))
    
  allocate(buffer(iedge(1),iedge(2)),stat=astatus)
  if (astatus /= 0) then
    print "(a,'Not enough memory to allocate 2d buffer.')",EXE_PROMPT
    stop
  endif
    
  istatus = sfrdata(sds_id, istart, istride, iedge, buffer)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    stop
  endif
  istatus = sfendacc(sds_id)
    
END SUBROUTINE read_hdf_sds_i1

SUBROUTINE read_hdf_sds_i2(sd_id, sds_name, istart, istride, iedge, buffer)
  INTEGER(kind=int4), intent(in) :: sd_id
  CHARACTER(*), intent(in) :: sds_name
  INTEGER(kind=int4), dimension(2), intent(in) :: istart, istride
  INTEGER(kind=int4), dimension(2), intent(inout) :: iedge
  INTEGER(kind=int2), dimension(:,:), allocatable, intent(out) :: buffer
  INTEGER(kind=int4) :: astatus, istatus, sds_id, sds_rank, sds_type, sds_nattr
  INTEGER(kind=int4), dimension(2) :: sds_dims
    
  integer :: sfselect, sfn2index, sfrdata, sfendacc, sfginfo
  
  sds_id = sfselect(sd_id, sfn2index(sd_id,trim(sds_name)))
    
  istatus = sfginfo(sds_id, sds_name, sds_rank, sds_dims, sds_type, sds_nattr)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    stop
  endif
  
  if (iedge(1) < 0) iedge(1) = sds_dims(1)
  if (iedge(2) < 0) iedge(2) = sds_dims(2)
  
  iedge(1) = min((sds_dims(1) - istart(1)),iedge(1))
  iedge(2) = min((sds_dims(2) - istart(2)),iedge(2))
    
  allocate(buffer(iedge(1),iedge(2)),stat=astatus)
  if (astatus /= 0) then
    print "(a,'Not enough memory to allocate 2d buffer.')",EXE_PROMPT
    stop
  endif
    
  istatus = sfrdata(sds_id, istart, istride, iedge, buffer)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    stop
  endif
  istatus = sfendacc(sds_id)
    
END SUBROUTINE read_hdf_sds_i2

END MODULE land_sfc_properties
