module CLAVRX_STATIC_NAV_MODULE

use NETCDF

integer:: num_elem_fd, elem_start, num_elem
integer:: num_line_fd, lin_start, num_line

character(len=120), parameter:: EXE_PROMPT_NAV = "CLAVR-x Static Navigation Module >> "

real, dimension(:,:), save, private:: latitude_fd
real, dimension(:,:), save, private:: longitude_fd
integer, save, private:: num_elem_fd, num_line_fd

implicit none

contains
!-------------------------------------------------------------------------------------
! begin executable code
!-------------------------------------------------------------------------------------
subroutine READ_CLAVRX_STATIC_NAV(nav_file_name)

   character (len=120, intent(in):: nav_file_name 
   integer:: ncid, status

   !--- these are input
   nav_file_name = 'clavrx_ahi_fulldisk_nav.nc'
   lat_north = 30.0
   lat_south = -30.0
   lon_west = 100.0
   lon_east = 130.0
   missing = -999.0
   


   !-----------------------------------------------------------------------------------
   !  Read data
   !-----------------------------------------------------------------------------------
   status = nf90_open(nav_file_name, mode = nf90_nowrite, ncid = ncid)

   if (status /= nf90_noerr) then
      print *, EXE_PROMPT_NAV , 'ERROR: Static Navigation File Read Failed'
      print*, EXE_PROMPT_NAV ,' filename is: ', Nav_File_Name
      return
   endif

   status = nf90_get_att(ncid, nf90_global, "number_of_elements", num_elem_fd)
   status = nf90_get_att(ncid, nf90_global, "number_of_lines", num_line_fd)
   status = nf90_get_att(ncid, nf90_global, "sub_satellite_latitude", sub_sat_lat)
   status = nf90_get_att(ncid, nf90_global, "sub_satellite_longitude", sub_sat_lon)


   allocate(Latitude_Fd(num_elem_fd, num_line_fd))
   allocate(Longitude_Fd(num_elem_fd, num_line_fd))

   call read_netcdf_2d_real(ncid, (/1,1/), (/num_elem, num_line/),"latitude",latitude_fd)
   call read_netcdf_2d_real(ncid, (/1,1/), (/num_elem, num_line/),"longitude",longitude_fd)

   status = nf90_close(ncid)

   !-----------------------------------------------------------------------------------
   ! Check to See if Sub Sat Position is as Expected, if not fail
   !-----------------------------------------------------------------------------------

end subroutine READ_CLAVRX_STATIC_NAV

!------------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------------
subroutine DETERMINE_BOUNDS_STATIC_NAV(lat_south, lon_west, lat_north, lon_east, $
                                       elem_start, num_elem, line_start, num_line)

   real, intent(in):: lat_south, lon_west, lat_north, lon_east
   integer, intent(out):: elem_start, num_elem, line_start, num_line
   real, dimension(:), allocatable:: lat_nadir, lat_diff, lon_nadir, lon_diff
   integer:: i_start, i_end, j_start, j_end

   allocate(lat_nadir(num_line_fd))
   allocate(lat_diff(num_line_fd))
   allocate(lon_nadir(num_line_fd))
   allocate(lon_diff(num_line_fd))

   print *, 'computing subset bounds'
   j_north = missing
   j_south = missing
   i_west = missing
   i_east = missing

   lat_nadir = lat(num_elem_fd/2,:)
   lon_nadir = lon(:,num_line_fd/2)

   lat_diff = abs(lat_nadir - lat_north)
   idx = where(lat_diff eq min(lat_diff),cc)

   lat_diff = abs(lat_nadir - lat_south)
   j_south = minloc(lat_diff)

   lon_diff = abs(lon_nadir - lon_west)
   i_west = minloc(lon_diff)

   lon_diff = abs(lon_nadir - lon_east)
   i_east = minloc(lon_diff)

   ;---- reorder as necessary
   i_start = min(i_west,i_east)
   i_end = max(i_west,i_east)
   j_start = min(j_south,j_north)
   j_end = max(j_south,j_north)

   if (i_start lt 1 .or. i_end lt 1) then
     print *, 'longitude subsetting failed'
     i_start = 1
     i_end = num_elem_fd
   endif

   if (j_start lt 1 .or. j_end lt 1) then 
     print *, 'latitude subsetting failed'
     j_start = 1
     j_end = num_line_fd
   endif


   !--- configure output
   elem_start = i_start
   num_elem = i_end - i_start + 1

   line_start = j_start
   num_line = j_end - j_start + 1

   !--- clean up memory
   deallocate(lat_nadir)
   deallocate(lat_diff)
   deallocate(lon_nadir)
   deallocate(lon_diff)

end subroutine DETERMINE_BOUNDS_STATIC_NAV

!------------------------------------------------------------------------------
!  return lat and lon for this segment
!------------------------------------------------------------------------------
subroutine GET_DATA_STATIC_NAV(segment_number, segment_size,  &
                               elem_start, num_elem, line_start, num_line, &
                               lat, lon) 
   integer, intent(in):: segment_number, segment_size, elem_start, num_elem, line_start, num_line
   real, dimension(:,:), intent(out):: lat, lon
   integer:: i_start, i_end, j_start, j_end


   i_start = elem_start
   i_end = num_elem

   j_start = (segment_number - 1) * segment_size + 1 
   j_end = min(num_line_fd, j_start + segment_size)

   lat = latitude_fd(i_start:i_end,j_start:j_end)
   lon = longitude_fd(i_start:i_end,j_start:j_end)

end subroutine GET_DATA_STATIC_NAV

!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
subroutine RESET_STATIC_NAV()
   deallocate(Latitude_Fd)
   deallocate(Longitude_Fd)
end subroutine RESET_STATIC_NAV


end module CLAVRX_STATIC_NAV_MODULE
