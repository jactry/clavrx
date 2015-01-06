!----------------------------------------------------------------------
! MODULE name: NETCDF_READ_MODULE
! 
! Routines for reading in netCDF files, such as AHI data
!
! Authors: Andrew Heidinger, NOAA/NESDIS
!          Andi Walther, CIMSS
!          Denis Botambekov, CIMSS
!          William Straka, CIMSS
!
! DEPENDENCIES: Constants, NETCDF library
!
!----------------------------------------------------------------------
module NETCDF_READ_MODULE

 use CONSTANTS
 use NETCDF
 
 public:: read_netcdf_1d_real
 public:: read_netcdf_1d_int
 public:: read_netcdf_2d_real
 public:: read_netcdf_2d_int
 public:: read_netcdf_2d_char
 public:: read_netcdf_3d
 public:: read_netcdf_attribute
 
 subroutine read_netcdf_attribute_real (nc_file_id, attribute_id, var_name, attr)
   integer:: nc_file_id
   integer:: attribute_id
   character(30):: var_name
   real, intent(out) :: attr
   integer:: status

   status = nf90_get_att(nc_file_id, attribute_id, trim(var_name), attr)
 
 
 end subroutine
 
 
 

   ! ----------------------------------------------------------
   ! Read in 1D arrays (used code from DCOMP reader
   ! ----------------------------------------------------------
   subroutine read_netcdf_1d_real (nc_file_id, var_dim, var_name, var_output)
        implicit none
      integer, intent(in) :: nc_file_id
      integer, intent(in) :: var_dim
      character(30), intent(in) :: var_name
      real, intent(out), dimension(:) :: var_output

      integer :: nc_var_id
      integer :: status

      Sds_Start_1D = 1
      Sds_Stride_1D = 1
      Sds_Edge_1D = var_dim

      status = nf90_inq_varid(nc_file_id,trim(var_name), nc_var_id)
      if (status /= nf90_noerr) then
            print *, "Error: Unable to get variable id for ", trim(var_name)
            return
      endif

      !get Variable
      status = nf90_get_var(nc_file_id, nc_var_id, var_output, start=Sds_Start_1D, count=Sds_Edge_1D)
      if (status /= nf90_noerr) THEN
            print *,'Error: ',  trim(nf90_strerror(status)),'   ', trim(var_name)
            return
      ENDIF

   end subroutine read_netcdf_1d_real                                                                                                                           

   ! ----------------------------------------------------------
   ! Read in 1D arrays (used code from DCOMP reader
   ! ----------------------------------------------------------
   subroutine read_netcdf_1d_int (nc_file_id, var_dim, var_name, var_output)
        implicit none
      integer, intent(in) :: nc_file_id
      integer, intent(in) :: var_dim
      character(len=*), intent(in) :: var_name
      integer, intent(out), dimension(:) :: var_output

      integer :: nc_var_id
      integer :: status

      Sds_Start_1D = 1
      Sds_Stride_1D = 1
      Sds_Edge_1D = var_dim

      status = nf90_inq_varid(nc_file_id,trim(var_name), nc_var_id)
      if (status /= nf90_noerr) then 
            print *, "Error: Unable to get variable id for ", trim(var_name)
            return
      endif

      !get Variable
      status = nf90_get_var(nc_file_id, nc_var_id, var_output, start=Sds_Start_1D, count=Sds_Edge_1D)
      if (status /= nf90_noerr) THEN
            print *,'Error: ',  trim(nf90_strerror(status)),'   ', trim(var_name)
            return
      ENDIF

   end subroutine read_netcdf_1d_int

   ! ----------------------------------------------------------
   ! Read in 2D arrays (used code from DCOMP reader)
   ! ----------------------------------------------------------
   subroutine read_netcdf_2d_real (nc_file_id, start_var, var_dim, var_name, var_output)
        implicit none
      integer, intent(in) :: nc_file_id
      integer, intent(in) :: start_var(:)
      integer, dimension(:), intent(in) :: var_dim

      character(len=*), intent(in) :: var_name
      real, intent(out), dimension(:,:) :: var_output

      integer :: nc_var_id
      integer :: status

      status = nf90_inq_varid(nc_file_id, trim(var_name), nc_var_id)
      if (status /= nf90_noerr) then
            print *, "Error: Unable to get variable id for ", trim(var_name)
            return
      endif

      !get Variable
      status = nf90_get_var(nc_file_id, nc_var_id, var_output, start=start_var, count=var_dim)
      if ((status /= nf90_noerr)) THEN
            print *,'Error: ',  trim(nf90_strerror(status)),'   ', trim(var_name)
            return
      ENDIF

   end subroutine read_netcdf_2d_real

   ! ----------------------------------------------------------
   ! Read in 2D arrays Integers
   ! ----------------------------------------------------------
   subroutine read_netcdf_2d_int (nc_file_id, start_var, var_dim, var_name, var_output)
        implicit none
      integer, intent(in) :: nc_file_id
      integer, intent(in) :: start_var(:)
      integer, dimension(:), intent(in) :: var_dim

      character(len=*), intent(in) :: var_name
      integer , intent(out), dimension(:,:) :: var_output

      integer :: nc_var_id
      integer :: status

      status = nf90_inq_varid(nc_file_id, trim(var_name), nc_var_id)
      if (status /= nf90_noerr) then
            print *, "Error: Unable to get variable id for ", trim(var_name)
            return
      endif

      !get Variable
      status = nf90_get_var(nc_file_id, nc_var_id, var_output, start=start_var, count=var_dim)
      if ((status /= nf90_noerr)) THEN
            print *,'Error: ',  trim(nf90_strerror(status)),'   ', trim(var_name)
            return
      ENDIF

   end subroutine read_netcdf_2d_int

   ! ----------------------------------------------------------
   ! Read in 2D arrays Characters
   ! ----------------------------------------------------------

   subroutine read_netcdf_2d_char (nc_file_id, start_var, var_dim, var_name, var_output)
        implicit none
      integer, intent(in) :: nc_file_id
      integer, intent(in) :: start_var(:)
      integer, dimension(:), intent(in) :: var_dim

      character(len=*), intent(in) :: var_name
      character(len=30) , intent(out), dimension(:,:) :: var_output
      character(len=30), allocatable, dimension(:,:) :: var

      integer :: nc_var_id
      integer :: status, tmp1, tmp2, i, j
      integer, dimension(2) ::dimIDs

      status = nf90_inq_varid(nc_file_id, trim(var_name), nc_var_id)
      if (status /= nf90_noerr) then
            print *, "Error: Unable to get variable id for ", trim(var_name)
            return
      endif

      !find dimentions
      status = nf90_inquire_variable(nc_file_id, nc_var_id, dimids = dimIDs)
      status = nf90_inquire_dimension(nc_file_id, dimIDs(1), len = tmp1)
      status = nf90_inquire_dimension(nc_file_id, dimIDs(2), len = tmp2)
      allocate (var(tmp1,tmp2))

      !get Variable
      status = nf90_get_var(nc_file_id, nc_var_id, var, start=(/1,1/), count=(/tmp1,tmp2/) )
      if ((status /= nf90_noerr)) THEN
            print *,'Error: ',  trim(nf90_strerror(status)),'   ', trim(var_name)
            return
      ENDIF

      !extract and save classifier names to the final array
      do i = 1, tmp2
        if ((var(i,1) .ge. 'a' .and. var(i,1) .le. 'z') &
        .or.(var(i,1) .ge. 'A' .and. var(i,1) .le. 'Z')) then 
           var_output(i,:) = trim(var(i,1))
        endif
      enddo

      if (allocated(var)) deallocate (var)


   end subroutine read_netcdf_2d_char

   ! ----------------------------------------------------------
   ! Read in 3D arrays (used code from DCOMP reader
   ! ----------------------------------------------------------
   subroutine read_netcdf_3d (nc_file_id, start_var, var_dim, var_name, var_output)
         implicit none
      integer, intent(in) :: nc_file_id
      integer, intent(in) :: start_var(:)
      integer, dimension(:), intent(in) :: var_dim

      character(len=30), intent(in) :: var_name
      real, intent(out), dimension(:,:,:) :: var_output

      integer :: nc_var_id
      integer :: status = 0

      status = nf90_inq_varid(nc_file_id, trim(var_name), nc_var_id)
      if (status /= nf90_noerr) then
            print *, "Error: Unable to get variable id for ", trim(var_name)
            return
      endif

      !get Variable
      status = nf90_get_var(nc_file_id, nc_var_id, var_output, start=start_var, count=var_dim)
      if ((status /= nf90_noerr)) THEN
            print *,'Error: ',  trim(nf90_strerror(status)),'   ', trim(var_name)
            return
      ENDIF


   end subroutine read_netcdf_3d
 
 


end module NETCDF_READ_MODULE
