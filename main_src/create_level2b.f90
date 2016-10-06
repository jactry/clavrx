! $Id$


program create_level2b
      use FILE_TOOLS, only: &
         getlun
      
      use cx_hdf_read_mod, only: &
         hdf_att &
      , MAXNCNAM &
      , hdf_get_file_att &
      , hdf_get_finfo &
      , hdf_get_file_sds &
      , hdf_sds &
      ,  hdf_data &
      , transform_sds_to_real
      
      
      use cx_grid_tools_mod
      
      
   implicit none
   
   
   
   
   integer :: N_Command_Line_Args
   integer :: ios
   integer :: ifile
   
   integer, PARAMETER :: N_FILES_MAX = 1200
  ! integer, parameter:: int4 = selected_int_kind(8)
   type l2b_config 
      character ( len = 3) :: node = 'asc'
      character ( len = 4) :: spatial = 'near'
      character ( len = 14) :: overlap ='nadir_overlap'
      character ( len = 10) :: geo = '1d_lat_lon'
      character(len=1020):: Dir_In
      character(len=1020):: Dir_Out
      integer :: Year
      integer :: jday
      integer :: time
      integer :: time_window
      character(len=12) :: sensor_name
      real :: lon_west
      real :: lon_east
      real ::lat_south
      real :: lat_north
      real :: dlon_output
      real :: dlat_output
      character(len=1020) :: file_input (N_FILES_MAX)
      integer :: n_files
      
   end type
   
   type ( l2b_config ) :: cnf
   
   type sds_out_type
      real(kind = real4), allocatable :: data(:,:)
      character(len = MAXNCNAM) :: name
   end type
   
   type(sds_out_type), allocatable :: sds_out(:)
   
   
   integer:: config_file_lun
   logical :: has_dateline
   
   integer :: nlon_out
   integer :: nlat_out
   real, allocatable :: lat_out(:,:)
   real, allocatable :: lon_out(:,:)
   
   
   integer :: num_sds_input
   integer :: num_sds_output
   character(len=240) :: file
   integer :: natt,nsds
   type(hdf_att), allocatable:: attrs(:)
   character ( len = MAXNCNAM), allocatable :: att_name(:)
   character ( len = MAXNCNAM), allocatable :: sds_name(:) 
   type(hdf_sds), dimension(:), target, allocatable :: sds_new  
   
   logical :: is_first_file
   
   integer :: i
   integer :: num_elem_inp
   integer :: num_line_inp
   
  
   type(hdf_sds), pointer :: ps                          
   type(hdf_data), pointer :: psd   
   real, dimension(:,:), allocatable :: Lon_Inp
   real, dimension(:,:), allocatable :: Lat_Inp
   logical, dimension(:,:), allocatable :: Gap_Inp
   integer, allocatable :: ielem_out(:,:)
   integer, allocatable :: iline_out(:,:)
   
   logical, allocatable :: idx_name(:)
   integer :: isds, isds2
   integer (kind = int4) :: num_points
   integer(kind=int4), dimension(:), allocatable:: Input_Update_Index
   integer(kind=int4), dimension(:), allocatable:: Output_Update_Index
   real(kind=real4), dimension(:), allocatable:: Input_Array_1d
   
   integer(kind=int4 ) :: ilon,ilat,ipoint
   
   logical, allocatable :: pix_to_update(:,:)
   
   integer :: ielem, iline
   real(kind=real4), dimension(:,:), allocatable:: Scaled_Sds_Data_Output
   real(kind=real4), dimension(:), allocatable::  Output_Array_1d
   
   real(kind=real4), dimension(:), allocatable:: data_real
   
   
   ! ++++++++
   
   !  check users desires
   N_Command_Line_Args = iargc()
   if (N_Command_Line_Args .GE. 1) call getarg(1, cnf % Node)
   if (N_Command_Line_Args .GE. 2) call getarg(2, cnf % Spatial)
   if (N_Command_Line_Args .GE. 3) call getarg(3, cnf % Overlap)
   if (N_Command_Line_Args .EQ. 4) call getarg(4, cnf % Geo)
   if (N_Command_Line_Args .GT. 4 ) then
      print *,  "COMP_ASC_DES_LEVEL2b ERROR:: Unexpected Number of Command Line Arguments"
   end if 
   
   
   config_file_lun = GETLUN()
   
   open(unit=config_file_lun &
      ,file="comp_asc_des_level2b_input" &
      ,status="old" &
      ,action="read" &
      ,position="rewind" &
      ,iostat=ios)
      
   if (ios /= 0) then
      print *,  'error opening comp_asc_dec control file, ios = ', ios
      stop 1
   else
      read(unit=config_file_lun,fmt="(A)") cnf % dir_in
      read(unit=config_file_lun,fmt="(A)") cnf % dir_out
      read(unit=config_file_lun,fmt=*) cnf % Year
      read(unit=config_file_lun,fmt=*) cnf % Jday
      read(unit=config_file_lun,fmt=*) cnf % Time, cnf % Time_Window
      read(unit=config_file_lun,fmt=*) cnf % Sensor_Name
      read(unit=config_file_lun,fmt=*) cnf % Lon_West,cnf %  Lon_East,cnf %  Dlon_output
      read(unit=config_file_lun,fmt=*) cnf % Lat_South, cnf % Lat_North, cnf % Dlat_output
   endif
   
   
   do ifile = 1, N_FILES_MAX
       read(unit=config_file_lun,fmt=*,iostat=ios) cnf % File_Input(ifile)
       if (ios /= 0) then
         cnf % n_files = ifile - 1
         exit
      end if
   
   end do
   

    ! -check input
    
   
   ! - confgure output arrays and structures   
   has_dateline = .false.
   
   if ( cnf % lon_west .gt. cnf % lon_east ) has_dateline = .true.
   if (has_dateline) then
      nlon_out = nint( ((360.0 + cnf % Lon_East) - cnf % Lon_West) / cnf % dlon_Output) + 1
   else 
      nlon_out = nint (( cnf % Lon_East - cnf % Lon_West ) / cnf % dlon_Output) + 1
   end if
   
   Nlat_Out = nint(( cnf % Lat_North - cnf % Lat_South) / cnf % Dlat_Output) + 1
   
   allocate ( lat_out (nlon_out,nlat_out))
   allocate ( lon_out (nlon_out,nlat_out))
   
   ! - computes the level-2b array lon/lat arrays
   call lonlat_edge_to_array (  &
         cnf % Lat_South &
        , cnf % dlat_output      &
        , Nlat_out      &
        , cnf % Lon_West  &
        , cnf % dlon_output      &
        , Nlon_out      &
        , lat_out &
        , lon_out )
   
   !  read input
   !---------------------------------------------------------------------------------
   ! loop through orbit files
   !---------------------------------------------------------------------------------
   
   allocate ( pix_to_update (nlon_out,nlat_out))
   allocate(Scaled_Sds_Data_Output(Nlon_Out,Nlat_Out))
   allocate(Output_Array_1d(Nlon_Out*Nlat_Out))
   
   pix_to_update = .true.
   is_first_file = .true.
   
   
   file_loop: do ifile = 1, cnf % n_files
      print *, "processing file ", ifile, " of ", cnf % n_files, " ", trim(cnf % File_Input(ifile))
      file = trim(trim(cnf%dir_in)//trim(cnf%File_Input(ifile)))
   
      if ( hdf_get_file_att (file,natt,attrs) < 0 ) then
         print*,'error reading attributes'
         stop
      end if
      
      
      do i = 1, natt
            
         if ( attrs(i) % name .eq. 'NUMBER_OF_ELEMENTS' ) Num_Elem_Inp = attrs(i) % data % i2values (1)
         if ( attrs(i) % name .eq. 'NUMBER_OF_SCANS_LEVEL2' ) Num_Line_Inp = attrs(i) % data % i2values (1)
         !if ( attrs(i) % name .eq. 'START_YEAR' ) Start_Year_Input = attrs(i) % data % i2values (1)
         !if ( attrs(i) % name .eq. 'END_YEAR' ) End_Year_Input = attrs(i) % data % i2values (1)
         !if ( attrs(i) % name .eq. 'START_DAY_OF_YEAR' ) Start_Day_Input = attrs(i) % data % i2values (1)
         !if ( attrs(i) % name .eq. 'END_DAY_OF_YEAR' ) End_Day_Input = attrs(i) % data % i2values (1)
         !if ( attrs(i) % name .eq. 'START_TIME' ) Start_Time_Input = attrs(i) % data % r4values (1)
         !if ( attrs(i) % name .eq. 'END_TIME' ) End_Time_Input = attrs(i) % data % r4values (1)
         !if ( attrs(i) % name .eq. 'WMO_SATELLITE_CODE' ) Sc_Id_Input = attrs(i) % data % i2values (1)
         !if ( attrs(i) % name .eq. 'L1B' ) L1b_input = attrs(i) % data % c1values (1)
         
      end do
      
      
      num_points = num_elem_inp * num_line_inp
      allocate(Input_Update_Index(Num_Points))
      allocate(Output_Update_Index(Num_Points))
      allocate(Input_Array_1d(Num_Points))
      
      ! - this retuens information about sds data and aglobal attributes
      if ( hdf_get_finfo(file, num_sds_input, sds_name, natt, att_name) < 0 ) then
         print*,'error reading attributes'
         stop
      end if
      
      if (is_first_file) then
         allocate (sds_out(Num_Sds_Input))        
         
      end if
      
      ! - read some inportant variables to check if we have right data before read all..
      
      allocate(Lon_Inp(Num_Elem_Inp,Num_Line_Inp))
      allocate(Lat_Inp(Num_Elem_Inp,Num_Line_Inp))
      allocate(Gap_Inp(Num_Elem_Inp,Num_Line_Inp))
      
      
      if (hdf_get_file_sds(file, nsds, sds_new, nsdsn = 3, sds_name = ['longitude','latitude','gap_pixel_mask']) < 0) then
         print*,'hdf file not readable ', trim(file)
         stop  
      end if 
      
      ps => sds_new(1) ; psd=> ps%data
      lon_inp =  reshape(psd%i2values,[Num_Elem_Inp,Num_Line_Inp])   * (ps % attr(7) % data % r4values(1) )   
      
      ps => sds_new(2) ; psd=> ps%data
      lat_inp =  reshape(psd%i2values,[Num_Elem_Inp,Num_Line_Inp])   * (ps % attr(7) % data % r4values(1) ) 
      
      ps => sds_new(3) ; psd=> ps%data
      Gap_Inp = (reshape(psd%i1values  ,[Num_Elem_Inp,Num_Line_Inp]) .gt. 0)
      
      
      allocate ( Ielem_out(nlon_out,nlat_out))
      allocate ( Iline_out(nlon_out,nlat_out))
      ! - compute indicies
      call INDEX_IN_REGULAR_GRID ( &
          cnf % Lat_South &
        , cnf % dlat_output      &
        , Nlat_out      &
        , cnf % Lon_West  &
        , cnf % dlon_output      &
        , Nlon_out      &
        , cnf % spatial == 'near' &
        , lon_inp &
        , lat_inp &
        , gap_inp &
        , Ielem_out &
        , Iline_out )
        
      pix_to_update = Ielem_out .gt. 0  
      
      ! - set flags that let us decide if we take this pixel or what we have already in
      
        !----- compute 1d input to output index mapping 
      ipoint = 0
      do Ilon = 1, Nlon_Out
         do Ilat = 1, Nlat_Out
            if (pix_to_update(Ilon,Ilat)) then
               ipoint = ipoint + 1
               ipoint = max(0,min(Ipoint,Num_Points))
               Ielem = Ielem_Out(Ilon,Ilat)
               Iline = Iline_Out(Ilon,Ilat)
               Input_Update_Index(Ipoint) = Ielem + (Iline-1)*num_elem_inp
               
               Output_Update_Index(Ipoint) = Ilon + (Ilat-1)*Nlon_Out
              ! Scan_Element_Number_Output(Ilon,Ilat) = Ielem
            end if
         end do
      end do
      
      
      ! - read science data
      if (hdf_get_file_sds(file, nsds, sds_new, nsdsn = num_sds_input , sds_name = sds_name, dclb = .true.) < 0) then
            print*,'hdf file not readable ', trim(file)
            stop  
      end if 
      
      num_sds_output = size ( sds_out)
      
      ! - loop over input variables
      do Isds = 1, Num_Sds_Input   
         ps => sds_new(Isds) ; psd=> ps%data
        
         if (is_First_file) then
            sds_out(isds) % name = ps % name
            allocate (sds_out(isds) % data (Nlon_Out, Nlat_Out) )
            sds_out(isds) % data = -999.
         end if 
         
         ! - loop over already defined output variables
         ! - this is done to get sure to match identical variables even 
         ! - the level2 files are not identical
         do isds2 = 1, num_sds_output
            if ( .NOT. ps % name .EQ. sds_out(isds2) % name ) cycle
            print*,isds2,trim(sds_out(isds2) % name)
            print*,psd % calibr
            
           ! do i = 1, ps%nattr
           !    print*,trim(ps%attr(i)%name)
           !    print*,ps%attr(i)%data%r4values
           ! end do
            
            
            call transform_sds_to_real ( psd, Input_Array_1d)
            
           ! Input_Array_1d = reshape(data_real, (/Num_Elem_Inp * Num_Line_Inp/)) 
             
            Scaled_Sds_Data_Output = sds_out(isds) % data 
            
            Output_Array_1d = reshape(Scaled_Sds_Data_Output, (/Nlon_Out*Nlat_Out/))  
            
            Output_Array_1d(Output_Update_Index(1:ipoint)) = Input_Array_1d(Input_Update_Index(1:ipoint))
           
            sds_out(isds) % data  = reshape(Output_Array_1d, (/Nlon_Out,Nlat_Out/))
            
            print*,sds_out(isds) % data (2170,644),maxval(sds_out(isds) % data ),minval(sds_out(isds) % data )
            
               
              
              
               
         end do
          
      end do
      
      !- normally number is the same, but we take the number of first file in case it is different
      
      print*,'success'
      
   
      
      deallocate(Lon_Inp)
      deallocate(Lat_Inp)
      deallocate(Gap_inp)
      
      deallocate ( Ielem_out)
      deallocate ( Iline_out)
      
      
      is_first_file = .false.
   end do file_loop
   
   
   
   
! ++++++++++++++ write in level2b  +++++++++++++++++


end program create_level2b
