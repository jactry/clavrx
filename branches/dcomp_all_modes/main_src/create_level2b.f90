! $Id$


program create_level2b

   use FILE_TOOLS, only: &
         get_lun
      
   use cx_sds_io_mod,only: &
        cx_sds_finfo &
      , cx_sds_read &
      , cx_sds_read_raw &
      , cx_sds_att &
      , cx_att_int &
      , cx_att_r4 
   
   use cx_sds_type_definitions_mod, only: &
      cx_sds_type &
      , cx_att_type &
      , cx_sds_data_type &
      , MAXNCNAM
      
   use cx_grid_tools_mod,only: &
      INDEX_IN_REGULAR_GRID &
      , lonlat_edge_to_array
     
   use cx_hdf_write_mod,only: &
      hdf_file_open &
      , create_sds &
      , write_sds &
      , add_att &
      , close_file &
      , compress_sds &
      , close_sds
   
    use cx_prd_mod, only: &
      prd_dtype &
      , prd_individual_dtype
    
    use constants, only: &
      int4, real4, int1, int2  &
      , MISSING_VALUE_INT1 &
      , MISSING_VALUE_INT2 &
      , MISSING_VALUE_INT4 &
      , MISSING_VALUE_REAL4 
    
    use cx_scale_tools_mod, only: &
      scale_i2_rank2  
    
   use date_tools_mod , only: &
        date_type &
        , time_is_in_window
      
      
   implicit none

   integer :: N_Command_Line_Args
   integer :: ios
   integer :: ifile
   
   integer, PARAMETER :: N_FILES_MAX = 1200
 
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
      logical :: set
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
   
   character ( len = MAXNCNAM), allocatable :: att_name(:)
   character ( len = MAXNCNAM), allocatable :: sds_name(:) 
 
   
   logical :: is_first_file
   
   integer :: i, j, idx
   integer :: num_elem_inp
   integer :: num_line_inp
   
   real, dimension(:,:), allocatable :: Lon_Inp
   real, dimension(:,:), allocatable :: Lat_Inp
   logical, dimension(:,:), allocatable :: Gap_Inp
   integer, allocatable :: ielem_out(:,:)
   integer, allocatable :: iline_out(:,:)
   
   logical, allocatable :: idx_name(:)
   integer :: isds, isds2
   integer (kind = int4) :: num_points
   integer(kind=int4), allocatable:: Input_Update_Index(:)
   integer(kind=int4), allocatable:: Output_Update_Index(:)
   real(kind=real4),  allocatable:: Input_Array_1d(:)
   
   integer(kind=int4 ) :: ilon,ilat,ipoint
   
   logical, allocatable :: pix_to_update(:,:)
   
   integer :: ielem, iline
   real(kind=real4), allocatable:: Scaled_Sds_Data_Output(:,:)
   real(kind=real4), allocatable::  Output_Array_1d(:)
   
   real(kind=real4), allocatable:: data_real(:)
   
   real(kind = 8), pointer :: pp(:,:)
   integer :: test
   type(cx_att_type) :: att
   integer :: ftype
   real , allocatable :: temp_i1_2d(:,:)
   
   integer :: id_file
   integer :: id 
   character (len =1024 ) :: file_level2b
   
   integer :: sds_dims_2d (2), sds_start_2d(2), sds_stride_2d(2)
   
   type(prd_dtype), target :: prd
   type(prd_individual_dtype), pointer:: prd_i
   
   integer :: istatus
   integer :: compress_flag
   
      integer(kind=int4), parameter :: TWO_BYTE_MAX = 32767, & !(2**15)/2 - 1
                                        TWO_BYTE_MIN = -32767   !-(2**15)/2
                                        
   integer(kind=int4), parameter :: one_byte_max = 127, & !(2**8)/2 - 1
                                        one_byte_min = -127   !-(2**8)/2
   
   real :: scale_factor
   real :: Add_Offset
   integer(kind=int2), dimension(:,:),allocatable :: Two_Byte_dummy
   
   
   integer :: start_year
   integer :: end_year
   integer :: start_day_of_year
   integer :: end_day_of_year
   real :: start_time
   real :: end_time
   integer :: wmo_sat
   integer :: l1b  
   
   type(date_type) :: date_file_start,date_file_end
   type(date_type) :: date_cnf_start,date_cnf_end
   
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
   
   compress_flag = 1
   
   config_file_lun = GET_LUN()
   
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
   
   call date_cnf_start.set_date_with_doy (cnf % Year, cnf % Jday &
         , 0,0 )
         
   call date_cnf_end.set_date_with_doy (cnf % Year, cnf % Jday &
         ,23, 59)
   
   do ifile = 1, N_FILES_MAX
       read(unit=config_file_lun,fmt=*,iostat=ios) cnf % File_Input(ifile)
       if (ios /= 0) then
         cnf % n_files = ifile - 1
         exit
      end if
   
   end do
   

    ! -check input
    
   
   ! - configure output arrays and structures   
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
   
      ! - read global variables, assume those are integers
      test = cx_att_int (file,'NUMBER_OF_ELEMENTS',num_elem_inp) 
      test = cx_att_int (file,'NUMBER_OF_SCANS_LEVEL2',Num_Line_Inp)     
      test = cx_att_int (file,'START_YEAR',start_year)
      test = cx_att_int (file,'END_YEAR',end_year)
      test = cx_att_int (file,'START_DAY_OF_YEAR',start_day_of_year)
      if ( test .NE. 0) test = cx_att_int (file,'START_DAY',start_day_of_year)
      test = cx_att_int (file,'END_DAY_OF_YEAR',end_day_of_year)
      if ( test .NE. 0) test = cx_att_int (file,'END_DAY',end_day_of_year)
      test = cx_att_r4 (file,'START_TIME',start_time)
      test = cx_att_r4 (file,'END_TIME',end_time)
      test = cx_att_int (file,'WMO_SATELLITE_CODE',wmo_sat)
      !test = cx_att_cha (file,'L1B',l1b)      
        
     
      
      call date_file_start.set_date_with_doy (start_year, start_day_of_year &
         , floor(start_time), floor(60*(start_time-floor(start_time))) )
         
      call date_file_end.set_date_with_doy (end_year, end_day_of_year &
         , floor(end_time) ,  floor(60*(end_time-floor(end_time))))
      
      
      ! - check time window
      
      if (.not. time_is_in_window ( date_file_start, date_cnf_start, date_cnf_end ) ) then
         print*,'file is not in window'
         print*,'next file'
         cycle
      
      end if 
      
      
      num_points = num_elem_inp * num_line_inp
      allocate(Input_Update_Index(Num_Points))
      allocate(Output_Update_Index(Num_Points))
      allocate(Input_Array_1d(Num_Points))

      ! - this returns information about sds data and aglobal attributes
      if ( cx_sds_finfo ( file, ftype, num_sds_input, sds_name, natt, att_name) < 0 ) then
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
      
      test = cx_sds_read(file, 'longitude',lon_inp)
      test = cx_sds_read(file, 'latitude',lat_inp)
      test = cx_sds_read(file, 'gap_pixel_mask',temp_i1_2d)
      
      gap_inp = temp_i1_2d .GT. 0
      

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
      
      
      deallocate ( Ielem_out)
      deallocate ( Iline_out)
      
      num_sds_output = size ( sds_out)
      
      ! - loop over input variables
      do Isds = 1, Num_Sds_Input   
         sds_out(isds) % set = .false.
         test = cx_sds_read ( file, sds_name (isds) , temp_i1_2d)
         if ( sds_name (isds) .EQ. 'scan_line_number' ) cycle
         if ( sds_name (isds) .EQ. 'scan_line_time' ) cycle
         if ( sds_name (isds) .EQ. 'bad_scan_line_flag' ) cycle
         if ( sds_name (isds) .EQ. 'asc_des_flag' ) cycle
         sds_out(isds) % set = .true.
         if (is_First_file) then
            
            sds_out(isds) % name = sds_name (isds)
            allocate (sds_out(isds) % data (Nlon_Out, Nlat_Out) )
            sds_out(isds) % data = -999.
         end if 
        
         ! - loop over already defined output variables
         ! - this is done to get sure to match identical variables even 
         ! - if the level2 files are not identical
         do isds2 = 1, num_sds_output
          
            if ( .NOT. sds_name (isds) .EQ. sds_out(isds2) % name ) cycle
            
            Input_Array_1d = reshape(temp_i1_2d, (/Num_Elem_Inp * Num_Line_Inp/)) 
                      
            Scaled_Sds_Data_Output = sds_out(isds) % data 
             
            Output_Array_1d = reshape(Scaled_Sds_Data_Output, (/Nlon_Out*Nlat_Out/))        
             
            Output_Array_1d(Output_Update_Index(1:ipoint)) = Input_Array_1d(Input_Update_Index(1:ipoint))
            
            sds_out(isds) % data  = reshape(Output_Array_1d, (/Nlon_Out,Nlat_Out/))
             

         end do
          
      end do
      
     
      deallocate(Input_Update_Index)
      deallocate(Output_Update_Index)
      deallocate(Input_Array_1d)
      
      deallocate(Lon_Inp)
      deallocate(Lat_Inp)
      deallocate(Gap_inp)
      
     
      
      
      is_first_file = .false.
   end do file_loop
   
   
   
   
! ++++++++++++++ write in level2b  +++++++++++++++++
! we have now all data stored in sds_out
! lets write this in level2b file
!  
!    open and create file
!    do i = 1, num_sds_output
!   File_Level2b= trim(File_2b_Root)//".level2b.hdf"   
   
   call prd % read_products()
         
   File_Level2b= 'test_level2b.hdf'
   id_file = hdf_file_open(trim(File_Level2b), create=.true.)
   
   sds_dims_2d = [nlon_out,nlat_out]
   sds_start_2d =[0,0]
   sds_stride_2d = [1,1]
  
   
   do i = 1,  num_sds_output 
      ! - check if this is in csv ( should be)
      idx = -999
      do j = 1 , prd % num_products
         if ( trim(prd % product(j) % name) .EQ. trim(sds_out(i) % name)) idx = j 
      end do  
      if (idx .eq. -999) cycle
      
      ! - when found set pointer 
      prd_i => prd % product(idx)
      
      ! first define and create sds variables
      select case ( prd_i % dim)
      case(1)
         prd_i % Sds_Id= create_sds (id_file, prd_i % name ,  nlon_out, prd_i % dtype)
      case(2)
         prd_i % Sds_Id= create_sds (id_file, prd_i % name , sds_dims_2d , prd_i % dtype)
         istatus = compress_sds ( prd_i % Sds_Id,Compress_Flag, sds_dims_2d) 
      end select
      
      call add_att( prd_i % Sds_Id, 'SCALED', prd_i % scaling)
      call add_att( prd_i % sds_id, 'unit', trim(prd_i % unit)) 
      call add_att( prd_i % sds_id, 'standard_name', trim(prd_i % standard_name))
      call add_att( prd_i % sds_id, 'long_name', trim(prd_i % long_name))
      
      if ( prd_i % scaling .gt. 0 ) then
            
         select case ( prd_i % dtype ) 
         case(1)
            scale_factor = (prd_i%act_max - prd_i%act_min )/(one_byte_max - one_byte_min)
            add_offset = prd_i%act_min - scale_factor * one_byte_min
            call add_att ( prd_i % sds_id, 'scaled_missing',Missing_Value_Int1 )
         case(2)
            scale_factor = (prd_i%act_max - prd_i%act_min )/(two_byte_max - two_byte_min)
            add_offset = prd_i%act_min - scale_factor * two_byte_min 
            call add_att ( prd_i % sds_id, 'scaled_missing',Missing_Value_Int2 )  
         end select   
              
            
         call add_att ( prd_i % sds_id, 'add_offset', add_offset)
         call add_att ( prd_i % sds_id, 'scale_factor',scale_factor)
         call add_att ( prd_i % sds_id, 'actual_min',prd_i%act_min )
         call add_att ( prd_i % sds_id, 'actual_max',prd_i%act_max )
      end if
      
      
      ! - write
      select case (prd_i % dim)
       
      case(2)
         select case (prd_i % dtype)
         
         case(1)
            if (prd_i % scaling == 1 ) then
                     
                 !    call SCALE_VECTOR_I1_RANK2(sds_out(i) % data,prd_i % scaling &
                 !       ,prd_i % act_min,prd_i % act_max,Missing_Value_Real4 &
                 !       ,sds_out(i) % data)
                 !    Istatus = write_sds ( prd_i % sds_id,[0,0],[1,1],sds_dims_2d, &
                 !       sds_out(i) % data) + Istatus 
            else
                     
                 !    Istatus = write_sds ( prd_i % sds_id, Sds_Start_2d, Sds_Stride_2d, Sds_dims_2d, &
                 !              data_dim2_dtype1) + Istatus
            end if
                  
         case(2)
            if (prd_i % scaling == 1) then
               call scale_i2_rank2(sds_out(i) % data  &
                     ,prd_i % act_min,prd_i % act_max,Missing_Value_Real4,two_byte_dummy)

               Istatus = write_sds ( prd_i % sds_id,Sds_Start_2d,Sds_Stride_2d,Sds_dims_2d, &
                        Two_Byte_Dummy ) + Istatus
               if ( allocated (two_byte_dummy) ) deallocate(two_byte_dummy)  
            else 
                 !    Istatus = write_sds ( prd_i % sds_id, Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                 !              data_dim2_dtype2) + Istatus
            end if 
                  
         case(4)
                  !- this is only diagnostic variables, we don't scale..
            Istatus = write_sds ( prd_i % sds_id, Sds_Start_2d, Sds_Stride_2d, Sds_dims_2d,  &
                      sds_out(i) % data ) + Istatus   
         end select   
      end select
       
      call close_sds(prd_i % sds_id)                 
   end do
   
   call close_file (id_file)    
   


!     end do 


   
   
  
end program create_level2b
