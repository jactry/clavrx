! $Id:$

module area_file_mod
   
   use file_tools
   
   type, public:: area_header_type
      integer :: Area_Status           ! Area status
      integer  :: Swap_Bytes !whether or not bytes need to be swapped
      integer :: Version_Num           ! Area version number
      integer :: Sat_Id_Num            ! Satellite ID (SSS)
      integer :: Img_Date              ! Year and Julian day (YYDDD)
      integer :: Img_Time              ! Time of image (HHMMSS)
      integer :: North_Bound           ! Upper left line in Sat coords (Y-Coord)
      integer :: West_Vis_Pixel        ! Upper left element in Sat coords (X-Coord)
      integer :: Z_Coor                ! Upper z-coord (Z-Coord)
      integer :: Num_Line              ! Number of lines in image (Y-SIZE)
      integer :: Num_Elem              ! Number of elememts in image (X-SIZE)
      integer :: Bytes_Per_Pixel       ! Number of bytes per data element
      integer :: Line_Res              ! Line resolution (Y-RES) 
      integer :: Elem_Res              ! Element resolution (X-RES)
      integer :: Num_Chan              ! Number of bands (Z-RES)
      integer :: Num_Byte_Ln_Prefix    ! Number of bytes in line prefix (multiple of 4)
      integer :: Proj_Num              ! Project number
      integer :: Creation_Date         ! Creation date (YYDDD)
      integer :: Creation_Time         ! Creation time (HHMMSS)
      integer :: Sndr_Filter_Map       ! Filter map for soundings
      integer :: img_id_num            ! Image ID number
      integer, dimension(4) :: id     ! Reserved for radar appications
      character(len=32):: comment     ! 32 char comments
      integer :: pri_key_calib         ! Calibration colicil (area number)
      integer :: pri_key_nav           ! Primary navigation codicil  (data)
      integer :: sec_key_nav           ! Secondary navigation codicil (nav)
      integer :: val_code              ! Validity code
      integer, dimension(8) :: pdl    ! PDL in packed-byte format
      integer :: band8                 ! Where band-8 came from
      integer :: act_Img_Date          ! Actual image start day (YYDDD)
      integer :: act_Img_Time          ! Actual image start time (HHMMSS)
      integer :: act_Start_Scan        ! Actual start scan
      integer :: len_prefix_doc        ! Length of prefix documentation (DOC)
      integer :: len_prefix_calib      ! Length of prefix calibration (CAL)
      integer :: len_prefix_lev        ! Length of prefix level (LEV)
      character(len=4)::src_Type      ! Source type
      character(len=4)::calib_Type    ! Calibration type
      integer :: avg_or_Sample         ! Data was averaged (0) or sampled (1)
      integer :: poes_Signal           ! LAC, GAC, HRPT
      integer :: poes_up_down          ! POES ascending/descending
      integer :: cal_offset               ! needed for SEVIRI
      character(len=4) :: orig_Src_Type   ! Original source type of data
      integer, dimension(7) :: reserved    ! Reserved (6 calib pointer)
   end type area_header_type

contains

   logical function is_area_file ( file , sat_id_num )
      character ( len = * ) :: file
      integer, intent(out), optional :: sat_id_num
      integer :: lun
      integer :: area_status , area_version_num 
      
      lun = getlun()
      open (unit = lun , file = trim(file) , form = 'unformatted' &
               , access = 'direct' , recl = 4 , status ='old' , action ='read')
      read ( unit =lun, rec = 1) area_status
      read ( unit =lun, rec = 2) area_version_num 
     
      if ( area_version_num == 4) then           
         if (present (sat_id_num)) read ( unit =lun, rec = 3) sat_id_num
         is_area_file = .true.
        
      else
         is_area_file = .false.
      end if
      close ( unit = lun)
      
      is_area_file = .true.
   end function is_area_file 
   
   
   !
   !
   !
   subroutine read_area_header ( file , num_elem , num_pix )
      character ( len = * ) , intent(in):: file
      integer, intent(out), optional :: num_elem , num_pix
    
      print*,'hallo'
      num_elem =0
      num_pix = 0
      print*,file
   
   end subroutine read_area_header 
   

end module area_file_mod
