! $Id$
!
!  ----  ------
module cx_sds_io_mod

   use cx_hdf_read_mod, only: &
    hdf_get_finfo &
    , hdf_get_file_sds &
    , hdf_get_file_att
   
   use  cx_sds_type_definitions_mod, only: &
      cx_sds_type &
      , cx_att_type &
      , cx_sds_data_type
   
   implicit none
   private
   
   include 'cx_sds_constants.F90'
   
   integer, parameter, public :: MAXNCDIM = 32
   integer, parameter, public :: MAXNCNAM = 128
   
   interface cx_sds_read
      module procedure cx_sds_read_2d_real &
          , cx_sds_read_1d_real &
          , cx_sds_read_3d_real &
          , cx_sds_read_5d_real
   
   end interface
   
   
   
   
   public :: cx_sds_finfo
   public :: cx_sds_read_raw
   public :: cx_sds_read    
   public :: cx_sds_att
   public :: cx_att_int
   public :: cx_att_r4 
   contains
   

   ! ------------------------------------------------------------------------------
   !
   ! ------------------------------------------------------------------------------
   function cx_sds_finfo ( file , ftype, nsds, sds_name, natt, att_name )
      implicit none
      
      integer :: cx_sds_finfo
   
      integer, intent(out) :: ftype
      character(len=*), intent(in) :: file
      integer, intent(out) :: nsds
      integer, intent(out) :: natt                 
      character ( len = MAXNCNAM), intent(out), allocatable :: sds_name(:)
      character ( len = MAXNCNAM), intent(out), allocatable :: att_name(:)
      
      
      ! - inialize
      cx_sds_finfo = 1
      
      !- first check which kind of file
      ! --   
      ftype = 1 ! todo define params for HDF4, HDF5(2) and NCDF(3), not_existing(-1), not defined(0)
      
      cx_sds_finfo = hdf_get_finfo(file, nsds, sds_name, natt, att_name)
         
   
   end function cx_sds_finfo
   ! ------------------------------------------------------------------------------
   !
   ! ------------------------------------------------------------------------------
   function cx_sds_att (file, att_name, att )
      implicit none
      integer :: cx_sds_att
      character(len =*), intent(in) :: file
      character(len =*), intent(in) :: att_name
      type (cx_att_type),intent(out) :: att
      
      integer :: natt
      type(cx_att_type),  allocatable :: attrs(:)
      
      integer :: i

      cx_sds_att = hdf_get_file_att(trim(file), natt, attrs)
      
      do i = 1, natt 
        
         if (trim(att_name) .EQ. trim(attrs(i) % name )) then
            att = attrs(i)
         end if
      
      end do
      
   
   end function cx_sds_att
   
   
   ! This is a wrapper for integer output
   ! plaeas use this if you know that attribute is an integer
   function cx_att_int ( file, att_name, att_int)
      implicit none
      integer :: cx_att_int
      character(len =*), intent(in) :: file
      character(len =*), intent(in) :: att_name
      integer ,intent(out) :: att_int
      integer :: test
      type (cx_att_type) :: att
      
      test = cx_sds_att (file,att_name,att) 
      cx_att_int = -1
      att_int = -999
      if (allocated (att % data % i4values)) then
         att_int = att % data % i4values(1)
         cx_att_int = 0
      end if
          
      if (allocated (att % data % i2values)) then
         att_int = att % data % i2values(1)
         cx_att_int = 0
      end if    
      
   
   end function cx_att_int
   
   
   ! This is a wrapper for integer output
   ! plaeas use this if you know that attribute is an integer
   function cx_att_cha ( file, att_name, att_cha)
      implicit none
      integer :: cx_att_cha
      character(len =*), intent(in) :: file
      character(len =*), intent(in) :: att_name
      character ,intent(out) :: att_cha(:)
      integer :: test
      type (cx_att_type) :: att
      
      test = cx_sds_att (file,att_name,att) 
      cx_att_cha = -1
      att_cha = 'none'
      if (allocated (att % data % c1values)) then
         att_cha(1:att % data %dimsize(1)) = att % data % c1values(1:att % data %dimsize(1))
         cx_att_cha = 0
      end if
          
      
      
   
   end function cx_att_cha
   
   
     ! This is a wrapper for integer output
   ! plaeas use this if you know that attribute is an integer
   function cx_att_r4 ( file, att_name, att_r4)
      implicit none
      integer :: cx_att_r4
      character(len =*), intent(in) :: file
      character(len =*), intent(in) :: att_name
      real ,intent(out) :: att_r4
      integer :: test
      type (cx_att_type) :: att
      
      test = cx_sds_att (file,att_name,att) 
      cx_att_r4 = -1
      att_r4 = -999.
      if (allocated (att % data % r4values)) then
         att_r4 = att % data % r4values(1)
         cx_att_r4 = 0
      end if
    end function cx_att_r4      
   
   
   ! ------------------------------------------------------------------------------
   !
   ! ------------------------------------------------------------------------------
   function cx_sds_read_raw ( file, sds_name, sds )   
      integer :: cx_sds_read_raw
      character (len = * ), intent(in) :: file
      character (len = * ), intent(in) :: sds_name
      type (cx_sds_type),  intent(out), allocatable, target :: sds(:)
     
      integer :: nsds
      
      
      cx_sds_read_raw = -1
      
      sds.nattr = 9
      
     
      cx_sds_read_raw = hdf_get_file_sds(file, nsds,sds,1, (/sds_name/), .false.,.false.) 
      
      ! sds.copy_from_
   
   end function cx_sds_read_raw
   
   ! ------------------------------------------------------------------------------
   !
   ! ------------------------------------------------------------------------------     
   function cx_sds_read_1d_real ( file, sds_name, out )
      integer :: cx_sds_read_1d_real
      character (len = * ), intent(in) :: file
      character (len = * ), intent(in) :: sds_name 
      real, intent(out), allocatable :: out(:)
      real, allocatable:: temp_1d(:)
      type ( cx_sds_type), allocatable, target :: sds(:)
      type ( cx_sds_data_type), pointer :: pd
      type ( cx_sds_type), pointer :: ps
      integer :: test
      integer :: dim1, dim2
      real :: add_offset(1)
      real :: slope (1)
      real :: missing(1)
      real :: scaled(1)
      
      
      if (  cx_sds_read_raw ( file, sds_name, sds) < 0 ) goto 9999 
      
      pd=>sds(1).data
      ps=>sds(1)
      
      add_offset = ps.get_att('add_offset')
      slope = ps.get_att('scale_factor')
      missing = ps.get_att('SCALED_MISSING')
      scaled = ps.get_att('SCALED')
     
      
      
      allocate(temp_1d(pd.nval))
      call pd.transform_to_real(temp_1d)
      
      
       
      
      dim1 = pd%dimsize(1)
     
     
      allocate(out(dim1))
      
      out = reshape (temp_1d,(/dim1/))
      if (scaled(1) .EQ. 1) then
         out = out * slope(1) + add_offset(1) 
      end if
      
     
      where (reshape (temp_1d,(/dim1/)) .EQ. missing(1))
         out = -999.
      end where
      
      
      
      cx_sds_read_1d_real = 0
9999 continue
      !cx_sds_read_2d_real = -1 
   
   end function cx_sds_read_1d_real   
   
   
   function cx_sds_read_2d_real ( file, sds_name, out )
      integer :: cx_sds_read_2d_real
      character (len = * ), intent(in) :: file
      character (len = * ), intent(in) :: sds_name 
      real, intent(out), allocatable :: out(:,:)
      real, allocatable:: temp_1d(:)
      type ( cx_sds_type), allocatable, target :: sds(:)
      type ( cx_sds_data_type), pointer :: pd
      type ( cx_sds_type), pointer :: ps
      integer :: test
      integer :: dim1, dim2
      real :: add_offset(1)
      real :: slope (1)
      real :: missing(1)
      real :: scaled(1)
      
      
      if (  cx_sds_read_raw ( file, sds_name, sds) < 0 ) goto 9999 
      pd=>sds(1).data
      ps=>sds(1)
      
      add_offset = ps.get_att('add_offset')
      slope = ps.get_att('scale_factor')
      missing = ps.get_att('SCALED_MISSING')
      scaled = ps.get_att('SCALED')
      
     
    
     
     
      allocate(temp_1d(pd.nval))
      call pd.transform_to_real(temp_1d)
      
     
      
      dim1 = pd%dimsize(1)
      dim2 = pd%dimsize(2)
     
      allocate(out(dim1,dim2))
      
      out = reshape (temp_1d,(/dim1,dim2/))
      if (scaled(1) .EQ. 1) then
         out = out * slope(1) + add_offset(1) 
      end if
     
      where (reshape (temp_1d,(/dim1,dim2/)) .EQ. missing(1))
         out = -999.
      end where
      
     
      
      cx_sds_read_2d_real = 0
9999 continue
      !cx_sds_read_2d_real = -1 
   
   end function cx_sds_read_2d_real
   
   function cx_sds_read_3d_real ( file, sds_name, out )
      integer :: cx_sds_read_3d_real
      character (len = * ), intent(in) :: file
      character (len = * ), intent(in) :: sds_name 
      real, intent(out), allocatable :: out(:,:,:)
      real, allocatable:: temp_1d(:)
      type ( cx_sds_type), allocatable, target :: sds(:)
      type ( cx_sds_data_type), pointer :: pd
      type ( cx_sds_type), pointer :: ps
      integer :: test
      integer :: dim1, dim2, dim3
      real :: add_offset(1)
      real :: slope (1)
      real :: missing(1)
      real :: scaled(1)
      
      
      if (  cx_sds_read_raw ( file, sds_name, sds) < 0 ) goto 9999 
      pd=>sds(1).data
      ps=>sds(1)
      
      add_offset = ps.get_att('add_offset')
      slope = ps.get_att('scale_factor')
      missing = ps.get_att('SCALED_MISSING')
      scaled = ps.get_att('SCALED')
      
      
      
      allocate(temp_1d(pd.nval))
      call pd.transform_to_real(temp_1d)
      
      
      
      
      dim1 = pd%dimsize(1)
      dim2 = pd%dimsize(2)
      dim3 = pd%dimsize(3)
      allocate(out(dim1,dim2,dim3))
      
      out = reshape (temp_1d,(/dim1,dim2,dim3/))
      if (scaled(1) .EQ. 1) then
         out = out * slope(1) + add_offset(1) 
      end if
      
     
      where (reshape (temp_1d,(/dim1,dim2,dim3/)) .EQ. missing(1))
         out = -999.
      end where
      
     
      
      cx_sds_read_3d_real = 0
9999 continue
      !cx_sds_read_2d_real = -1 
   
   end function cx_sds_read_3d_real
   
   function cx_sds_read_5d_real ( file, sds_name, out )
      integer :: cx_sds_read_5d_real
      character (len = * ), intent(in) :: file
      character (len = * ), intent(in) :: sds_name 
      real, intent(out), allocatable :: out(:,:,:,:,:)
      real, allocatable:: temp_1d(:)
      type ( cx_sds_type), allocatable, target :: sds(:)
      type ( cx_sds_data_type), pointer :: pd
      type ( cx_sds_type), pointer :: ps
      integer :: test
      integer :: dim1, dim2, dim3, dim4 , dim5
      real :: add_offset(1)
      real :: slope (1)
      real :: missing(1)
      real :: scaled(1)
      
      
      if (  cx_sds_read_raw ( file, sds_name, sds) < 0 ) goto 9999 
      pd=>sds(1).data
      ps=>sds(1)
      
      add_offset = ps.get_att('add_offset')
      slope = ps.get_att('scale_factor')
      missing = ps.get_att('SCALED_MISSING')
      scaled = ps.get_att('SCALED')
      
      
      
      allocate(temp_1d(pd.nval))
      call pd.transform_to_real(temp_1d)
      
      
      
      
      dim1 = pd%dimsize(1)
      dim2 = pd%dimsize(2)
      dim3 = pd%dimsize(3)
      dim4 = pd%dimsize(4)
      dim5 = pd%dimsize(5)
     
      allocate(out(dim1,dim2,dim3,dim4,dim5))
      
      out = reshape (temp_1d,(/dim1,dim2,dim3,dim4,dim5/))
      if (scaled(1) .EQ. 1) then
         out = out * slope(1) + add_offset(1) 
      end if 
      
     
      where (reshape (temp_1d,(/dim1,dim2,dim3,dim4,dim5/)) .EQ. missing(1))
         out = -999.
      end where
      
     
      
      cx_sds_read_5d_real = 0
9999 continue
      !cx_sds_read_2d_real = -1 
   
   end function cx_sds_read_5d_real

end module
