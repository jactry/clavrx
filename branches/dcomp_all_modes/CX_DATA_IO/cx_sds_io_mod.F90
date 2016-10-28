! $Id$
!
!  ----  ------
module cx_sds_io_mod

   use cx_hdf_read_mod, only: &
    hdf_get_finfo &
    , hdf_get_file_sds 
   
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
      module procedure cx_sds_read_2d_real , cx_sds_read_1d_real, cx_sds_read_3d_real, cx_sds_read_5d_real
   
   end interface
   
   
   public :: cx_sds_finfo
   public :: cx_sds_read_raw
   public :: cx_sds_read    
   contains
   

   ! ------------------------------------------------------------------------------
   !
   ! ------------------------------------------------------------------------------
   function cx_sds_finfo ( file , ftype, nsds, sds_name, natt, att_name )
      integer :: cx_sds_finfo
      
      integer, intent(out) :: ftype
      character(len=*), intent(in) :: file
      integer, intent(out) :: nsds
      integer, intent(out) :: natt                 
      character ( len = MAXNCNAM), intent(out), allocatable :: sds_name(:)
      character ( len = MAXNCNAM), intent(out), allocatable :: att_name(:)
      
      
      cx_sds_finfo = 1
      !- first check which kind of file
      ! --   
      ftype = 1 ! todo define params for HDF4, HDF5(2) and NCDF(3), not_existing(-1), not defined(0)
      
      cx_sds_finfo = hdf_get_finfo(file, nsds, sds_name, natt, att_name)
         
   
   end function cx_sds_finfo
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
         !out = out * slope(1) + add_offset(1) 
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
