! $Id$
!
!
!
module data_write_tools
   implicit none
   

   type cr_data_type
      character ( len = 3 ) :: key
      character ( len = 50 ) :: name
      character ( len = 100) :: longname
      logical :: is_scaled
      integer :: offset (2)
      double precision:: slope 
      integer :: nx
      integer :: ny
      integer :: nx_full
      integer :: ny_full
      real    , allocatable :: data (:,:)
      integer , allocatable :: data_scaled (:,:)
      character (len=256) :: filename = 'test2.h5'
      
      integer :: h5_offset (2)
      integer :: h5_count  (2) 
      
      logical :: h5_created = .false.
        
   contains
      procedure :: init
      procedure :: get_property
      procedure :: set_property
      procedure :: update
      procedure :: write_to_h5
      procedure :: deallocate_data
   end type cr_data_type

contains

   !-------------------------------------------------------------------------------------
   !
   !-------------------------------------------------------------------------------------
   subroutine init ( self  , nx_full , ny_full )
      class ( cr_data_type ) :: self
      integer , intent (in) :: nx_full , ny_full
      
      self % h5_offset = [0,0]
      self % nx_full = nx_full  
      self % ny_full = ny_full
      self % h5_created = .false.
   
   end subroutine init
   
   !-------------------------------------------------------------------------------------
   !
   !-------------------------------------------------------------------------------------
   subroutine get_property ( self &
      , nx , ny , data , h5_offset )
      class ( cr_data_type ) :: self
      integer , intent(out) , optional :: nx
      integer , intent(out) , optional :: ny
      real , intent(out) , optional  :: data(:,:)
      integer , optional :: h5_offset ( 2)
      
      
      if ( present ( nx ) ) nx = self % nx
      if ( present ( ny ) ) ny = self % ny
      if ( present ( data ) ) data = self % data
      if ( present ( h5_offset) ) h5_offset = self % h5_offset
      
   end subroutine
   
   !-------------------------------------------------------------------------------------
   !
   !-------------------------------------------------------------------------------------
   subroutine set_property ( self &
      , nx , ny , data ,file , name , h5_offset )
      class ( cr_data_type ) :: self
      integer , intent(in) , optional :: nx
      integer , intent(in) , optional :: ny
      real    , intent(in) , optional :: data(:,:)
      character ( len = *) , optional :: file
      character ( len = *) , optional :: name
      integer , optional :: h5_offset ( 2)
      integer :: n_data (2)
      
      if ( present ( nx ) )  self % nx = nx
      if ( present ( ny ) )  self % ny = ny
      if ( present( file) ) self % filename = file
      if ( present( name ) ) self % name = name
      if ( present ( h5_offset )) self % h5_offset = h5_offset
      if ( present ( data ) )  then
         n_data = shape ( data )
         if (allocated (  self % data )) deallocate (  self % data )
         allocate ( self % data (n_data(1),n_data(2)), source = data )
         
         call self % update ()
      end if   
   end subroutine
   
   !-------------------------------------------------------------------------------------
   !
   !-------------------------------------------------------------------------------------   
   subroutine update ( self )
      class ( cr_data_type ) :: self   
      integer :: dim_data (2)
      dim_data = shape ( self % data)
      
      self % nx = dim_data (1)
      self % ny = dim_data (2)
   
   end subroutine update
   
   !-------------------------------------------------------------------------------------
   !
   !-------------------------------------------------------------------------------------      
   subroutine deallocate_data ( self )
      class ( cr_data_type ) :: self
      deallocate ( self % data )
   
   end subroutine deallocate_data
   
   !-------------------------------------------------------------------------------------
   !
   !-------------------------------------------------------------------------------------
   subroutine write_to_h5 ( self  )
#include "cx_h5_attr_add_char.interface"   
      
      use hdf5
      
      implicit none
      class ( cr_data_type ) :: self
      integer :: hdferr
      integer (HID_T)   :: file_id       ! File identifier
      integer (HSIZE_T) :: dims(1:2)
      integer (HSIZE_T) :: dims1(1:2)
      integer (HSIZE_T) :: dims_full(1:2)
      integer (HID_T)   :: file_h5 , dataspace , dset, memspace , cparms,filespace
      integer (HSIZE_T) :: data_dims(7)
      integer (HSIZE_T) :: data_dims_a(1)
      integer (HSIZE_T) :: count(1:2)
      integer (HSIZE_T) :: offset(1:2)
      integer (HSIZE_T) :: stride(1:2)
      integer (HSIZE_T) :: block(1:2)
      integer (HSIZE_T) :: max_dims(1:2)
      integer :: arank = 1
      

      ! INTEGER(HID_T) :: file_id       ! File identifier
      INTEGER(HID_T) :: dset_id       ! Dataset identifier
      INTEGER(HID_T) :: attr_id       ! Attribute identifier
      INTEGER(HID_T) :: aspace_id     ! Attribute Dataspace identifier
      INTEGER(HID_T) :: atype_id      ! Attribute Dataspace identifier
      INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/2/) ! Attribute dimension
      integer :: rank = 2
      INTEGER(SIZE_T) :: attrlen    ! Length of the attribute string

      CHARACTER(LEN=80), DIMENSION(2) ::  attr_data  ! Attribute data

      INTEGER     ::   error ! Error flag
      ! INTEGER(HSIZE_T), DIMENSION(1) :: data_dims

      !
      ! Initialize attribute's data
      !
      attr_data(1) = 'Dataset character attribute'
      attr_data(2) = "Some other string here     "
      attrlen = 80
      
      call h5open_f(hdferr)
    
      dims   = [ self%nx , self%ny ]
      block  = [ 1,1 ]
      stride = [ 1,1 ]
      
      data_dims(:) = 0
      data_dims(1) = dims(1)
      data_dims(2) = dims(2)
      
      if ( .not. self % h5_created ) then 
           
         call h5fopen_f (trim(self % filename), H5F_ACC_RDWR_F, file_id, hdferr)
        
         dims_full = [ self % nx_full, self % ny_full ]
         call h5screate_simple_f ( 2 , dims_full, dataspace ,hdferr )
         call h5dcreate_f(file_id,  self % name , H5T_NATIVE_REAL,dataspace,dset_id,hdferr)
         
         
         call cx_h5_attr_add_char ( dset_id, 'UNITS','km')
         
         !call cx_h5_attr_add_char ( file_id, 'UNITS','km')
         
         call h5sclose_f(dataspace, hdferr)
         call h5dclose_f(dset_id, hdferr)
         call h5fclose_f(file_id, hdferr)
         
         self % h5_created = .true.
              
      end if
             
      call h5fopen_f (trim(self % filename), H5F_ACC_RDWR_F, file_id, hdferr)
      call h5dopen_f(file_id, self % name , dset, hdferr)
      call h5dget_space_f(dset, dataspace, hdferr)
        
      offset = self % h5_offset
      count = dims
      
      
               
      call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
            offset, count , hdferr, stride, BLOCK) 
            
      call h5screate_simple_f(rank, dims, memspace, hdferr)  
        
      call h5dwrite_f (dset, H5T_NATIVE_REAL, self% data, data_dims, hdferr &
            & ,  memspace, dataspace) 
                       
      call h5sclose_f(dataspace, hdferr)
      call h5dclose_f(dset, hdferr)
      
      call h5fclose_f(file_id, hdferr)
      self % h5_created = .true.
      self % h5_offset(2) = offset(2) + dims(2)    
      
      
            
   end subroutine write_to_h5
   

end module data_write_tools


