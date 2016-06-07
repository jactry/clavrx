! $Id$
!
!

module level2_hdf_mod
   use constants
   use hdf
   
   implicit none
  
   integer :: file_id = -1
  
   interface add_att
      module procedure add_att_char, add_att_int &
         , add_att_real, add_att_int2
   end interface
   
   interface write_sds 
      module procedure write_sds_1d_dt1, write_sds_1d_dt3, write_sds_1d_dt4 &
               ,write_sds_2d_dt1, write_sds_2d_dt2, write_sds_2d_dt4 
   end interface
   
   type file_dims_type
      logical :: is_set = .false.
      integer :: dim1
      integer :: dim2
      integer :: dim3
      
      contains
      procedure:: set => file_dims__set
   
   end type file_dims_type
   
   type(file_dims_type) :: dims_file
   
   
   integer :: istatus
   
contains
   !
   !
   !
   subroutine file_dims__set (this, dim1, dim2, dim3)
      class(file_dims_type) :: this
      integer, intent(in) :: dim1, dim2, dim3
      
      this % dim1 = dim1
      this % dim2 = dim2
      this % dim3 = dim3
      this%is_set= .true.
      
   end subroutine file_dims__set
   
   !
   !
   !
   integer function file_open ( file)
      character(len=*) , intent(in) :: file
      integer :: dum
      integer :: sfstart
      
      dum = sfstart(trim(file),DFACC_CREATE)
      print*,'dum==> ',dum
      file_open = dum
   end function file_open
   
   !
   !
   !
   subroutine add_att_char (id,name, value)
      integer, intent(in) :: id
      character (len= *), intent(in) :: name
      character (len= *), intent(in) :: value
      integer :: sfscatt
      istatus = sfscatt ( id,name,DFNT_CHAR8,len_trim(value),trim(value))
   end subroutine add_att_char
   
   subroutine add_att_int (id, name, value)
      integer, intent(in) :: id
      character (len= *), intent(in) :: name
      integer, intent(in) :: value
      integer :: sfsnatt
      istatus = sfsnatt ( id,name,DFNT_INT32,1,value)   
   end subroutine add_att_int
   
   subroutine add_att_int2 (id, name, value)
      integer, intent(in) :: id
      character (len= *), intent(in) :: name
      integer(kind = int2), intent(in) :: value
      integer :: sfsnatt
      istatus = sfsnatt ( id,name,DFNT_INT32,1,value)   
   end subroutine add_att_int2
   
   subroutine add_att_real (id,name, value)
      integer, intent(in) :: id
      character (len= *), intent(in) :: name
      real, intent(in) :: value
      integer :: sfsnatt
      istatus = sfsnatt ( id,name,DFNT_FLOAT32,1,value)   
   end subroutine add_att_real
   
   !
   !
   !
   integer (kind=int4) function create_sds ( id, name, dim, dtype )
      integer, intent(in) :: id
      character (len=*), intent(in) :: name
      integer,intent(in) :: dim
      integer,intent(in) :: dtype
      integer :: dim_hdf(dim)
      integer :: dtype_hdf
      integer::sfcreate
      
      select case (dtype)
      case(1)
         dtype_hdf = DFNT_INT8
      case(2)
         dtype_hdf = DFNT_INT16
      case(3)
         dtype_hdf = DFNT_INT32
      case(4)
         dtype_hdf = DFNT_FLOAT32
      end select
      
      select case ( dim )
         case(1)
         dim_hdf = dims_file % dim1
         case(2)
         dim_hdf = (/ dims_file % dim1, dims_file % dim2 /)
      end select
      
      create_sds = sfcreate(id,name,dtype_hdf,dim, dim_hdf)
      
   end function create_sds
   
   !
   !
   !
   integer function compress_sds ( id, cmpr_flag, chunk)
      integer, intent(in) :: id
      integer, intent(in) :: cmpr_flag
      integer, intent(in) :: chunk(2)
      integer :: sfschnk
      
      integer :: comp_type
      integer :: comp_prm(2)
      
      if ( cmpr_flag == 0) return
      if ( cmpr_flag == 1 ) then
         comp_type = 4
         comp_prm(1) = 6
         comp_prm(2) = 0
      end if
      if ( cmpr_flag == 2 ) then
         comp_type = 5
         comp_prm(1) = 32
         comp_prm(2) = 2
      end if
      
      compress_sds = sfschnk(Id,chunk,Comp_Type,Comp_Prm)
   
   end function compress_sds
   
   !
   !
   !
   integer function write_sds_1d_dt1 ( sds_id, start, stride, edge , values )
      integer ::sds_id
      integer :: start
      integer :: stride
      integer :: edge
      integer (kind =1 ) :: values(:)  
      integer :: sfwdata
     
      write_sds_1d_dt1 = sfwdata ( sds_id,start,stride,edge,values)
   end function write_sds_1d_dt1
   
   integer function write_sds_1d_dt3 ( sds_id, start, stride, edge , values )
      integer ::sds_id
      integer :: start
      integer :: stride
      integer :: edge
      integer (kind = 4 ) :: values(:)  
      integer :: sfwdata
     
      write_sds_1d_dt3 = sfwdata ( sds_id,start,stride,edge,values)
    
   end function write_sds_1d_dt3
   
   integer function write_sds_1d_dt4 ( sds_id, start, stride, edge , values )
      integer ::sds_id
      integer :: start
      integer :: stride
      integer :: edge
      real :: values(:) 
      integer :: sfwdata 
     
      write_sds_1d_dt4 = sfwdata ( sds_id,start,stride,edge,values)
   end function write_sds_1d_dt4
   
   integer function write_sds_2d_dt1 ( sds_id, start, stride, edge , values )
      integer ::sds_id
      integer :: start(2)
      integer :: stride(2)
      integer :: edge(2)
      integer (kind =1 ) :: values(:,:) 
      integer :: sfwdata 
      
      write_sds_2d_dt1 = sfwdata ( sds_id,start,stride,edge,values)
   end function write_sds_2d_dt1
   
   integer function write_sds_2d_dt2 ( sds_id, start, stride, edge , values )
     integer ::sds_id
      integer :: start(2)
      integer :: stride(2)
      integer :: edge(2)
      integer (kind =2 ) :: values(:,:) 
      integer :: sfwdata 
      
      write_sds_2d_dt2 = sfwdata ( sds_id,start,stride,edge,values)
   end function write_sds_2d_dt2
   
   integer function write_sds_2d_dt4 ( sds_id, start, stride, edge , values )
       integer ::sds_id
      integer :: start(2)
      integer :: stride(2)
      integer :: edge(2)
      real  :: values(:,:)     
      integer :: sfwdata   
      write_sds_2d_dt4 = sfwdata ( sds_id,start,stride,edge,values)
   end function write_sds_2d_dt4
   
   
   !
   !
   !
   subroutine close_sds (id_sds)
      integer :: id_sds
      integer :: istatus
      integer :: sfendacc
      istatus = sfendacc ( id_sds )
      
   end subroutine close_sds
   
   subroutine close_file(id_file)
      integer :: id_file
      integer :: istatus
      integer :: sfend
      istatus = sfend ( id_file)
      if ( istatus /= 0 ) then
         print*,'level-2 file error whil closing ', istatus,id_file
         stop
      end if
   end subroutine close_file
   

end module level2_hdf_mod
