! $Id$
!
!

module cx_hdf_write_mod
   use constants
   use hdf
   
   implicit none
  
   integer :: file_id = -1
  
   interface add_att
      module procedure add_att_char, add_att_int &
         , add_att_real, add_att_int1, add_att_int2
   end interface
   
   interface create_sds
      module procedure create_sds_1d, create_sds_2d
   end interface
   
   interface write_sds 
      module procedure write_sds_1d_dt1, write_sds_1d_dt3, write_sds_1d_dt4 &
               ,write_sds_2d_dt1, write_sds_2d_dt2, write_sds_2d_dt4 
   end interface
  
   
   integer :: istatus
   
   public  hdf_file_open 
contains
   !
   !
   !
   !
   !
   !
   integer function hdf_file_open ( file, create)
      character(len=*) , intent(in) :: file
      integer :: dum
      integer :: sfstart
      logical, optional, intent(in) :: create
      logical :: create_loc
      
      create_loc = .false.
      if (present(create)) create_loc = create
      
      if ( create_loc ) then
         dum = sfstart(trim(file),DFACC_CREATE)
      else
         dum = sfstart(trim(file),DFACC_WRITE)
      end if
      hdf_file_open = dum
   end function hdf_file_open
   
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
   
   subroutine add_att_int1 (id, name, value)
      integer, intent(in) :: id
      character (len= *), intent(in) :: name
      integer(kind=int1), intent(in) :: value
      integer :: sfsnatt
      
      istatus = sfsnatt ( id,name,DFNT_INT8,1,value)   
   end subroutine add_att_int1
   
   subroutine add_att_int (id, name, value)
      integer, intent(in) :: id
      character (len= *), intent(in) :: name
      integer, intent(in) :: value
      integer :: sfsnatt
      
      istatus = sfsnatt ( id,name,DFNT_INT16,1,value)   
   end subroutine add_att_int
   
   
   subroutine add_att_int2 (id, name, value)
      integer, intent(in) :: id
      character (len= *), intent(in) :: name
      integer(kind = int2), intent(in) :: value
      integer :: sfsnatt
      
      istatus = sfsnatt ( id,name,DFNT_INT16,1,value)   
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
   integer (kind=int4) function create_sds_1d ( id, name, dim, dtype )
      integer, intent(in) :: id
      character (len=*), intent(in) :: name
      integer,intent(in) :: dim
      integer,intent(in) :: dtype
      
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
      
   
      
      create_sds_1d = sfcreate(id,name,dtype_hdf,1, dim)
      
   end function create_sds_1d
   
   
      !
   !
   !
   integer (kind=int4) function create_sds_2d ( id, name, dim, dtype )
      integer, intent(in) :: id
      character (len=*), intent(in) :: name
      integer,intent(in) :: dim(:)
      integer,intent(in) :: dtype
      integer :: rank
     
      integer :: dtype_hdf
      integer::sfcreate
      
      rank = size(dim)
      
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
      
      
      create_sds_2d = sfcreate(id,name,dtype_hdf,rank, dim)
      
   end function create_sds_2d
   

   
   
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
      
      if ( istatus /= 0 ) then
         !print*,'level-2 sds closing warning: ', istatus, id_sds
        ! stop
      end if
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
   

end module cx_hdf_write_mod
