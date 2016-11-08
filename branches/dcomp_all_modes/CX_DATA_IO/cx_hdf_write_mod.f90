! $Id$
!
!

module cx_hdf_write_mod
   
   
   
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
	
	include 'hdf.f90'
   
   public  hdf_file_open 
   public copy_global_attributes
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
      integer(kind=1), intent(in) :: value
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
      integer(kind = 2), intent(in) :: value
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
   integer (kind=4) function create_sds_1d ( id, name, dim, dtype )
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
   integer (kind=4) function create_sds_2d ( id, name, dim, dtype )
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
         print*,'level-2 file error while closing ', istatus,id_file
         stop
      end if
   end subroutine close_file
   
   !------------------------------------------------------------------------------------------------
! SUBROUTINE Name: COPY_GLOBAL_ATTRIBUTES
!
! Function:
!    Copies the global attributes from one file and puts them in a different file
!-----------------------------------------------------------------------------------------------

   subroutine COPY_GLOBAL_ATTRIBUTES(File_in,file_out, exclude)
      
      character (len = * ) :: file_in
      character (len = * ) :: file_out
      character (len =100), optional :: exclude (:)
      integer:: Sd_Id_Input
      integer:: Sd_Id_Output
      integer:: Istatus
      integer:: num_global_attrs
      integer:: Num_Sds_Input
      integer:: iattr
      character(len=100):: Attr_Name
      integer:: Data_Type
      integer:: count
      integer:: Attr_Buffer_i4
      integer:: Attr_Buffer_i2
      integer:: Attr_Buffer_i1
      real:: Attr_Buffer_r4
      character(len=500):: Attr_Buffer_char

      integer:: sffinfo
      integer:: sfgainfo
      integer:: sfrnatt
      integer:: sfsnatt
   
   sd_id_input = hdf_file_open(trim(File_in))
   sd_id_output = hdf_file_open(trim(File_out))
   
   Istatus = sffinfo(Sd_Id_Input,Num_Sds_Input,num_global_attrs)

   do iattr = 0, num_global_attrs-1

       Istatus = sfgainfo(Sd_Id_Input,iattr,Attr_Name,Data_Type,count)

       !--- skip certain level2 attributes in level2b
       if (Attr_Name == "FILENAME" .or. &
           Attr_Name == "L1B" .or. &
           Attr_Name == "NUMBER_OF_SCANS_LEVEL1B" .or. &
           Attr_Name == "NUMBER_OF_SCANS_LEVEL2" .or. &
           Attr_Name == "START_YEAR" .or.  &
           Attr_Name == "END_YEAR" .or.  &
           Attr_Name == "START_DAY" .or.  &
           Attr_Name == "END_DAY" .or.  &
           Attr_Name == "START_TIME" .or.  &
           Attr_Name == "END_TIME") then
           cycle
        endif

       if (Data_Type == DFNT_INT8) then
           Istatus = sfrnatt(Sd_Id_Input,iattr,Attr_Buffer_i1)
           Istatus = sfsnatt(Sd_Id_Output,trim(Attr_Name),Data_Type,count,Attr_Buffer_i1)
       endif

       if (Data_Type == DFNT_INT16) then
           Istatus = sfrnatt(Sd_Id_Input,iattr,Attr_Buffer_i2)
           Istatus = sfsnatt(Sd_Id_Output,trim(Attr_Name),Data_Type,count,Attr_Buffer_i2)
       endif

       if (Data_Type == DFNT_INT32) then
           Istatus = sfrnatt(Sd_Id_Input,iattr,Attr_Buffer_i4)
           Istatus = sfsnatt(Sd_Id_Output,trim(Attr_Name),Data_Type,count,Attr_Buffer_i4)
       endif

       if (Data_Type == DFNT_FLOAT32) then
           Istatus = sfrnatt(Sd_Id_Input,iattr,Attr_Buffer_r4)
           Istatus = sfsnatt(Sd_Id_Output,trim(Attr_Name),Data_Type,count,Attr_Buffer_r4)
       endif

       if (Data_Type == DFNT_CHAR) then
           Istatus = sfrnatt(Sd_Id_Input,iattr,Attr_Buffer_char)
           Istatus = sfsnatt(Sd_Id_Output,trim(Attr_Name),Data_Type,count,trim(Attr_Buffer_char))
       endif
   enddo
   
   call close_file (sd_id_input)   
   call close_file (sd_id_output)   

   end subroutine 

   

end module cx_hdf_write_mod
