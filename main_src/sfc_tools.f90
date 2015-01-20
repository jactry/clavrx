! $Header: https://svn.ssec.wisc.edu/repos/aw_clavrx/trunk/sfc_tools.f90 4 2014-03-20 21:14:10Z awalther $

module sfc_tools
   integer, parameter, public:: int1 = selected_int_kind(1)
   integer, parameter, public:: int2 = selected_int_kind(3)
   integer, parameter, public:: int4 = selected_int_kind(8)
   integer, parameter, public:: int8 = selected_int_kind(10)
   integer, parameter, public:: real4 = selected_real_kind(6,37)
   integer, parameter, public:: real8 = selected_real_kind(15,307)
   integer, parameter, public:: ipre = real4
   
   
   interface read_hdf_sds
      module procedure  &
         read_hdf_sds_i1,  &
         read_hdf_sds_i2, &
         read_hdf_sds_i2_scaled, &
         read_hdf_sds_r4
   end interface
   character ( len =20):: EXE_PROMPT ='clavrx13> '
   
contains
   
   !
   !
   !
  
   
   !
   !
   !
   function read_hdf_global_attribute_float64(id, attr_name) result(attr)
      integer(kind=int4), intent(in) :: id
      character(len=*), intent(in) :: attr_name
  
      real(kind=real8), dimension(1) :: buffer  
      real(kind=real8) :: attr
      integer(kind=int4) :: istatus, attr_index
      integer :: sffattr, sfrnatt
  
      attr_index = sffattr(id, trim(attr_name))
      istatus = sfrnatt(id, attr_index, buffer)
      if (istatus /= 0) then
         print "(a,'attribute ',a,' reading error from id: ',i0)",'clavrx>>',trim(attr_name),id
         buffer(1) = missing_value_real8
      end if
  
      attr = buffer(1)
      return

   end function read_hdf_global_attribute_float64
   
   !
   !
   !
   subroutine read_hdf_sds_dimensions(id, sds_name, num_lat, num_lon)
      integer(kind=int4), intent(in) :: id
      character(len=*), intent(in) :: sds_name
      integer(kind=int4), intent(out) :: num_lat, num_lon
  
      integer(kind=int4) :: sds_id, istatus, sds_rank, sds_type, sds_nattr
      integer(kind=int4), dimension(2) :: sds_dims
      integer :: sfselect, sfn2index, sfginfo, sfendacc
  
      sds_id = sfselect(id, sfn2index(id,trim(sds_name)))
    
      istatus = sfginfo(sds_id, sds_name, sds_rank, sds_dims, sds_type, sds_nattr)
      if (istatus /= 0) then
         print "(a,'error reading ',a,' from sd_id: ',i0)",' CLAVRX >>  ' ,trim(sds_name),id    
         stop
      end if
 
      num_lon = sds_dims(1)
      num_lat = sds_dims(2)
  
      istatus = sfendacc(sds_id)
  
   end subroutine read_hdf_sds_dimensions
   
   
   !-----------------------------------------------------------------------------------------------
   ! Subroutine to read hdf data.
   !--------------------------------------------------------------------------------------------
   subroutine read_hdf_sds_i1(sd_id, sds_name, istart, istride, iedge, buffer)
      integer(kind=int4), intent(in) :: sd_id
      character(*), intent(in) :: sds_name
      integer(kind=int4), dimension(2), intent(in) :: istart, istride
      integer(kind=int4), dimension(2), intent(inout) :: iedge
      integer(kind=int1), dimension(:,:), intent(out), allocatable :: buffer
      integer(kind=int4) :: astatus, istatus, sds_id, sds_rank, sds_type, sds_nattr
      integer(kind=int4), dimension(2) :: sds_dims
    
      integer :: sfselect, sfn2index, sfrdata, sfendacc, sfginfo
  
      sds_id = sfselect(sd_id, sfn2index(sd_id,trim(sds_name)))
    
      istatus = sfginfo(sds_id, sds_name, sds_rank, sds_dims, sds_type, sds_nattr)
      if (istatus /= 0) then
         print "(a,'error reading1 ',a,' from sd_id: ',i0)",exe_prompt,trim(sds_name),sd_id
         stop
      end if
  
      if (iedge(1) < 0) iedge(1) = sds_dims(1)
      if (iedge(2) < 0) iedge(2) = sds_dims(2)
  
      iedge(1) = min((sds_dims(1) - istart(1)),iedge(1))
      iedge(2) = min((sds_dims(2) - istart(2)),iedge(2))
      
      allocate(buffer(iedge(1),iedge(2)),stat=astatus)
      if (astatus /= 0) then
         print "(a,'Not enough memory to allocate 2d buffer.')",EXE_PROMPT
         stop
      end if
      
      istatus = sfrdata(sds_id, istart, istride, iedge, buffer)
      if (istatus /= 0) then
         print "(a,'error reading2 ',a,' from sd_id: ',i0)",exe_prompt,trim(sds_name),sd_id
         stop
      end if
      istatus = sfendacc(sds_id)
    
   end subroutine read_hdf_sds_i1
!
!
!
   subroutine read_hdf_sds_i2(sd_id, sds_name, istart, istride, iedge, buffer)
      integer(kind=int4), intent(in) :: sd_id
      character(*), intent(in) :: sds_name
      integer(kind=int4), dimension(2), intent(in) :: istart, istride
      integer(kind=int4), dimension(2), intent(inout) :: iedge
      integer(kind=int2), dimension(:,:), allocatable, intent(out) :: buffer
      integer(kind=int4) :: astatus, istatus, sds_id, sds_rank, sds_type, sds_nattr
      integer(kind=int4), dimension(2) :: sds_dims
    
      integer :: sfselect, sfn2index, sfrdata, sfendacc, sfginfo
  
      sds_id = sfselect(sd_id, sfn2index(sd_id,trim(sds_name)))
    
      istatus = sfginfo(sds_id, sds_name, sds_rank, sds_dims, sds_type, sds_nattr)
      if (istatus /= 0) then
         print "(a,'error reading ',a,' from sd_id: ',i0)",exe_prompt,trim(sds_name),sd_id
         stop
      end if
  
      if (iedge(1) < 0) iedge(1) = sds_dims(1)
      if (iedge(2) < 0) iedge(2) = sds_dims(2)
  
      iedge(1) = min((sds_dims(1) - istart(1)),iedge(1))
      iedge(2) = min((sds_dims(2) - istart(2)),iedge(2))
    
      allocate(buffer(iedge(1),iedge(2)),stat=astatus)
      if (astatus /= 0) then
         print "(a,'not enough memory to allocate 2d buffer.')",exe_prompt
         stop
      endif
    
      istatus = sfrdata(sds_id, istart, istride, iedge, buffer)
      if (istatus /= 0) then
         print "(a,'error reading ',a,' from sd_id: ',i0)",exe_prompt,trim(sds_name),sd_id
         stop
      end if
      istatus = sfendacc(sds_id)
    
   end subroutine read_hdf_sds_i2


!
!
!
   SUBROUTINE read_hdf_sds_i2_scaled(sd_id, sds_name, istart, istride, iedge,name_scale, name_offset, buffer_r4)
	   implicit none
	
      INTEGER(kind=int4), intent(in) :: sd_id
      CHARACTER(*), intent(in) :: sds_name
      INTEGER(kind=int4), dimension(2), intent(in) :: istart, istride
      INTEGER(kind=int4), dimension(2), intent(inout) :: iedge
      character ( len =30 ), intent(in) :: name_scale, name_offset 
      real       (kind=real4), dimension(:,:), allocatable, intent(out) :: buffer_r4
      INTEGER(kind=int2), dimension(:,:), allocatable:: buffer
      INTEGER(kind=int4) :: astatus, istatus, sds_id, sds_rank, sds_type, sds_nattr
      INTEGER(kind=int4), dimension(2) :: sds_dims
      REAL(kind=real4), dimension(1) :: scale_fac, offset
      integer :: sfselect, sfn2index, sfrdata, sfendacc, sfginfo
      integer :: attr_index , sffattr , sfrnatt
  
  
     
      sds_id = sfselect(sd_id, sfn2index(sd_id,trim(sds_name)))
      istatus = sfginfo(sds_id, sds_name, sds_rank, sds_dims, sds_type, sds_nattr)
      if (istatus /= 0) then
         print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
         stop
      endif
  
      attr_index = sffattr(sds_id, "scale_factor")
      istatus = sfrnatt(sds_id, attr_index, scale_fac)
      if (istatus /= 0) then
         print "(a,'Attribute (scale_factor) reading error ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
         stop
      endif
  
      attr_index = sffattr(sds_id, "add_offset")
      istatus = sfrnatt(sds_id, attr_index, offset)
      if (istatus /= 0) then
         print "(a,'Attribute (add_offset) reading error ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
         stop
      endif
      istatus = sfginfo(sds_id, sds_name, sds_rank, sds_dims, sds_type, sds_nattr)
      if (istatus /= 0) then
         print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
         stop
      endif
  
      if (iedge(1) < 0) iedge(1) = sds_dims(1)
      if (iedge(2) < 0) iedge(2) = sds_dims(2)
  
      iedge(1) = min((sds_dims(1) - istart(1)),iedge(1))
      iedge(2) = min((sds_dims(2) - istart(2)),iedge(2))
    
      allocate(buffer(iedge(1),iedge(2)),buffer_r4(iedge(1),iedge(2)) ,stat=astatus)
      if (astatus /= 0) then
         print "(a,'Not enough memory to allocate 2d buffer.')",EXE_PROMPT
         stop
      endif
    
      istatus = sfrdata(sds_id, istart, istride, iedge, buffer)
      buffer_r4 = buffer*scale_fac(1) + offset(1)
  
      if (istatus /= 0) then
         print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
         stop
      endif
      istatus = sfendacc(sds_id)
      
    
   END SUBROUTINE read_hdf_sds_i2_scaled
!
!
!
   SUBROUTINE read_hdf_sds_r4(sd_id, sds_name, istart, istride, iedge, buffer)
      INTEGER(kind=int4), intent(in) :: sd_id
      CHARACTER(*), intent(in) :: sds_name
      INTEGER(kind=int4), dimension(2), intent(in) :: istart, istride
      INTEGER(kind=int4), dimension(2), intent(inout) :: iedge
      real (kind=real4), dimension(:,:), allocatable, intent(out) :: buffer
      INTEGER(kind=int4) :: astatus, istatus, sds_id, sds_rank, sds_type, sds_nattr
      INTEGER(kind=int4), dimension(2) :: sds_dims
    
      integer :: sfselect, sfn2index, sfrdata, sfendacc, sfginfo
  
  
      
      sds_id = sfselect(sd_id, sfn2index(sd_id,trim(sds_name)))
    
      istatus = sfginfo(sds_id, sds_name, sds_rank, sds_dims, sds_type, sds_nattr)
      if (istatus /= 0) then
         print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
         stop
      endif
  
      if (iedge(1) < 0) iedge(1) = sds_dims(1)
      if (iedge(2) < 0) iedge(2) = sds_dims(2)
  
      iedge(1) = min((sds_dims(1) - istart(1)),iedge(1))
      iedge(2) = min((sds_dims(2) - istart(2)),iedge(2))
    
      allocate(buffer(iedge(1),iedge(2)),stat=astatus)
      if (astatus /= 0) then
         print "(a,'Not enough memory to allocate 2d buffer.')",EXE_PROMPT
         stop
      endif
    
      istatus = sfrdata(sds_id, istart, istride, iedge, buffer)
      if (istatus /= 0) then
         print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
         stop
      end if
      istatus = sfendacc(sds_id)
    
   END SUBROUTINE read_hdf_sds_r4
end module  sfc_tools
