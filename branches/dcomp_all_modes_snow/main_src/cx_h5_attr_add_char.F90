
SUBROUTINE cx_h5_attr_add_char ( grp_id , aname , avalue )
   use hdf5

   INTEGER(HID_T) , intent(in) :: grp_id      
   character ( len = * ) :: aname
   character ( len = *), target :: avalue
!INTF_END   
   
      integer :: arank
      integer (hsize_t) , dimension(1) :: adims = 1
      integer ( hid_t) :: aspace_id
      integer :: hdferr
      integer ( HID_T) :: atype_id
      integer ( size_t) :: attrlen
      
      TYPE(C_PTR) :: f_ptr
      integer ( hid_t) :: attr_id
      INTEGER, DIMENSION(7) :: data_dims
      
      
      
       attrlen = len_trim(avalue)
      
      arank = 1
      CALL h5screate_simple_f(arank, adims, aspace_id, hdferr)
      CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, hdferr)
      CALL h5tset_size_f(atype_id, attrlen, hdferr)
      
      CALL h5acreate_f( grp_id , aname, atype_id, aspace_id, &
                  attr_id, hdferr)
      
      data_dims(1) = 2
      f_ptr = C_LOC(avalue)
      CALL h5awrite_f(attr_id, atype_id, f_ptr,  hdferr)
      CALL h5aclose_f(attr_id, hdferr) 
      

END SUBROUTINE cx_h5_attr_add_char
