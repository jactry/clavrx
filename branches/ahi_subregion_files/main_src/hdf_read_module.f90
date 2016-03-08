!$Id$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: hdf_read_module.f90 (src)
!       hdf_read_module(program)
!
! PURPOSE: 
!
! DESCRIPTION: 
!
! AUTHORS:
!
! COPYRIGHT
!   Copyright (C) 2006  Michael J. Pavolonis
!   National Oceanic and Atmospheric Administration
!
! GEOCAT
!
! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
! 02110-1301, USA.
!--------------------------------------------------------------------------------------

MODULE hdf_read_module

  USE CONSTANTS
  USE HDF
  !USE FILE_UTILITY
  !USE COMPARE_FLOAT_NUMBERS
  !USE MESSAGE_HANDLER

  implicit none
  
  interface hdf_sds_reader
    module procedure read_hdf_sds_int8_1d, &
                     read_hdf_sds_int16_1d, &
                     read_hdf_sds_int32_1d, &
                     read_hdf_sds_float32_1d, &
                     read_hdf_sds_float64_1d, &
                     read_hdf_sds_int8_2d, &
                     read_hdf_sds_int16_2d, &
                     read_hdf_sds_int32_2d, &
     read_hdf_sds_float32_2d, &
     read_hdf_sds_float64_2d, &
     read_hdf_sds_int8_3d, &
                     read_hdf_sds_int16_3d, &
                     read_hdf_sds_int32_3d, &
     read_hdf_sds_float32_3d, &
     read_hdf_sds_float64_3d, &
     read_hdf_sds_int8_4d, &
                     read_hdf_sds_int16_4d, &
                     read_hdf_sds_int32_4d, &
     read_hdf_sds_float32_4d, &
     read_hdf_sds_float64_4d, &
                     read_hdf_sds_int8_5d, &
                     read_hdf_sds_int16_5d, &
                     read_hdf_sds_int32_5d, &
     read_hdf_sds_float32_5d, &
     read_hdf_sds_float64_5d, &
                     read_hdf_sds_float64_6d
  end interface
  
  interface hdf_sds_attribute_reader
    module procedure read_hdf_attribute_char8_scalar,   &
                     read_hdf_attribute_int8_scalar,    &
                     read_hdf_attribute_int16_scalar,   &
     read_hdf_attribute_int32_scalar,   &
     read_hdf_attribute_float32_scalar, &
     read_hdf_attribute_float64_scalar, &
                     read_hdf_attribute_int8_vector,    &
                     read_hdf_attribute_int16_vector,   &
                     read_hdf_attribute_int32_vector,   &
                     read_hdf_attribute_float32_vector, &
                     read_hdf_attribute_float64_vector
  end interface
  
  
  interface unscale_cat
    module procedure unscale_cat_int8_2d, &
                     unscale_cat_int16_2d
  end interface
  
  CONTAINS

  
!---------------------------------------------------------------------
! This routine is used to read navigation and other constants fields
! from a GEOCAT NAV file.
!---------------------------------------------------------------------

SUBROUTINE read_navigation_file(file_id, &
                                xstart, &
ystart, &
xsize, &
ysize, &
xstride, &
        lat, &
lon, &
satzen, &
sataz, &
zsfc, &
space_mask, &
sfc_type, &
eco_type, &
land_mask, &
coast_mask, &
volcano_mask, &
desert_mask)
  INTEGER(kind=int4), intent(in) :: file_id
  INTEGER(kind=int4), intent(in) :: xstart, ystart
  INTEGER(kind=int4), intent(in) :: xsize, ysize
  INTEGER(kind=int4), intent(in) :: xstride
  REAL(kind=real4), dimension(:,:), intent(inout) :: lat
  REAL(kind=real4), dimension(:,:), intent(inout) :: lon
  REAL(kind=real4), dimension(:,:), intent(inout) :: satzen
  REAL(kind=real4), dimension(:,:), intent(inout) :: sataz
  REAL(kind=real4), dimension(:,:), intent(inout) :: zsfc
  INTEGER(kind=int1), dimension(:,:), intent(inout) :: space_mask
  INTEGER(kind=int1), dimension(:,:), intent(inout) :: sfc_type
  INTEGER(kind=int1), dimension(:,:), intent(inout) :: eco_type
  INTEGER(kind=int1), dimension(:,:), intent(inout) :: land_mask
  INTEGER(kind=int1), dimension(:,:), intent(inout) :: coast_mask
  INTEGER(kind=int1), dimension(:,:), intent(inout) :: volcano_mask
  INTEGER(kind=int1), dimension(:,:), intent(inout) :: desert_mask
  
  INTEGER(kind=int4), dimension(2) :: istart, istride, iedge
  INTEGER(kind=int4) :: error_status
  INTEGER(kind=int1), dimension(:,:), allocatable :: int1_buffer
  REAL(kind=real4), dimension(:,:), allocatable :: real4_buffer
  
  istart = (/xstart-1,ystart-1/)
  iedge = (/xsize, ysize/)
  istride = (/xstride, 1/) 

  error_status = hdf_sds_reader(file_id, "pixel_latitude", istart, istride, iedge, real4_buffer)
  lat(1:xsize,1:ysize) = real4_buffer
  if (error_status == sym%FAILURE) stop
  
  error_status = hdf_sds_reader(file_id, "pixel_longitude", istart, istride, iedge, real4_buffer)
  lon(1:xsize,1:ysize) = real4_buffer
  if (error_status == sym%FAILURE) stop
  
  error_status = hdf_sds_reader(file_id, "pixel_satellite_zenith_angle", istart, istride, iedge, real4_buffer)
  satzen(1:xsize,1:ysize) = real4_buffer
  if (error_status == sym%FAILURE) stop
  
  error_status = hdf_sds_reader(file_id, "pixel_satellite_azimuth_angle", istart, istride, iedge, real4_buffer)
  sataz(1:xsize,1:ysize) = real4_buffer
  if (error_status == sym%FAILURE) stop
  
  error_status = hdf_sds_reader(file_id, "pixel_surface_elevation", istart, istride, iedge, real4_buffer)
  zsfc(1:xsize,1:ysize) = real4_buffer
  if (error_status == sym%FAILURE) stop
  
  error_status = hdf_sds_reader(file_id, "pixel_space_mask", istart, istride, iedge, int1_buffer)
  space_mask(1:xsize,1:ysize) = int1_buffer
  if (error_status == sym%FAILURE) stop
  
  error_status = hdf_sds_reader(file_id, "pixel_surface_type", istart, istride, iedge, int1_buffer)
  sfc_type(1:xsize,1:ysize) = int1_buffer
  if (error_status == sym%FAILURE) stop
  
  error_status = hdf_sds_reader(file_id, "pixel_ecosystem_type", istart, istride, iedge, int1_buffer)
  eco_type(1:xsize,1:ysize) = int1_buffer
  if (error_status == sym%FAILURE) stop
  
  error_status = hdf_sds_reader(file_id, "pixel_land_mask", istart, istride, iedge, int1_buffer)
  land_mask(1:xsize,1:ysize) = int1_buffer
  if (error_status == sym%FAILURE) stop
  
  error_status = hdf_sds_reader(file_id, "pixel_coast_mask", istart, istride, iedge, int1_buffer)
  coast_mask(1:xsize,1:ysize) = int1_buffer
  if (error_status == sym%FAILURE) stop
  
  error_status = hdf_sds_reader(file_id, "pixel_volcano_mask", istart, istride, iedge, int1_buffer)
  volcano_mask(1:xsize,1:ysize) = int1_buffer
  if (error_status == sym%FAILURE) stop
  
  error_status = hdf_sds_reader(file_id, "pixel_desert_mask", istart, istride, iedge, int1_buffer)
  desert_mask(1:xsize,1:ysize) = int1_buffer
  if (error_status == sym%FAILURE) stop
  
  deallocate(int1_buffer, real4_buffer)

END SUBROUTINE read_navigation_file

!---------------------------------------------------------------------
! This routine is used to unscale int8 2d GEOCAT format data.
!---------------------------------------------------------------------

SUBROUTINE unscale_cat_int8_2d(fillval, factor, offset, &
                               nx, ny, scaled_data, unscaled_data)
  INTEGER(kind=int1), intent(in) :: fillval
  REAL(kind=real8), intent(in) :: factor
  REAL(kind=real8), intent(in) :: offset
  INTEGER(kind=int4), intent(in) :: nx, ny
  INTEGER(kind=int1), dimension(:,:), intent(in) :: scaled_data
  REAL(kind=real4), dimension(:,:), intent(out) :: unscaled_data
  
  INTEGER(kind=int4) :: ielem, iline
  
  do iline = 1, ny
    do ielem = 1, nx
      if (scaled_data(ielem,iline) /= fillval) then
        unscaled_data(ielem,iline) = scaled_data(ielem,iline)*factor + offset
      else
        unscaled_data(ielem,iline) = missing_value_real4
      endif
    end do
  end do
  
END SUBROUTINE unscale_cat_int8_2d

!---------------------------------------------------------------------
! This routine is used to unscale int8 2d GEOCAT format data.
!---------------------------------------------------------------------

SUBROUTINE unscale_cat_int16_2d(fillval, factor, offset, &
                                nx, ny, scaled_data, unscaled_data)
  INTEGER(kind=int2), intent(in) :: fillval
  REAL(kind=real8), intent(in) :: factor
  REAL(kind=real8), intent(in) :: offset
  INTEGER(kind=int4), intent(in) :: nx, ny
  INTEGER(kind=int2), dimension(:,:), intent(in) :: scaled_data
  REAL(kind=real4), dimension(:,:), intent(out) :: unscaled_data
  
  INTEGER(kind=int4) :: ielem, iline
  
  do iline = 1, ny
    do ielem = 1, nx
      if (scaled_data(ielem,iline) /= fillval) then
        unscaled_data(ielem,iline) = scaled_data(ielem,iline)*factor + offset
      else
        unscaled_data(ielem,iline) = missing_value_real4
      endif
    end do
  end do
  
END SUBROUTINE unscale_cat_int16_2d
  
!---------------------------------------------------------------------
! This routine is used to open an HDF file.
!---------------------------------------------------------------------

FUNCTION open_file_hdf_read(filename,file_id) result(error_status)
  CHARACTER(*), intent(in) :: filename  
  INTEGER(kind=int4), intent(out) :: file_id  
  
  INTEGER(kind=int4) :: error_status
  INTEGER :: sfstart
  
  error_status = sym%SUCCESS
  
  file_id = sfstart(trim(filename),DFACC_READ)
  if (file_id == FAIL) then
    print "(a,'Cannot Open HDF file, ',a)",&
      EXE_PROMPT,trim(filename)
    error_status = sym%FAILURE
  endif
  
  return

END FUNCTION open_file_hdf_read

!---------------------------------------------------------------------
! This routine is used to close an HDF file.
!---------------------------------------------------------------------

SUBROUTINE close_file_hdf_read(file_id,filename)
  INTEGER(kind=int4), intent(in) :: file_id
  CHARACTER(*), intent(in) :: filename
  
  INTEGER :: istatus
  INTEGER :: sfend
  
  istatus = sfend(file_id)
  if (istatus == FAIL) then
    print "(a,'Cannot close HDF file, ',a,' - aborting')",EXE_PROMPT,trim(filename)
    stop
  endif

END SUBROUTINE close_file_hdf_read

!-------------------------------------------------------------------
! This routine is used to read char8 HDF SDS attributes.
!-------------------------------------------------------------------

FUNCTION read_hdf_attribute_char8_scalar(id, sds_name, attr_name, attr, type) result(error_status)
  INTEGER(kind=int4), intent(in) :: id
  CHARACTER(len=*), intent(in) :: sds_name, attr_name
  CHARACTER(len=*), intent(inout) :: attr
  CHARACTER(len=1000), dimension(1) :: buffer
  INTEGER(kind=int4), intent(out), optional :: type
    
  CHARACTER(len=1020) :: name
  INTEGER(kind=int4) :: istatus, attr_index, sds_type, count, sds_id, error_status
  INTEGER :: sffattr, sfrcatt, sfgainfo, sfn2index, sfselect, sfendacc
  
  error_status = sym%SUCCESS
  
  sds_id = sfselect(id, sfn2index(id,trim(sds_name)))
  if (sds_id == FAIL) sds_id = id
  
  attr_index = sffattr(sds_id, trim(attr_name))
  if (attr_index == FAIL) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(attr_name),sds_id
    error_status = sym%FAILURE
    return
  endif
  
  istatus = sfgainfo(sds_id, attr_index, name, sds_type, count)
  if (istatus /= 0) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(attr_name),sds_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (present(type)) then
    type = sds_type
    istatus = sfendacc(sds_id)
    error_status = sym%SUCCESS
    return
  endif
  
  if (sds_type /= DFNT_INT8 .and. &
      sds_type /= DFNT_CHAR8 .and. &
      sds_type /= DFNT_CHAR) then
    print "(a,'Error reading attribute (type mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(attr_name),sds_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
    
  istatus = sfrcatt(sds_id, attr_index, buffer)
  if (istatus /= 0) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(attr_name),sds_id
    attr = " "
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  attr = TRIM(buffer(1)(1:count))

  istatus = sfendacc(sds_id)
  
  return

END FUNCTION read_hdf_attribute_char8_scalar

!-------------------------------------------------------------------
! This routine is used to read int8 HDF SDS attributes.
!-------------------------------------------------------------------

FUNCTION read_hdf_attribute_int8_scalar(id, sds_name, attr_name, attr, type) result(error_status)
  INTEGER(kind=int4), intent(in) :: id
  CHARACTER(len=*), intent(in) :: sds_name, attr_name
  INTEGER(kind=int1), intent(inout) :: attr
  INTEGER(kind=int4), intent(out), optional :: type
  
  INTEGER(kind=int1), dimension(1) :: buffer  
  CHARACTER(len=1020) :: name
  INTEGER(kind=int4) :: istatus, attr_index, sds_type, count, sds_id, error_status
  INTEGER :: sffattr, sfrnatt, sfgainfo, sfn2index, sfselect, sfendacc
  
  error_status = sym%SUCCESS
  attr = missing_value_int1
  
  sds_id = sfselect(id, sfn2index(id,trim(sds_name)))
  if (sds_id == FAIL) sds_id = id
  
  attr_index = sffattr(sds_id, trim(attr_name))
  if (attr_index == FAIL) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(attr_name),sds_id
    error_status = sym%FAILURE
    return
  endif
  
  istatus = sfgainfo(sds_id, attr_index, name, sds_type, count)
  if (istatus /= 0) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(attr_name),sds_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (present(type)) then
    type = sds_type
    istatus = sfendacc(sds_id)
    error_status = sym%SUCCESS
    return
  endif
  
  if (sds_type /= DFNT_INT8 .and. &
      sds_type /= DFNT_CHAR8 .and. &
      sds_type /= DFNT_CHAR) then
    print "(a,'Error reading attribute (type mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(attr_name),sds_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (count /= 1) then
    print "(a,'Error reading attribute (count mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(attr_name),sds_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
    
  istatus = sfrnatt(sds_id, attr_index, buffer)
  if (istatus /= 0) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(attr_name),sds_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  istatus = sfendacc(sds_id)
  
  attr = buffer(1)
  
  return

END FUNCTION read_hdf_attribute_int8_scalar

!-------------------------------------------------------------------
! This routine is used to read int16 HDF SDS attributes.
!-------------------------------------------------------------------

FUNCTION read_hdf_attribute_int16_scalar(id, sds_name, attr_name, attr, type) result(error_status)
  INTEGER(kind=int4), intent(in) :: id
  CHARACTER(len=*), intent(in) :: sds_name, attr_name
  INTEGER(kind=int2), intent(inout) :: attr
  INTEGER(kind=int4), intent(out), optional :: type
  
  INTEGER(kind=int2), dimension(1) :: buffer  
  CHARACTER(len=1020) :: name
  INTEGER(kind=int4) :: istatus, attr_index, sds_type, count, sds_id, error_status
  INTEGER :: sffattr, sfrnatt, sfgainfo, sfn2index, sfselect, sfendacc
  
  error_status = sym%SUCCESS
  attr = missing_value_int2
  
  sds_id = sfselect(id, sfn2index(id,trim(sds_name)))
  if (sds_id == FAIL) sds_id = id
  
  attr_index = sffattr(sds_id, trim(attr_name))
  if (attr_index == FAIL) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(attr_name),sds_id
    error_status = sym%FAILURE
    return
  endif
  
  istatus = sfgainfo(sds_id, attr_index, name, sds_type, count)
  if (istatus /= 0) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(attr_name),sds_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (present(type)) then
    type = sds_type
    istatus = sfendacc(sds_id)
    error_status = sym%SUCCESS
    return
  endif
  
  if (sds_type /= DFNT_INT16) then
    print*,trim(sds_name),type,DFNT_INT16
    print "(a,'Error reading attribute (type mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(attr_name),sds_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (count /= 1) then
    print "(a,'Error reading attribute (count mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(attr_name),sds_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
    
  istatus = sfrnatt(sds_id, attr_index, buffer)
  if (istatus /= 0) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(attr_name),sds_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  istatus = sfendacc(sds_id)
  
  attr = buffer(1)
  
  return
  
END FUNCTION read_hdf_attribute_int16_scalar

!-------------------------------------------------------------------
! This routine is used to read int32 HDF SDS attributes.
!-------------------------------------------------------------------

FUNCTION read_hdf_attribute_int32_scalar(id, sds_name, attr_name, attr, type) result(error_status)
  INTEGER(kind=int4), intent(in) :: id
  CHARACTER(len=*), intent(in) :: sds_name, attr_name
  INTEGER(kind=int4), intent(inout) :: attr
  INTEGER(kind=int4), intent(out), optional :: type
  
  INTEGER(kind=int4), dimension(1) :: buffer  
  CHARACTER(len=1020) :: name
  INTEGER(kind=int4) :: istatus, attr_index, sds_type, count, sds_id, error_status
  INTEGER :: sffattr, sfrnatt, sfgainfo, sfn2index, sfselect, sfendacc
  
  error_status = sym%SUCCESS
  attr = missing_value_int4
  
  sds_id = sfselect(id, sfn2index(id,trim(sds_name)))
  if (sds_id == FAIL) sds_id = id
  
  attr_index = sffattr(sds_id, trim(attr_name))
  if (attr_index == FAIL) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(attr_name),sds_id
    error_status = sym%FAILURE
    return
  endif
  
  istatus = sfgainfo(sds_id, attr_index, name, sds_type, count)
  if (istatus /= 0) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(attr_name),sds_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (present(type)) then
    type = sds_type
    istatus = sfendacc(sds_id)
    error_status = sym%SUCCESS
    return
  endif
  
  if (sds_type /= DFNT_INT32) then
    print "(a,'Error reading attribute (type mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(attr_name),sds_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (count /= 1) then
    print "(a,'Error reading attribute (count mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(attr_name),sds_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
    
  istatus = sfrnatt(sds_id, attr_index, buffer)
  if (istatus /= 0) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(attr_name),sds_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  istatus = sfendacc(sds_id)
  
  attr = buffer(1)
  
  return

END FUNCTION read_hdf_attribute_int32_scalar

!-------------------------------------------------------------------
! This routine is used to read float32 HDF SDS attributes.
!-------------------------------------------------------------------

FUNCTION read_hdf_attribute_float32_scalar(id, sds_name, attr_name, attr, type) result(error_status)
  INTEGER(kind=int4), intent(in) :: id
  CHARACTER(len=*), intent(in) :: sds_name, attr_name
  REAL(kind=real4), intent(inout) :: attr
  INTEGER(kind=int4), intent(out), optional :: type
  
  REAL(kind=real4), dimension(1) :: buffer  
  CHARACTER(len=1020) :: name
  INTEGER(kind=int4) :: istatus, attr_index, sds_type, count, sds_id, error_status
  INTEGER :: sffattr, sfrnatt, sfgainfo, sfn2index, sfselect, sfendacc
  
  error_status = sym%SUCCESS
  attr = missing_value_real4
  
  sds_id = sfselect(id, sfn2index(id,trim(sds_name)))
  if (sds_id == FAIL) sds_id = id
  
  attr_index = sffattr(sds_id, trim(attr_name))
  if (attr_index == FAIL) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(attr_name),sds_id
    error_status = sym%FAILURE
    return
  endif
  
  istatus = sfgainfo(sds_id, attr_index, name, sds_type, count)
  if (istatus /= 0) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(attr_name),sds_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (present(type)) then
    type = sds_type
    istatus = sfendacc(sds_id)
    error_status = sym%SUCCESS
    return
  endif
  
  if (sds_type /= DFNT_FLOAT32) then
    print "(a,'Error reading attribute (type mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(attr_name),sds_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (count /= 1) then
    print "(a,'Error reading attribute (count mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(attr_name),sds_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
    
  istatus = sfrnatt(sds_id, attr_index, buffer)
  if (istatus /= 0) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(attr_name),sds_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  istatus = sfendacc(sds_id)
  
  attr = buffer(1)
  
  return

END FUNCTION read_hdf_attribute_float32_scalar

!-------------------------------------------------------------------
! This routine is used to read float64 HDF SDS attributes.
!-------------------------------------------------------------------

FUNCTION read_hdf_attribute_float64_scalar(id, sds_name, attr_name, attr, type) result(error_status)
  INTEGER(kind=int4), intent(in) :: id
  CHARACTER(len=*), intent(in) :: sds_name, attr_name
  REAL(kind=real8), intent(inout) :: attr
  INTEGER(kind=int4), intent(out), optional :: type
  
  REAL(kind=real8), dimension(1) :: buffer  
  CHARACTER(len=1020) :: name
  INTEGER(kind=int4) :: istatus, attr_index, sds_type, count, sds_id, error_status
  INTEGER :: sffattr, sfrnatt, sfgainfo, sfn2index, sfselect, sfendacc
  
  error_status = sym%SUCCESS
  attr = missing_value_real8
  
  sds_id = sfselect(id, sfn2index(id,trim(sds_name)))
  if (sds_id == FAIL) sds_id = id
  
  attr_index = sffattr(sds_id, trim(attr_name))
  if (attr_index == FAIL) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(attr_name),sds_id
    error_status = sym%FAILURE
    return
  endif
  
  istatus = sfgainfo(sds_id, attr_index, name, sds_type, count)
  if (istatus /= 0) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(attr_name),sds_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (present(type)) then
    type = sds_type
    istatus = sfendacc(sds_id)
    error_status = sym%SUCCESS
    return
  endif
  
  if (sds_type /= DFNT_FLOAT64) then
    print "(a,'Error reading attribute (type mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(attr_name),sds_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (count /= 1) then
    print "(a,'Error reading attribute (count mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(attr_name),sds_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
    
  istatus = sfrnatt(sds_id, attr_index, buffer)
  if (istatus /= 0) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(attr_name),sds_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  istatus = sfendacc(sds_id)
  
  attr = buffer(1)
  
  return

END FUNCTION read_hdf_attribute_float64_scalar

!-------------------------------------------------------------------
! This routine is used to read int8 HDF SDS attributes.
!-------------------------------------------------------------------

FUNCTION read_hdf_attribute_int8_vector(id, sds_name, attr_name, attr, type) result(error_status)
  INTEGER(kind=int4), intent(in) :: id
  CHARACTER(len=*), intent(in) :: sds_name, attr_name
  INTEGER(kind=int1), dimension(:), allocatable, intent(inout) :: attr
  INTEGER(kind=int4), intent(out), optional :: type
  
  CHARACTER(len=1020) :: name
  INTEGER(kind=int4) :: istatus, attr_index, sds_type, count, sds_id, error_status, astatus
  INTEGER :: sffattr, sfrnatt, sfgainfo, sfn2index, sfselect, sfendacc
  
  error_status = sym%SUCCESS
  !attr = missing_value_int1  ! FIXME: what if it hasn't been allocated yet? 
  
  sds_id = sfselect(id, sfn2index(id,trim(sds_name)))
  if (sds_id == FAIL) sds_id = id
  
  attr_index = sffattr(sds_id, trim(attr_name))
  if (attr_index == FAIL) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(attr_name),sds_id
    error_status = sym%FAILURE
    return
  endif
  
  istatus = sfgainfo(sds_id, attr_index, name, sds_type, count)
  if (istatus /= 0) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(attr_name),sds_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (present(type)) then
    type = sds_type
    istatus = sfendacc(sds_id)
    error_status = sym%SUCCESS
    return
  endif
  
  if (sds_type /= DFNT_INT8 .and. &
      sds_type /= DFNT_UINT8 .and. &
      sds_type /= DFNT_CHAR8 .and. &
      sds_type /= DFNT_CHAR) then
    print "(a,'Error reading attribute (type mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(attr_name),sds_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (allocated(attr)) then
    if (size(attr,1) < count) then
      deallocate(attr,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 1d HDF attr buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif
  
  if ((.not. allocated(attr))) then
    allocate(attr(count),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 1d attr buffer.')",EXE_PROMPT
      stop
    endif
  endif
    
  istatus = sfrnatt(sds_id, attr_index, attr)
  if (istatus /= 0) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(attr_name),sds_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  istatus = sfendacc(sds_id)
  
  return

END FUNCTION read_hdf_attribute_int8_vector


!-------------------------------------------------------------------
! This routine is used to read int16 HDF SDS attributes.
!-------------------------------------------------------------------

FUNCTION read_hdf_attribute_int16_vector(id, sds_name, attr_name, attr, type) result(error_status)
  INTEGER(kind=int4), intent(in) :: id
  CHARACTER(len=*), intent(in) :: sds_name, attr_name
  INTEGER(kind=int2), dimension(:), allocatable, intent(inout) :: attr
  INTEGER(kind=int4), intent(out), optional :: type
  
  CHARACTER(len=1020) :: name
  INTEGER(kind=int4) :: istatus, attr_index, sds_type, count, sds_id, error_status, astatus
  INTEGER :: sffattr, sfrnatt, sfgainfo, sfn2index, sfselect, sfendacc
  
  error_status = sym%SUCCESS

  ! attr = missing_value_real4
  
  sds_id = sfselect(id, sfn2index(id,trim(sds_name)))
  if (sds_id == FAIL) sds_id = id
  
  attr_index = sffattr(sds_id, trim(attr_name))
  if (attr_index == FAIL) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(attr_name),sds_id
    error_status = sym%FAILURE
    return
  endif
  
  istatus = sfgainfo(sds_id, attr_index, name, sds_type, count)
  if (istatus /= 0) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(attr_name),sds_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (present(type)) then
    type = sds_type
    istatus = sfendacc(sds_id)
    error_status = sym%SUCCESS
    return
  endif
  
  if (sds_type /= DFNT_INT16 .AND. &
      sds_type /= DFNT_UINT16) then
    print "(a,'Error reading attribute (type mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(attr_name),sds_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (allocated(attr)) then
    if (size(attr,1) < count) then
      deallocate(attr,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 1d HDF attr buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif
  
  if ((.not. allocated(attr))) then
    allocate(attr(count),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 1d attr buffer.')",EXE_PROMPT
      stop
    endif
  endif

  attr = missing_value_int2
    
  istatus = sfrnatt(sds_id, attr_index, attr)
  if (istatus /= 0) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(attr_name),sds_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  istatus = sfendacc(sds_id)
  
  return

END FUNCTION read_hdf_attribute_int16_vector

!-------------------------------------------------------------------
! This routine is used to read int32 HDF SDS attributes.
!-------------------------------------------------------------------

FUNCTION read_hdf_attribute_int32_vector(id, sds_name, attr_name, attr, type) result(error_status)
  INTEGER(kind=int4), intent(in) :: id
  CHARACTER(len=*), intent(in) :: sds_name, attr_name
  INTEGER(kind=int4), dimension(:), allocatable, intent(inout) :: attr
  INTEGER(kind=int4), intent(out), optional :: type
  
  CHARACTER(len=1020) :: name
  INTEGER(kind=int4) :: istatus, attr_index, sds_type, count, sds_id, error_status, astatus
  INTEGER :: sffattr, sfrnatt, sfgainfo, sfn2index, sfselect, sfendacc
  
  error_status = sym%SUCCESS

  ! attr = missing_value_real4
  
  sds_id = sfselect(id, sfn2index(id,trim(sds_name)))
  if (sds_id == FAIL) sds_id = id
  
  attr_index = sffattr(sds_id, trim(attr_name))
  if (attr_index == FAIL) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(attr_name),sds_id
    error_status = sym%FAILURE
    return
  endif
  
  istatus = sfgainfo(sds_id, attr_index, name, sds_type, count)
  if (istatus /= 0) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(attr_name),sds_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (present(type)) then
    type = sds_type
    istatus = sfendacc(sds_id)
    error_status = sym%SUCCESS
    return
  endif
  
  if (sds_type /= DFNT_INT32 .AND. &
      sds_type /= DFNT_UINT32) then
    print "(a,'Error reading attribute (type mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(attr_name),sds_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (allocated(attr)) then
    if (size(attr,1) < count) then
      deallocate(attr,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 1d HDF attr buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif
  
  if ((.not. allocated(attr))) then
    allocate(attr(count),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 1d attr buffer.')",EXE_PROMPT
      stop
    endif
  endif

  attr = missing_value_int4
    
  istatus = sfrnatt(sds_id, attr_index, attr)
  if (istatus /= 0) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(attr_name),sds_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  istatus = sfendacc(sds_id)
  
  return

END FUNCTION read_hdf_attribute_int32_vector

!-------------------------------------------------------------------
! This routine is used to read float32 HDF SDS attributes.
!-------------------------------------------------------------------

FUNCTION read_hdf_attribute_float32_vector(id, sds_name, attr_name, attr, type) result(error_status)
  INTEGER(kind=int4), intent(in) :: id
  CHARACTER(len=*), intent(in) :: sds_name, attr_name
  REAL(kind=real4), dimension(:), allocatable, intent(inout) :: attr
  INTEGER(kind=int4), intent(out), optional :: type
  
  CHARACTER(len=1020) :: name
  INTEGER(kind=int4) :: istatus, attr_index, sds_type, count, sds_id, error_status, astatus
  INTEGER :: sffattr, sfrnatt, sfgainfo, sfn2index, sfselect, sfendacc
  
  error_status = sym%SUCCESS
  !attr = missing_value_real4
  
  sds_id = sfselect(id, sfn2index(id,trim(sds_name)))
  if (sds_id == FAIL) sds_id = id
  
  attr_index = sffattr(sds_id, trim(attr_name))
  if (attr_index == FAIL) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(attr_name),sds_id
    error_status = sym%FAILURE
    return
  endif
  
  istatus = sfgainfo(sds_id, attr_index, name, sds_type, count)
  if (istatus /= 0) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(attr_name),sds_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (present(type)) then
    type = sds_type
    istatus = sfendacc(sds_id)
    error_status = sym%SUCCESS
    return
  endif
  
  if (sds_type /= DFNT_FLOAT32) then
    print "(a,'Error reading attribute (type mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(attr_name),sds_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (allocated(attr)) then
    if (size(attr,1) < count) then
      deallocate(attr,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 1d HDF attr buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif
  
  if ((.not. allocated(attr))) then
    allocate(attr(count),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 1d attr buffer.')",EXE_PROMPT
      stop
    endif
  endif
    
  attr = missing_value_real4
  
  istatus = sfrnatt(sds_id, attr_index, attr)
  if (istatus /= 0) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(attr_name),sds_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  istatus = sfendacc(sds_id)
  
  return

END FUNCTION read_hdf_attribute_float32_vector

!-------------------------------------------------------------------
! This routine is used to read float32 HDF SDS attributes.
!-------------------------------------------------------------------

FUNCTION read_hdf_attribute_float64_vector(id, sds_name, attr_name, attr, type) result(error_status)
  INTEGER(kind=int4), intent(in) :: id
  CHARACTER(len=*), intent(in) :: sds_name, attr_name
  REAL(kind=real8), dimension(:), allocatable, intent(inout) :: attr
  INTEGER(kind=int4), intent(out), optional :: type
  
  CHARACTER(len=1020) :: name
  INTEGER(kind=int4) :: istatus, attr_index, sds_type, count, sds_id, error_status, astatus
  INTEGER :: sffattr, sfrnatt, sfgainfo, sfn2index, sfselect, sfendacc
  
  error_status = sym%SUCCESS
  !attr = missing_value_real8
  
  sds_id = sfselect(id, sfn2index(id,trim(sds_name)))
  if (sds_id == FAIL) sds_id = id
  
  attr_index = sffattr(sds_id, trim(attr_name))
  if (attr_index == FAIL) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(attr_name),sds_id
    error_status = sym%FAILURE
    return
  endif
  
  istatus = sfgainfo(sds_id, attr_index, name, sds_type, count)
  if (istatus /= 0) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(attr_name),sds_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (present(type)) then
    type = sds_type
    istatus = sfendacc(sds_id)
    error_status = sym%SUCCESS
    return
  endif
  
  if (sds_type /= DFNT_FLOAT64) then
    print "(a,'Error reading attribute (type mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(attr_name),sds_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (allocated(attr)) then
    if (size(attr,1) < count) then
      deallocate(attr,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 1d HDF attr buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif
  
  if ((.not. allocated(attr))) then
    allocate(attr(count),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 1d attr buffer.')",EXE_PROMPT
      stop
    endif
  endif
    
  attr = missing_value_real8

  istatus = sfrnatt(sds_id, attr_index, attr)
  if (istatus /= 0) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(attr_name),sds_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  istatus = sfendacc(sds_id)
  
  return

END FUNCTION read_hdf_attribute_float64_vector













!-------------------------------------------------------------------
! Subroutine to read 1D int8 hdf data.
!-------------------------------------------------------------------
  
FUNCTION read_hdf_sds_int8_1d(sd_id, sds_name, istart, istride, iedge, buffer, type) result(error_status)
  INTEGER(kind=int4), intent(in) :: sd_id
  CHARACTER(*), intent(in) :: sds_name
  INTEGER(kind=int4), dimension(1), intent(inout) :: istart, istride
  INTEGER(kind=int4), dimension(1), intent(inout) :: iedge
  INTEGER(kind=int1), dimension(:), allocatable, intent(inout) :: buffer
  INTEGER(kind=int4), intent(out), optional :: type
  INTEGER(kind=int4) :: astatus, istatus, sds_id, sds_rank, sds_type, sds_nattr, error_status
  INTEGER(kind=int4), dimension(1) :: sds_dims, max_iedge
    
  integer :: sfselect, sfn2index, sfrdata, sfendacc, sfginfo
  
  error_status = sym%SUCCESS
  
  sds_id = sfselect(sd_id, sfn2index(sd_id,trim(sds_name)))
  if (sds_id == FAIL) then
    print "(a,'Error selecting ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    error_status = sym%FAILURE
    return
  endif
    
  istatus = sfginfo(sds_id, trim(sds_name), sds_rank, sds_dims, sds_type, sds_nattr)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (present(type)) then
    type = sds_type
    istatus = sfendacc(sds_id)
    error_status = sym%SUCCESS
    return
  endif
  
  if (sds_type /= DFNT_INT8 .and. &
      sds_type /= DFNT_UINT8 .and. &
      sds_type /= DFNT_CHAR8 .and. &
      sds_type /= DFNT_CHAR) then
    print "(a,'Error reading (type mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (sds_rank /= 1) then
    print "(a,'Error reading (rank mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (istart(1) < 0) istart(1) = 0
  
  if (istride(1) < 1) istride(1) = 1
  
  if (iedge(1) < 0) iedge(1) = int(ceiling(real(sds_dims(1))/real(istride(1))))
  
  max_iedge(1) = int(ceiling(real(sds_dims(1) - istart(1))/real(istride(1))))
  iedge(1) = min(max_iedge(1),iedge(1))
    
  if (allocated(buffer)) then
    if (size(buffer,1) < iedge(1)) then
      deallocate(buffer,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 1d HDF buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif
  
  if ((.not. allocated(buffer))) then
    allocate(buffer(iedge(1)),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 1d buffer.')",EXE_PROMPT
      stop
    endif
  endif
    
  istatus = sfrdata(sds_id, istart, istride, iedge, buffer)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  istatus = sfendacc(sds_id)
  
  return
    
END FUNCTION read_hdf_sds_int8_1d

!-------------------------------------------------------------------
! Subroutine to read 1D int16 hdf data.
!-------------------------------------------------------------------
  
FUNCTION read_hdf_sds_int16_1d(sd_id, sds_name, istart, istride, iedge, buffer, type) result(error_status)
  INTEGER(kind=int4), intent(in) :: sd_id
  CHARACTER(*), intent(in) :: sds_name
  INTEGER(kind=int4), dimension(1), intent(inout) :: istart, istride
  INTEGER(kind=int4), dimension(1), intent(inout) :: iedge
  INTEGER(kind=int2), dimension(:), allocatable, intent(inout) :: buffer
  INTEGER(kind=int4), intent(out), optional :: type
  INTEGER(kind=int4) :: astatus, istatus, sds_id, sds_rank, sds_type, sds_nattr, error_status
  INTEGER(kind=int4), dimension(1) :: sds_dims, max_iedge
    
  integer :: sfselect, sfn2index, sfrdata, sfendacc, sfginfo
  
  error_status = sym%SUCCESS
  
  sds_id = sfselect(sd_id, sfn2index(sd_id,trim(sds_name)))
  if (sds_id == FAIL) then
    print "(a,'Error selecting ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    error_status = sym%FAILURE
    return
  endif
    
  istatus = sfginfo(sds_id, trim(sds_name), sds_rank, sds_dims, sds_type, sds_nattr)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (present(type)) then
    type = sds_type
    istatus = sfendacc(sds_id)
    error_status = sym%SUCCESS
    return
  endif
  
  if (sds_type /= DFNT_INT16) then
    print "(a,'Error reading (type mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (sds_rank /= 1) then
    print "(a,'Error reading (rank mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (istart(1) < 0) istart(1) = 0
  
  if (istride(1) < 1) istride(1) = 1
  
  if (iedge(1) < 0) iedge(1) = int(ceiling(real(sds_dims(1))/real(istride(1))))
  
  max_iedge(1) = int(ceiling(real(sds_dims(1) - istart(1))/real(istride(1))))
  iedge(1) = min(max_iedge(1),iedge(1))
    
  if (allocated(buffer)) then
    if (size(buffer,1) < iedge(1)) then
      deallocate(buffer,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 1d HDF buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif
  
  if ((.not. allocated(buffer))) then
    allocate(buffer(iedge(1)),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 1d buffer.')",EXE_PROMPT
      stop
    endif
  endif
    
  istatus = sfrdata(sds_id, istart, istride, iedge, buffer)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  istatus = sfendacc(sds_id)
  
  return
    
END FUNCTION read_hdf_sds_int16_1d

!-------------------------------------------------------------------
! Subroutine to read 1D int32 hdf data.
!-------------------------------------------------------------------
  
FUNCTION read_hdf_sds_int32_1d(sd_id, sds_name, istart, istride, iedge, buffer, type) result(error_status)
  INTEGER(kind=int4), intent(in) :: sd_id
  CHARACTER(*), intent(in) :: sds_name
  INTEGER(kind=int4), dimension(1), intent(inout) :: istart, istride
  INTEGER(kind=int4), dimension(1), intent(inout) :: iedge
  INTEGER(kind=int4), dimension(:), allocatable, intent(inout) :: buffer
  INTEGER(kind=int4), intent(out), optional :: type
  INTEGER(kind=int4) :: astatus, istatus, sds_id, sds_rank, sds_type, sds_nattr, error_status
  INTEGER(kind=int4), dimension(1) :: sds_dims, max_iedge
    
  integer :: sfselect, sfn2index, sfrdata, sfendacc, sfginfo
  
  error_status = sym%SUCCESS
  
  sds_id = sfselect(sd_id, sfn2index(sd_id,trim(sds_name)))
  if (sds_id == FAIL) then
    print "(a,'Error selecting ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    error_status = sym%FAILURE
    return
  endif
    
  istatus = sfginfo(sds_id, trim(sds_name), sds_rank, sds_dims, sds_type, sds_nattr)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (present(type)) then
    type = sds_type
    istatus = sfendacc(sds_id)
    error_status = sym%SUCCESS
    return
  endif
  
  if (sds_type /= DFNT_INT32) then
    print "(a,'Error reading (type mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (sds_rank /= 1) then
    print "(a,'Error reading (rank mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (istart(1) < 0) istart(1) = 0
  
  if (istride(1) < 1) istride(1) = 1
  
  if (iedge(1) < 0) iedge(1) = int(ceiling(real(sds_dims(1))/real(istride(1))))
  
  max_iedge(1) = int(ceiling(real(sds_dims(1) - istart(1))/real(istride(1))))
  iedge(1) = min(max_iedge(1),iedge(1))
    
  if (allocated(buffer)) then
    if (size(buffer,1) < iedge(1)) then
      deallocate(buffer,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 1d HDF buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif
  
  if ((.not. allocated(buffer))) then
    allocate(buffer(iedge(1)),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 1d buffer.')",EXE_PROMPT
      stop
    endif
  endif
    
  istatus = sfrdata(sds_id, istart, istride, iedge, buffer)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  istatus = sfendacc(sds_id)
  
  return
    
END FUNCTION read_hdf_sds_int32_1d

!-------------------------------------------------------------------
! Subroutine to read 1D float32 hdf data.
!-------------------------------------------------------------------
  
FUNCTION read_hdf_sds_float32_1d(sd_id, sds_name, istart, istride, iedge, buffer, type) result(error_status)
  INTEGER(kind=int4), intent(in) :: sd_id
  CHARACTER(*), intent(in) :: sds_name
  INTEGER(kind=int4), dimension(1), intent(inout) :: istart, istride
  INTEGER(kind=int4), dimension(1), intent(inout) :: iedge
  REAL(kind=real4), dimension(:), allocatable, intent(inout) :: buffer
  INTEGER(kind=int4), intent(out), optional :: type
  INTEGER(kind=int4) :: astatus, istatus, sds_id, sds_rank, sds_type, sds_nattr, error_status
  INTEGER(kind=int4), dimension(1) :: sds_dims, max_iedge
    
  integer :: sfselect, sfn2index, sfrdata, sfendacc, sfginfo
  
  error_status = sym%SUCCESS
  
  sds_id = sfselect(sd_id, sfn2index(sd_id,trim(sds_name)))
  if (sds_id == FAIL) then
    print "(a,'Error selecting ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    error_status = sym%FAILURE
    return
  endif
    
  istatus = sfginfo(sds_id, trim(sds_name), sds_rank, sds_dims, sds_type, sds_nattr)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (present(type)) then
    type = sds_type
    istatus = sfendacc(sds_id)
    error_status = sym%SUCCESS
    return
  endif
  
  if (sds_type /= DFNT_FLOAT32) then
    print "(a,'Error reading (type mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (sds_rank /= 1) then
    print "(a,'Error reading (rank mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (istart(1) < 0) istart(1) = 0
  
  if (istride(1) < 1) istride(1) = 1
  
  if (iedge(1) < 0) iedge(1) = int(ceiling(real(sds_dims(1))/real(istride(1))))
  
  max_iedge(1) = int(ceiling(real(sds_dims(1) - istart(1))/real(istride(1))))
  iedge(1) = min(max_iedge(1),iedge(1))
  
  if (allocated(buffer)) then
    if (size(buffer,1) < iedge(1)) then
      deallocate(buffer,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 1d HDF buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif
  
  if ((.not. allocated(buffer))) then
    allocate(buffer(iedge(1)),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 1d buffer.')",EXE_PROMPT
      stop
    endif
  endif
    
  istatus = sfrdata(sds_id, istart, istride, iedge, buffer)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  istatus = sfendacc(sds_id)
  
  return
    
END FUNCTION read_hdf_sds_float32_1d

!-------------------------------------------------------------------
! Subroutine to read 1D float64 hdf data.
!-------------------------------------------------------------------
  
FUNCTION read_hdf_sds_float64_1d(sd_id, sds_name, istart, istride, iedge, buffer, type) result(error_status)
  INTEGER(kind=int4), intent(in) :: sd_id
  CHARACTER(*), intent(in) :: sds_name
  INTEGER(kind=int4), dimension(1), intent(inout) :: istart, istride
  INTEGER(kind=int4), dimension(1), intent(inout) :: iedge
  REAL(kind=real8), dimension(:), allocatable, intent(inout) :: buffer
  INTEGER(kind=int4), intent(out), optional :: type
  INTEGER(kind=int4) :: astatus, istatus, sds_id, sds_rank, sds_type, sds_nattr, error_status
  INTEGER(kind=int4), dimension(1) :: sds_dims, max_iedge
    
  integer :: sfselect, sfn2index, sfrdata, sfendacc, sfginfo
  
  error_status = sym%SUCCESS
  
  sds_id = sfselect(sd_id, sfn2index(sd_id,trim(sds_name)))
  if (sds_id == FAIL) then
    print "(a,'Error selecting ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    error_status = sym%FAILURE
    return
  endif
    
  istatus = sfginfo(sds_id, trim(sds_name), sds_rank, sds_dims, sds_type, sds_nattr)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (present(type)) then
    type = sds_type
    istatus = sfendacc(sds_id)
    error_status = sym%SUCCESS
    return
  endif
  
  if (sds_type /= DFNT_FLOAT64) then
    print "(a,'Error reading (type mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (sds_rank /= 1) then
    print "(a,'Error reading (rank mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (istart(1) < 0) istart(1) = 0
  
  if (istride(1) < 1) istride(1) = 1
  
  if (iedge(1) < 0) iedge(1) = int(ceiling(real(sds_dims(1))/real(istride(1))))
  
  max_iedge(1) = int(ceiling(real(sds_dims(1) - istart(1))/real(istride(1))))
  iedge(1) = min(max_iedge(1),iedge(1))
    
  if (allocated(buffer)) then
    if (size(buffer,1) < iedge(1)) then
      deallocate(buffer,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 1d HDF buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif
  
  if ((.not. allocated(buffer))) then
    allocate(buffer(iedge(1)),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 1d buffer.')",EXE_PROMPT
      stop
    endif
  endif
    
  istatus = sfrdata(sds_id, istart, istride, iedge, buffer)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  istatus = sfendacc(sds_id)
  
  return
    
END FUNCTION read_hdf_sds_float64_1d
  
!-------------------------------------------------------------------
! Subroutine to read 2D int8 hdf data.
!-------------------------------------------------------------------
  
FUNCTION read_hdf_sds_int8_2d(sd_id, sds_name, istart, istride, iedge, buffer, type) result(error_status)
  INTEGER(kind=int4), intent(in) :: sd_id
  CHARACTER(*), intent(in) :: sds_name
  INTEGER(kind=int4), dimension(2), intent(inout) :: istart, istride
  INTEGER(kind=int4), dimension(2), intent(inout) :: iedge
  INTEGER(kind=int1), dimension(:,:), allocatable, intent(inout) :: buffer
  INTEGER(kind=int4), intent(out), optional :: type
  INTEGER(kind=int4) :: astatus, istatus, sds_id, sds_rank, sds_type, sds_nattr, error_status
  INTEGER(kind=int4), dimension(2) :: sds_dims, max_iedge
    
  integer :: sfselect, sfn2index, sfrdata, sfendacc, sfginfo
  
  error_status = sym%SUCCESS
  
  sds_id = sfselect(sd_id, sfn2index(sd_id,trim(sds_name)))
  if (sds_id == FAIL) then
    print "(a,'Error selecting ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    error_status = sym%FAILURE
    return
  endif
    
  istatus = sfginfo(sds_id, trim(sds_name), sds_rank, sds_dims, sds_type, sds_nattr)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (present(type)) then
    type = sds_type
    istatus = sfendacc(sds_id)
    error_status = sym%SUCCESS
    return
  endif
  
  if (sds_type /= DFNT_INT8 .and. &
      sds_type /= DFNT_UINT8 .and. &
      sds_type /= DFNT_CHAR8 .and. &
      sds_type /= DFNT_CHAR) then
    print "(a,'Error reading (type mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (sds_rank /= 2) then
    print "(a,'Error reading (rank mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (istart(1) < 0) istart(1) = 0
  if (istart(2) < 0) istart(2) = 0
  
  if (istride(1) < 1) istride(1) = 1
  if (istride(2) < 1) istride(2) = 1
  
  if (iedge(1) < 0) iedge(1) = int(ceiling(real(sds_dims(1))/real(istride(1))))
  if (iedge(2) < 0) iedge(2) = int(ceiling(real(sds_dims(2))/real(istride(2))))
  
  max_iedge(1) = int(ceiling(real(sds_dims(1) - istart(1))/real(istride(1))))
  max_iedge(2) = int(ceiling(real(sds_dims(2) - istart(2))/real(istride(2))))
  iedge(1) = min(max_iedge(1),iedge(1))
  iedge(2) = min(max_iedge(2),iedge(2))
    
  if (allocated(buffer)) then
    if (size(buffer,1) < iedge(1) .or. &
        size(buffer,2) < iedge(2)) then
      deallocate(buffer,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 2d HDF buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif
  
  if ((.not. allocated(buffer))) then
    allocate(buffer(iedge(1),iedge(2)),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 2d buffer.')",EXE_PROMPT
      stop
    endif
  endif
    
  istatus = sfrdata(sds_id, istart, istride, iedge, buffer)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  istatus = sfendacc(sds_id)
  
  return
    
END FUNCTION read_hdf_sds_int8_2d

!-------------------------------------------------------------------
! Subroutine to read 2D int16 hdf data.
!-------------------------------------------------------------------
  
FUNCTION read_hdf_sds_int16_2d(sd_id, sds_name, istart, istride, iedge, buffer, type) result(error_status)
  INTEGER(kind=int4), intent(in) :: sd_id
  CHARACTER(*), intent(in) :: sds_name
  INTEGER(kind=int4), dimension(2), intent(inout) :: istart, istride
  INTEGER(kind=int4), dimension(2), intent(inout) :: iedge
  INTEGER(kind=int2), dimension(:,:), allocatable, intent(inout) :: buffer
  INTEGER(kind=int4), intent(out), optional :: type
  INTEGER(kind=int4) :: astatus, istatus, sds_id, sds_rank, sds_type, sds_nattr, error_status
  INTEGER(kind=int4), dimension(2) :: sds_dims, max_iedge
    
  integer :: sfselect, sfn2index, sfrdata, sfendacc, sfginfo
  
  error_status = sym%SUCCESS
  
  sds_id = sfselect(sd_id, sfn2index(sd_id,trim(sds_name)))
  if (sds_id == FAIL) then
    print "(a,'Error selecting ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    error_status = sym%FAILURE
    return
  endif
    
  istatus = sfginfo(sds_id, trim(sds_name), sds_rank, sds_dims, sds_type, sds_nattr)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (present(type)) then
    type = sds_type
    istatus = sfendacc(sds_id)
    error_status = sym%SUCCESS
    return
  endif
  
  if (sds_type /= DFNT_INT16) then
    print "(a,'Error reading (type mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (sds_rank /= 2) then
    print "(a,'Error reading (rank mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (istart(1) < 0) istart(1) = 0
  if (istart(2) < 0) istart(2) = 0
  
  if (istride(1) < 1) istride(1) = 1
  if (istride(2) < 1) istride(2) = 1
  
  if (iedge(1) < 0) iedge(1) = int(ceiling(real(sds_dims(1))/real(istride(1))))
  if (iedge(2) < 0) iedge(2) = int(ceiling(real(sds_dims(2))/real(istride(2))))
  
  max_iedge(1) = int(ceiling(real(sds_dims(1) - istart(1))/real(istride(1))))
  max_iedge(2) = int(ceiling(real(sds_dims(2) - istart(2))/real(istride(2))))
  iedge(1) = min(max_iedge(1),iedge(1))
  iedge(2) = min(max_iedge(2),iedge(2))
    
  if (allocated(buffer)) then
    if (size(buffer,1) < iedge(1) .or. &
        size(buffer,2) < iedge(2)) then
      deallocate(buffer,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 2d HDF buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif
  
  if ((.not. allocated(buffer))) then
    allocate(buffer(iedge(1),iedge(2)),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 2d buffer.')",EXE_PROMPT
      stop
    endif
  endif
    
  istatus = sfrdata(sds_id, istart, istride, iedge, buffer)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  istatus = sfendacc(sds_id)
  
  return
    
END FUNCTION read_hdf_sds_int16_2d

!-------------------------------------------------------------------
! Subroutine to read 2D int32 hdf data.
!-------------------------------------------------------------------
  
FUNCTION read_hdf_sds_int32_2d(sd_id, sds_name, istart, istride, iedge, buffer, type) result(error_status)
  INTEGER(kind=int4), intent(in) :: sd_id
  CHARACTER(*), intent(in) :: sds_name
  INTEGER(kind=int4), dimension(2), intent(inout) :: istart, istride
  INTEGER(kind=int4), dimension(2), intent(inout) :: iedge
  INTEGER(kind=int4), dimension(:,:), allocatable, intent(inout) :: buffer
  INTEGER(kind=int4), intent(out), optional :: type
  INTEGER(kind=int4) :: astatus, istatus, sds_id, sds_rank, sds_type, sds_nattr, error_status
  INTEGER(kind=int4), dimension(2) :: sds_dims, max_iedge
    
  integer :: sfselect, sfn2index, sfrdata, sfendacc, sfginfo
  
  error_status = sym%SUCCESS
  
  sds_id = sfselect(sd_id, sfn2index(sd_id,trim(sds_name)))
  if (sds_id == FAIL) then
    print "(a,'Error selecting ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    error_status = sym%FAILURE
    return
  endif
    
  istatus = sfginfo(sds_id, trim(sds_name), sds_rank, sds_dims, sds_type, sds_nattr)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (present(type)) then
    type = sds_type
    istatus = sfendacc(sds_id)
    error_status = sym%SUCCESS
    return
  endif
  
  if (sds_type /= DFNT_INT32) then
    print "(a,'Error reading (type mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (sds_rank /= 2) then
    print "(a,'Error reading (rank mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (istart(1) < 0) istart(1) = 0
  if (istart(2) < 0) istart(2) = 0
  
  if (istride(1) < 1) istride(1) = 1
  if (istride(2) < 1) istride(2) = 1
  
  if (iedge(1) < 0) iedge(1) = int(ceiling(real(sds_dims(1))/real(istride(1))))
  if (iedge(2) < 0) iedge(2) = int(ceiling(real(sds_dims(2))/real(istride(2))))
  
  max_iedge(1) = int(ceiling(real(sds_dims(1) - istart(1))/real(istride(1))))
  max_iedge(2) = int(ceiling(real(sds_dims(2) - istart(2))/real(istride(2))))
  iedge(1) = min(max_iedge(1),iedge(1))
  iedge(2) = min(max_iedge(2),iedge(2))
    
  if (allocated(buffer)) then
    if (size(buffer,1) < iedge(1) .or. &
        size(buffer,2) < iedge(2)) then
      deallocate(buffer,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 2d HDF buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif
  
  if ((.not. allocated(buffer))) then
    allocate(buffer(iedge(1),iedge(2)),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 2d buffer.')",EXE_PROMPT
      stop
    endif
  endif
    
  istatus = sfrdata(sds_id, istart, istride, iedge, buffer)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  istatus = sfendacc(sds_id)
  
  return
    
END FUNCTION read_hdf_sds_int32_2d

!-------------------------------------------------------------------
! Subroutine to read 2D float32 hdf data.
!-------------------------------------------------------------------
  
FUNCTION read_hdf_sds_float32_2d(sd_id, sds_name, istart, istride, iedge, buffer, type) result(error_status)
  INTEGER(kind=int4), intent(in) :: sd_id
  CHARACTER(*), intent(in) :: sds_name
  INTEGER(kind=int4), dimension(2), intent(inout) :: istart, istride
  INTEGER(kind=int4), dimension(2), intent(inout) :: iedge
  REAL(kind=real4), dimension(:,:), allocatable, intent(inout) :: buffer
  INTEGER(kind=int4), intent(out), optional :: type
  INTEGER(kind=int4) :: astatus, istatus, sds_id, sds_rank, sds_type, sds_nattr, error_status
  INTEGER(kind=int4), dimension(2) :: sds_dims, max_iedge
    
  integer :: sfselect, sfn2index, sfrdata, sfendacc, sfginfo
  
  error_status = sym%SUCCESS
  
  sds_id = sfselect(sd_id, sfn2index(sd_id,trim(sds_name)))
  if (sds_id == FAIL) then
    print "(a,'Error selecting ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    error_status = sym%FAILURE
    return
  endif
    
  istatus = sfginfo(sds_id, trim(sds_name), sds_rank, sds_dims, sds_type, sds_nattr)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (present(type)) then
    type = sds_type
    istatus = sfendacc(sds_id)
    error_status = sym%SUCCESS
    return
  endif
  
  if (sds_type /= DFNT_FLOAT32) then
    print "(a,'Error reading (type mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (sds_rank /= 2) then
    print "(a,'Error reading (rank mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (istart(1) < 0) istart(1) = 0
  if (istart(2) < 0) istart(2) = 0
  
  if (istride(1) < 1) istride(1) = 1
  if (istride(2) < 1) istride(2) = 1
  
  if (iedge(1) < 0) iedge(1) = int(ceiling(real(sds_dims(1))/real(istride(1))))
  if (iedge(2) < 0) iedge(2) = int(ceiling(real(sds_dims(2))/real(istride(2))))
  
  max_iedge(1) = int(ceiling(real(sds_dims(1) - istart(1))/real(istride(1))))
  max_iedge(2) = int(ceiling(real(sds_dims(2) - istart(2))/real(istride(2))))
  iedge(1) = min(max_iedge(1),iedge(1))
  iedge(2) = min(max_iedge(2),iedge(2))
  
  if (allocated(buffer)) then
    if (size(buffer,1) < iedge(1) .or. &
        size(buffer,2) < iedge(2)) then
      deallocate(buffer,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 2d HDF buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif
  
  if ((.not. allocated(buffer))) then
    allocate(buffer(iedge(1),iedge(2)),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 2d buffer.')",EXE_PROMPT
      stop
    endif
  endif
    
  istatus = sfrdata(sds_id, istart, istride, iedge, buffer)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  istatus = sfendacc(sds_id)
  
  return
    
END FUNCTION read_hdf_sds_float32_2d

!-------------------------------------------------------------------
! Subroutine to read 2D float64 hdf data.
!-------------------------------------------------------------------
  
FUNCTION read_hdf_sds_float64_2d(sd_id, sds_name, istart, istride, iedge, buffer, type) result(error_status)
  INTEGER(kind=int4), intent(in) :: sd_id
  CHARACTER(*), intent(in) :: sds_name
  INTEGER(kind=int4), dimension(2), intent(inout) :: istart, istride
  INTEGER(kind=int4), dimension(2), intent(inout) :: iedge
  REAL(kind=real8), dimension(:,:), allocatable, intent(inout) :: buffer
  INTEGER(kind=int4), intent(out), optional :: type
  INTEGER(kind=int4) :: astatus, istatus, sds_id, sds_rank, sds_type, sds_nattr, error_status
  INTEGER(kind=int4), dimension(2) :: sds_dims, max_iedge
    
  integer :: sfselect, sfn2index, sfrdata, sfendacc, sfginfo
  
  error_status = sym%SUCCESS
  
  sds_id = sfselect(sd_id, sfn2index(sd_id,trim(sds_name)))
  if (sds_id == FAIL) then
    print "(a,'Error selecting ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    error_status = sym%FAILURE
    return
  endif
    
  istatus = sfginfo(sds_id, trim(sds_name), sds_rank, sds_dims, sds_type, sds_nattr)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (present(type)) then
    type = sds_type
    istatus = sfendacc(sds_id)
    error_status = sym%SUCCESS
    return
  endif
  
  if (sds_type /= DFNT_FLOAT64) then
    print "(a,'Error reading (type mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (sds_rank /= 2) then
    print "(a,'Error reading (rank mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (istart(1) < 0) istart(1) = 0
  if (istart(2) < 0) istart(2) = 0
  
  if (istride(1) < 1) istride(1) = 1
  if (istride(2) < 1) istride(2) = 1
  
  if (iedge(1) < 0) iedge(1) = int(ceiling(real(sds_dims(1))/real(istride(1))))
  if (iedge(2) < 0) iedge(2) = int(ceiling(real(sds_dims(2))/real(istride(2))))
  
  max_iedge(1) = int(ceiling(real(sds_dims(1) - istart(1))/real(istride(1))))
  max_iedge(2) = int(ceiling(real(sds_dims(2) - istart(2))/real(istride(2))))
  iedge(1) = min(max_iedge(1),iedge(1))
  iedge(2) = min(max_iedge(2),iedge(2))
    
  if (allocated(buffer)) then
    if (size(buffer,1) < iedge(1) .or. &
        size(buffer,2) < iedge(2)) then
      deallocate(buffer,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 2d HDF buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif
  
  if ((.not. allocated(buffer))) then
    allocate(buffer(iedge(1),iedge(2)),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 2d buffer.')",EXE_PROMPT
      stop
    endif
  endif
    
  istatus = sfrdata(sds_id, istart, istride, iedge, buffer)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  istatus = sfendacc(sds_id)
  
  return
    
END FUNCTION read_hdf_sds_float64_2d


!-------------------------------------------------------------------
! Subroutine to read 3D int8 hdf data.
!-------------------------------------------------------------------
  
FUNCTION read_hdf_sds_int8_3d(sd_id, sds_name, istart, istride, iedge, buffer, type) result(error_status)
  INTEGER(kind=int4), intent(in) :: sd_id
  CHARACTER(*), intent(in) :: sds_name
  INTEGER(kind=int4), dimension(3), intent(inout) :: istart, istride
  INTEGER(kind=int4), dimension(3), intent(inout) :: iedge
  INTEGER(kind=int1), dimension(:,:,:), allocatable, intent(inout) :: buffer
  INTEGER(kind=int4), intent(out), optional :: type
  INTEGER(kind=int4) :: astatus, istatus, sds_id, sds_rank, sds_type, sds_nattr, error_status, idim
  INTEGER(kind=int4), dimension(3) :: sds_dims, max_iedge
    
  integer :: sfselect, sfn2index, sfrdata, sfendacc, sfginfo
  
  error_status = sym%SUCCESS
  
  sds_id = sfselect(sd_id, sfn2index(sd_id,trim(sds_name)))
  if (sds_id == FAIL) then
    print "(a,'Error selecting ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    error_status = sym%FAILURE
    return
  endif
    
  istatus = sfginfo(sds_id, trim(sds_name), sds_rank, sds_dims, sds_type, sds_nattr)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (present(type)) then
    type = sds_type
    istatus = sfendacc(sds_id)
    error_status = sym%SUCCESS
    return
  endif
  
  if (sds_type /= DFNT_INT8 .and. &
      sds_type /= DFNT_UINT8 .and. &
      sds_type /= DFNT_CHAR8 .and. &
      sds_type /= DFNT_CHAR) then
    print "(a,'Error reading (type mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (sds_rank /= 3) then
    print "(a,'Error reading (rank mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  do idim=1, sds_rank
  
    if (istart(idim) < 0) istart(idim) = 0
    if (istride(idim) < 1) istride(idim) = 1
  
    if (iedge(idim) < 0) iedge(idim) = int(ceiling(real(sds_dims(idim))/real(istride(idim))))
    max_iedge(idim) = int(ceiling(real(sds_dims(idim) - istart(idim))/real(istride(idim))))
    iedge(idim) = min(max_iedge(idim),iedge(idim))
    
  end do
    
  if (allocated(buffer)) then
    if (size(buffer,1) < iedge(1) .or. &
        size(buffer,2) < iedge(2) .or. &
size(buffer,3) < iedge(3)) then
      deallocate(buffer,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 3d HDF buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif
  
  if ((.not. allocated(buffer))) then
    allocate(buffer(iedge(1),iedge(2),iedge(3)),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 3d buffer.')",EXE_PROMPT
      stop
    endif
  endif
    
  istatus = sfrdata(sds_id, istart, istride, iedge, buffer)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  istatus = sfendacc(sds_id)
  
  return
    
END FUNCTION read_hdf_sds_int8_3d

!-------------------------------------------------------------------
! Subroutine to read 3D int16 hdf data.
!-------------------------------------------------------------------
  
FUNCTION read_hdf_sds_int16_3d(sd_id, sds_name, istart, istride, iedge, buffer, type) result(error_status)
  INTEGER(kind=int4), intent(in) :: sd_id
  CHARACTER(*), intent(in) :: sds_name
  INTEGER(kind=int4), dimension(3), intent(inout) :: istart, istride
  INTEGER(kind=int4), dimension(3), intent(inout) :: iedge
  INTEGER(kind=int2), dimension(:,:,:), allocatable, intent(inout) :: buffer
  INTEGER(kind=int4), intent(out), optional :: type
  INTEGER(kind=int4) :: astatus, istatus, sds_id, sds_rank, sds_type, sds_nattr, error_status, idim
  INTEGER(kind=int4), dimension(3) :: sds_dims, max_iedge
    
  integer :: sfselect, sfn2index, sfrdata, sfendacc, sfginfo
  
  error_status = sym%SUCCESS
  
  sds_id = sfselect(sd_id, sfn2index(sd_id,trim(sds_name)))
  if (sds_id == FAIL) then
    print "(a,'Error selecting ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    error_status = sym%FAILURE
    return
  endif
    
  istatus = sfginfo(sds_id, trim(sds_name), sds_rank, sds_dims, sds_type, sds_nattr)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (present(type)) then
    type = sds_type
    istatus = sfendacc(sds_id)
    error_status = sym%SUCCESS
    return
  endif
  
  if (sds_type /= DFNT_INT16) then
    print "(a,'Error reading (type mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (sds_rank /= 3) then
    print "(a,'Error reading (rank mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  do idim=1, sds_rank
  
    if (istart(idim) < 0) istart(idim) = 0
    if (istride(idim) < 1) istride(idim) = 1
  
    if (iedge(idim) < 0) iedge(idim) = int(ceiling(real(sds_dims(idim))/real(istride(idim))))
    max_iedge(idim) = int(ceiling(real(sds_dims(idim) - istart(idim))/real(istride(idim))))
    iedge(idim) = min(max_iedge(idim),iedge(idim))
    
  end do
    
  if (allocated(buffer)) then
    if (size(buffer,1) < iedge(1) .or. &
        size(buffer,2) < iedge(2) .or. &
       size(buffer,3) < iedge(3)) then
      deallocate(buffer,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 3d HDF buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif
  
  if ((.not. allocated(buffer))) then
    allocate(buffer(iedge(1),iedge(2),iedge(3)),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 3d buffer.')",EXE_PROMPT
      stop
    endif
  endif
    
  istatus = sfrdata(sds_id, istart, istride, iedge, buffer)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  istatus = sfendacc(sds_id)
  
  return
    
END FUNCTION read_hdf_sds_int16_3d

!-------------------------------------------------------------------
! Subroutine to read 3D int32 hdf data.
!-------------------------------------------------------------------
  
FUNCTION read_hdf_sds_int32_3d(sd_id, sds_name, istart, istride, iedge, buffer, type) result(error_status)
  INTEGER(kind=int4), intent(in) :: sd_id
  CHARACTER(*), intent(in) :: sds_name
  INTEGER(kind=int4), dimension(3), intent(inout) :: istart, istride
  INTEGER(kind=int4), dimension(3), intent(inout) :: iedge
  INTEGER(kind=int4), dimension(:,:,:), allocatable, intent(inout) :: buffer
  INTEGER(kind=int4), intent(out), optional :: type
  INTEGER(kind=int4) :: astatus, istatus, sds_id, sds_rank, sds_type, sds_nattr, error_status, idim
  INTEGER(kind=int4), dimension(3) :: sds_dims, max_iedge
    
  integer :: sfselect, sfn2index, sfrdata, sfendacc, sfginfo
  
  error_status = sym%SUCCESS
  
  sds_id = sfselect(sd_id, sfn2index(sd_id,trim(sds_name)))
  if (sds_id == FAIL) then
    print "(a,'Error selecting ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    error_status = sym%FAILURE
    return
  endif
    
  istatus = sfginfo(sds_id, trim(sds_name), sds_rank, sds_dims, sds_type, sds_nattr)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (present(type)) then
    type = sds_type
    istatus = sfendacc(sds_id)
    error_status = sym%SUCCESS
    return
  endif
  
  if (sds_type /= DFNT_INT32) then
    print "(a,'Error reading (type mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (sds_rank /= 3) then
    print "(a,'Error reading (rank mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  do idim=1, sds_rank
  
    if (istart(idim) < 0) istart(idim) = 0
    if (istride(idim) < 1) istride(idim) = 1
  
    if (iedge(idim) < 0) iedge(idim) = int(ceiling(real(sds_dims(idim))/real(istride(idim))))
    max_iedge(idim) = int(ceiling(real(sds_dims(idim) - istart(idim))/real(istride(idim))))
    iedge(idim) = min(max_iedge(idim),iedge(idim))
    
  end do
    
  if (allocated(buffer)) then
    if (size(buffer,1) < iedge(1) .or. &
        size(buffer,2) < iedge(2) .or. &
size(buffer,3) < iedge(3)) then
      deallocate(buffer,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 3d HDF buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif
  
  if ((.not. allocated(buffer))) then
    allocate(buffer(iedge(1),iedge(2),iedge(3)),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 3d buffer.')",EXE_PROMPT
      stop
    endif
  endif
    
  istatus = sfrdata(sds_id, istart, istride, iedge, buffer)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  istatus = sfendacc(sds_id)
  
  return
    
END FUNCTION read_hdf_sds_int32_3d

!-------------------------------------------------------------------
! Subroutine to read 3D float32 hdf data.
!-------------------------------------------------------------------
  
FUNCTION read_hdf_sds_float32_3d(sd_id, sds_name, istart, istride, iedge, buffer, type) result(error_status)
  INTEGER(kind=int4), intent(in) :: sd_id
  CHARACTER(*), intent(in) :: sds_name
  INTEGER(kind=int4), dimension(3), intent(inout) :: istart, istride
  INTEGER(kind=int4), dimension(3), intent(inout) :: iedge
  REAL(kind=real4), dimension(:,:,:), allocatable, intent(inout) :: buffer
  INTEGER(kind=int4), intent(out), optional :: type
  INTEGER(kind=int4) :: astatus, istatus, sds_id, sds_rank, sds_type, sds_nattr, error_status, idim
  INTEGER(kind=int4), dimension(3) :: sds_dims, max_iedge
    
  integer :: sfselect, sfn2index, sfrdata, sfendacc, sfginfo
  
  error_status = sym%SUCCESS
  
  sds_id = sfselect(sd_id, sfn2index(sd_id,trim(sds_name)))
  if (sds_id == FAIL) then
    print "(a,'Error selecting ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    error_status = sym%FAILURE
    return
  endif
    
  istatus = sfginfo(sds_id, trim(sds_name), sds_rank, sds_dims, sds_type, sds_nattr)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (present(type)) then
    type = sds_type
    istatus = sfendacc(sds_id)
    error_status = sym%SUCCESS
    return
  endif
  
  if (sds_type /= DFNT_FLOAT32) then
    print "(a,'Error reading (type mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (sds_rank /= 3) then
    print "(a,'Error reading (rank mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  do idim=1, sds_rank
  
    if (istart(idim) < 0) istart(idim) = 0
    if (istride(idim) < 1) istride(idim) = 1
  
    if (iedge(idim) < 0) iedge(idim) = int(ceiling(real(sds_dims(idim))/real(istride(idim))))
    max_iedge(idim) = int(ceiling(real(sds_dims(idim) - istart(idim))/real(istride(idim))))
    iedge(idim) = min(max_iedge(idim),iedge(idim))
    
  end do
  
  if (allocated(buffer)) then
    if (size(buffer,1) < iedge(1) .or. &
        size(buffer,2) < iedge(2) .or. &
size(buffer,3) < iedge(3)) then
      deallocate(buffer,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 3d HDF buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif
  
  if ((.not. allocated(buffer))) then
    allocate(buffer(iedge(1),iedge(2),iedge(3)),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 3d buffer.')",EXE_PROMPT
      stop
    endif
  endif
    
  istatus = sfrdata(sds_id, istart, istride, iedge, buffer)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  istatus = sfendacc(sds_id)
  
  return
    
END FUNCTION read_hdf_sds_float32_3d

!-------------------------------------------------------------------
! Subroutine to read 3D float64 hdf data.
!-------------------------------------------------------------------
  
FUNCTION read_hdf_sds_float64_3d(sd_id, sds_name, istart, istride, iedge, buffer, type) result(error_status)
  INTEGER(kind=int4), intent(in) :: sd_id
  CHARACTER(*), intent(in) :: sds_name
  INTEGER(kind=int4), dimension(3), intent(inout) :: istart, istride
  INTEGER(kind=int4), dimension(3), intent(inout) :: iedge
  REAL(kind=real8), dimension(:,:,:), allocatable, intent(inout) :: buffer
  INTEGER(kind=int4), intent(out), optional :: type
  INTEGER(kind=int4) :: astatus, istatus, sds_id, sds_rank, sds_type, sds_nattr, error_status, idim
  INTEGER(kind=int4), dimension(3) :: sds_dims, max_iedge
    
  integer :: sfselect, sfn2index, sfrdata, sfendacc, sfginfo
  
  error_status = sym%SUCCESS
  
  sds_id = sfselect(sd_id, sfn2index(sd_id,trim(sds_name)))
  if (sds_id == FAIL) then
    print "(a,'Error selecting ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    error_status = sym%FAILURE
    return
  endif
    
  istatus = sfginfo(sds_id, trim(sds_name), sds_rank, sds_dims, sds_type, sds_nattr)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (present(type)) then
    type = sds_type
    istatus = sfendacc(sds_id)
    error_status = sym%SUCCESS
    return
  endif
  
  if (sds_type /= DFNT_FLOAT64) then
    print "(a,'Error reading (type mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (sds_rank /= 3) then
    print "(a,'Error reading (rank mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  do idim=1, sds_rank
  
    if (istart(idim) < 0) istart(idim) = 0
    if (istride(idim) < 1) istride(idim) = 1
  
    if (iedge(idim) < 0) iedge(idim) = int(ceiling(real(sds_dims(idim))/real(istride(idim))))
    max_iedge(idim) = int(ceiling(real(sds_dims(idim) - istart(idim))/real(istride(idim))))
    iedge(idim) = min(max_iedge(idim),iedge(idim))
    
  end do
    
  if (allocated(buffer)) then
    if (size(buffer,1) < iedge(1) .or. &
        size(buffer,2) < iedge(2) .or. &
size(buffer,3) < iedge(3)) then
      deallocate(buffer,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 3d HDF buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif
  
  if ((.not. allocated(buffer))) then
    allocate(buffer(iedge(1),iedge(2),iedge(3)),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 3d buffer.')",EXE_PROMPT
      stop
    endif
  endif
    
  istatus = sfrdata(sds_id, istart, istride, iedge, buffer)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  istatus = sfendacc(sds_id)
  
  return
    
END FUNCTION read_hdf_sds_float64_3d

!-------------------------------------------------------------------
! Subroutine to read 4D int8 hdf data.
!-------------------------------------------------------------------
  
FUNCTION read_hdf_sds_int8_4d(sd_id, sds_name, istart, istride, iedge, buffer, type) result(error_status)
  INTEGER(kind=int4), intent(in) :: sd_id
  CHARACTER(*), intent(in) :: sds_name
  INTEGER(kind=int4), dimension(4), intent(inout) :: istart, istride
  INTEGER(kind=int4), dimension(4), intent(inout) :: iedge
  INTEGER(kind=int1), dimension(:,:,:,:), allocatable, intent(inout) :: buffer
  INTEGER(kind=int4), intent(out), optional :: type
  INTEGER(kind=int4) :: astatus, istatus, sds_id, sds_rank, sds_type, sds_nattr, error_status, idim
  INTEGER(kind=int4), dimension(4) :: sds_dims, max_iedge
    
  integer :: sfselect, sfn2index, sfrdata, sfendacc, sfginfo
  
  error_status = sym%SUCCESS
  
  sds_id = sfselect(sd_id, sfn2index(sd_id,trim(sds_name)))
  if (sds_id == FAIL) then
    print "(a,'Error selecting ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    error_status = sym%FAILURE
    return
  endif
    
  istatus = sfginfo(sds_id, trim(sds_name), sds_rank, sds_dims, sds_type, sds_nattr)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (present(type)) then
    type = sds_type
    istatus = sfendacc(sds_id)
    error_status = sym%SUCCESS
    return
  endif
  
  if (sds_type /= DFNT_INT8 .and. &
      sds_type /= DFNT_UINT8 .and. &
      sds_type /= DFNT_CHAR8 .and. &
      sds_type /= DFNT_CHAR) then
    print "(a,'Error reading (type mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (sds_rank /= 4) then
    print "(a,'Error reading (rank mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  do idim=1, sds_rank
  
    if (istart(idim) < 0) istart(idim) = 0
    if (istride(idim) < 1) istride(idim) = 1
  
    if (iedge(idim) < 0) iedge(idim) = int(ceiling(real(sds_dims(idim))/real(istride(idim))))
    max_iedge(idim) = int(ceiling(real(sds_dims(idim) - istart(idim))/real(istride(idim))))
    iedge(idim) = min(max_iedge(idim),iedge(idim))
    
  end do
    
  if (allocated(buffer)) then
    if (size(buffer,1) < iedge(1) .or. &
        size(buffer,2) < iedge(2) .or. &
size(buffer,3) < iedge(3) .or. &
size(buffer,4) < iedge(4)) then
      deallocate(buffer,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 4d HDF buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif
  
  if ((.not. allocated(buffer))) then
    allocate(buffer(iedge(1),iedge(2),iedge(3),iedge(4)),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 4d buffer.')",EXE_PROMPT
      stop
    endif
  endif
    
  istatus = sfrdata(sds_id, istart, istride, iedge, buffer)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  istatus = sfendacc(sds_id)
  
  return
    
END FUNCTION read_hdf_sds_int8_4d

!-------------------------------------------------------------------
! Subroutine to read 4D int16 hdf data.
!-------------------------------------------------------------------
  
FUNCTION read_hdf_sds_int16_4d(sd_id, sds_name, istart, istride, iedge, buffer, type) result(error_status)
  INTEGER(kind=int4), intent(in) :: sd_id
  CHARACTER(*), intent(in) :: sds_name
  INTEGER(kind=int4), dimension(4), intent(inout) :: istart, istride
  INTEGER(kind=int4), dimension(4), intent(inout) :: iedge
  INTEGER(kind=int2), dimension(:,:,:,:), allocatable, intent(inout) :: buffer
  INTEGER(kind=int4), intent(out), optional :: type
  INTEGER(kind=int4) :: astatus, istatus, sds_id, sds_rank, sds_type, sds_nattr, error_status, idim
  INTEGER(kind=int4), dimension(4) :: sds_dims, max_iedge
    
  integer :: sfselect, sfn2index, sfrdata, sfendacc, sfginfo
  
  error_status = sym%SUCCESS
  
  sds_id = sfselect(sd_id, sfn2index(sd_id,trim(sds_name)))
  if (sds_id == FAIL) then
    print "(a,'Error selecting ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    error_status = sym%FAILURE
    return
  endif
    
  istatus = sfginfo(sds_id, trim(sds_name), sds_rank, sds_dims, sds_type, sds_nattr)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (present(type)) then
    type = sds_type
    istatus = sfendacc(sds_id)
    error_status = sym%SUCCESS
    return
  endif
  
  if (sds_type /= DFNT_INT16) then
    print "(a,'Error reading (type mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (sds_rank /= 4) then
    print "(a,'Error reading (rank mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  do idim=1, sds_rank
  
    if (istart(idim) < 0) istart(idim) = 0
    if (istride(idim) < 1) istride(idim) = 1
  
    if (iedge(idim) < 0) iedge(idim) = int(ceiling(real(sds_dims(idim))/real(istride(idim))))
    max_iedge(idim) = int(ceiling(real(sds_dims(idim) - istart(idim))/real(istride(idim))))
    iedge(idim) = min(max_iedge(idim),iedge(idim))
    
  end do
    
  if (allocated(buffer)) then
    if (size(buffer,1) < iedge(1) .or. &
        size(buffer,2) < iedge(2) .or. &
size(buffer,3) < iedge(3) .or. &
size(buffer,4) < iedge(4)) then
      deallocate(buffer,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 4d HDF buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif
  
  if ((.not. allocated(buffer))) then
    allocate(buffer(iedge(1),iedge(2),iedge(3),iedge(4)),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 4d buffer.')",EXE_PROMPT
      stop
    endif
  endif
    
  istatus = sfrdata(sds_id, istart, istride, iedge, buffer)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  istatus = sfendacc(sds_id)
  
  return
    
END FUNCTION read_hdf_sds_int16_4d

!-------------------------------------------------------------------
! Subroutine to read 4D int32 hdf data.
!-------------------------------------------------------------------
  
FUNCTION read_hdf_sds_int32_4d(sd_id, sds_name, istart, istride, iedge, buffer, type) result(error_status)
  INTEGER(kind=int4), intent(in) :: sd_id
  CHARACTER(*), intent(in) :: sds_name
  INTEGER(kind=int4), dimension(4), intent(inout) :: istart, istride
  INTEGER(kind=int4), dimension(4), intent(inout) :: iedge
  INTEGER(kind=int4), dimension(:,:,:,:), allocatable, intent(inout) :: buffer
  INTEGER(kind=int4), intent(out), optional :: type
  INTEGER(kind=int4) :: astatus, istatus, sds_id, sds_rank, sds_type, sds_nattr, error_status, idim
  INTEGER(kind=int4), dimension(4) :: sds_dims, max_iedge
    
  integer :: sfselect, sfn2index, sfrdata, sfendacc, sfginfo
  
  error_status = sym%SUCCESS
  
  sds_id = sfselect(sd_id, sfn2index(sd_id,trim(sds_name)))
  if (sds_id == FAIL) then
    print "(a,'Error selecting ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    error_status = sym%FAILURE
    return
  endif
    
  istatus = sfginfo(sds_id, trim(sds_name), sds_rank, sds_dims, sds_type, sds_nattr)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (present(type)) then
    type = sds_type
    istatus = sfendacc(sds_id)
    error_status = sym%SUCCESS
    return
  endif
  
  if (sds_type /= DFNT_INT32) then
    print "(a,'Error reading (type mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (sds_rank /= 4) then
    print "(a,'Error reading (rank mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  do idim=1, sds_rank
  
    if (istart(idim) < 0) istart(idim) = 0
    if (istride(idim) < 1) istride(idim) = 1
  
    if (iedge(idim) < 0) iedge(idim) = int(ceiling(real(sds_dims(idim))/real(istride(idim))))
    max_iedge(idim) = int(ceiling(real(sds_dims(idim) - istart(idim))/real(istride(idim))))
    iedge(idim) = min(max_iedge(idim),iedge(idim))
    
  end do
    
  if (allocated(buffer)) then
    if (size(buffer,1) < iedge(1) .or. &
        size(buffer,2) < iedge(2) .or. &
size(buffer,3) < iedge(3) .or. &
size(buffer,4) < iedge(4)) then
      deallocate(buffer,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 4d HDF buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif
  
  if ((.not. allocated(buffer))) then
    allocate(buffer(iedge(1),iedge(2),iedge(3),iedge(4)),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 4d buffer.')",EXE_PROMPT
      stop
    endif
  endif
    
  istatus = sfrdata(sds_id, istart, istride, iedge, buffer)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  istatus = sfendacc(sds_id)
  
  return
    
END FUNCTION read_hdf_sds_int32_4d

!-------------------------------------------------------------------
! Subroutine to read 4D float32 hdf data.
!-------------------------------------------------------------------
  
FUNCTION read_hdf_sds_float32_4d(sd_id, sds_name, istart, istride, iedge, buffer, type) result(error_status)
  INTEGER(kind=int4), intent(in) :: sd_id
  CHARACTER(*), intent(in) :: sds_name
  INTEGER(kind=int4), dimension(4), intent(inout) :: istart, istride
  INTEGER(kind=int4), dimension(4), intent(inout) :: iedge
  REAL(kind=real4), dimension(:,:,:,:), allocatable, intent(inout) :: buffer
  INTEGER(kind=int4), intent(out), optional :: type
  INTEGER(kind=int4) :: astatus, istatus, sds_id, sds_rank, sds_type, sds_nattr, error_status, idim
  INTEGER(kind=int4), dimension(4) :: sds_dims, max_iedge
    
  integer :: sfselect, sfn2index, sfrdata, sfendacc, sfginfo
  
  error_status = sym%SUCCESS
  
  sds_id = sfselect(sd_id, sfn2index(sd_id,trim(sds_name)))
  if (sds_id == FAIL) then
    print "(a,'Error selecting ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    error_status = sym%FAILURE
    return
  endif
    
  istatus = sfginfo(sds_id, trim(sds_name), sds_rank, sds_dims, sds_type, sds_nattr)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (present(type)) then
    type = sds_type
    istatus = sfendacc(sds_id)
    error_status = sym%SUCCESS
    return
  endif
  
  if (sds_type /= DFNT_FLOAT32) then
    print "(a,'Error reading (type mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (sds_rank /= 4) then
    print "(a,'Error reading (rank mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  do idim=1, sds_rank
  
    if (istart(idim) < 0) istart(idim) = 0
    if (istride(idim) < 1) istride(idim) = 1
  
    if (iedge(idim) < 0) iedge(idim) = int(ceiling(real(sds_dims(idim))/real(istride(idim))))
    max_iedge(idim) = int(ceiling(real(sds_dims(idim) - istart(idim))/real(istride(idim))))
    iedge(idim) = min(max_iedge(idim),iedge(idim))
    
  end do
  
  if (allocated(buffer)) then
    if (size(buffer,1) < iedge(1) .or. &
        size(buffer,2) < iedge(2) .or. &
size(buffer,3) < iedge(3) .or. &
size(buffer,4) < iedge(4)) then
      deallocate(buffer,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 4d HDF buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif
  
  if ((.not. allocated(buffer))) then
    allocate(buffer(iedge(1),iedge(2),iedge(3),iedge(4)),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 4d buffer.')",EXE_PROMPT
      stop
    endif
  endif
    
  istatus = sfrdata(sds_id, istart, istride, iedge, buffer)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  istatus = sfendacc(sds_id)
  
  return
    
END FUNCTION read_hdf_sds_float32_4d

!-------------------------------------------------------------------
! Subroutine to read 4D float64 hdf data.
!-------------------------------------------------------------------
  
FUNCTION read_hdf_sds_float64_4d(sd_id, sds_name, istart, istride, iedge, buffer, type) result(error_status)
  INTEGER(kind=int4), intent(in) :: sd_id
  CHARACTER(*), intent(in) :: sds_name
  INTEGER(kind=int4), dimension(4), intent(inout) :: istart, istride
  INTEGER(kind=int4), dimension(4), intent(inout) :: iedge
  REAL(kind=real8), dimension(:,:,:,:), allocatable, intent(inout) :: buffer
  INTEGER(kind=int4), intent(out), optional :: type
  INTEGER(kind=int4) :: astatus, istatus, sds_id, sds_rank, sds_type, sds_nattr, error_status, idim
  INTEGER(kind=int4), dimension(4) :: sds_dims, max_iedge
    
  integer :: sfselect, sfn2index, sfrdata, sfendacc, sfginfo
  
  error_status = sym%SUCCESS
  
  sds_id = sfselect(sd_id, sfn2index(sd_id,trim(sds_name)))
  if (sds_id == FAIL) then
    print "(a,'Error selecting ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    error_status = sym%FAILURE
    return
  endif
    
  istatus = sfginfo(sds_id, trim(sds_name), sds_rank, sds_dims, sds_type, sds_nattr)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (present(type)) then
    type = sds_type
    istatus = sfendacc(sds_id)
    error_status = sym%SUCCESS
    return
  endif
  
  if (sds_type /= DFNT_FLOAT64) then
    print "(a,'Error reading (type mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (sds_rank /= 4) then
    print "(a,'Error reading (rank mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  do idim=1, sds_rank
  
    if (istart(idim) < 0) istart(idim) = 0
    if (istride(idim) < 1) istride(idim) = 1
  
    if (iedge(idim) < 0) iedge(idim) = int(ceiling(real(sds_dims(idim))/real(istride(idim))))
    max_iedge(idim) = int(ceiling(real(sds_dims(idim) - istart(idim))/real(istride(idim))))
    iedge(idim) = min(max_iedge(idim),iedge(idim))
    
  end do
    
  if (allocated(buffer)) then
    if (size(buffer,1) < iedge(1) .or. &
        size(buffer,2) < iedge(2) .or. &
size(buffer,3) < iedge(3) .or. &
size(buffer,4) < iedge(4)) then
      deallocate(buffer,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 4d HDF buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif
  
  if ((.not. allocated(buffer))) then
    allocate(buffer(iedge(1),iedge(2),iedge(3),iedge(4)),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 4d buffer.')",EXE_PROMPT
      stop
    endif
  endif
    
  istatus = sfrdata(sds_id, istart, istride, iedge, buffer)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  istatus = sfendacc(sds_id)
  
  return
    
END FUNCTION read_hdf_sds_float64_4d





!-------------------------------------------------------------------
! Subroutine to read 5D int8 hdf data.
!-------------------------------------------------------------------
  
FUNCTION read_hdf_sds_int8_5d(sd_id, sds_name, istart, istride, iedge, buffer, type) result(error_status)
  INTEGER(kind=int4), intent(in) :: sd_id
  CHARACTER(*), intent(in) :: sds_name
  INTEGER(kind=int4), dimension(5), intent(inout) :: istart, istride
  INTEGER(kind=int4), dimension(5), intent(inout) :: iedge
  INTEGER(kind=int1), dimension(:,:,:,:,:), allocatable, intent(inout) :: buffer
  INTEGER(kind=int4), intent(out), optional :: type
  INTEGER(kind=int4) :: astatus, istatus, sds_id, sds_rank, sds_type, sds_nattr, error_status, idim
  INTEGER(kind=int4), dimension(5) :: sds_dims, max_iedge
    
  integer :: sfselect, sfn2index, sfrdata, sfendacc, sfginfo
  
  error_status = sym%SUCCESS
  
  sds_id = sfselect(sd_id, sfn2index(sd_id,trim(sds_name)))
  if (sds_id == FAIL) then
    print "(a,'Error selecting ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    error_status = sym%FAILURE
    return
  endif
    
  istatus = sfginfo(sds_id, trim(sds_name), sds_rank, sds_dims, sds_type, sds_nattr)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (present(type)) then
    type = sds_type
    istatus = sfendacc(sds_id)
    error_status = sym%SUCCESS
    return
  endif
  
  if (sds_type /= DFNT_INT8 .and. &
      sds_type /= DFNT_UINT8 .and. &
      sds_type /= DFNT_CHAR8 .and. &
      sds_type /= DFNT_CHAR) then
    print "(a,'Error reading (type mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (sds_rank /= 5) then
    print "(a,'Error reading (rank mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  do idim=1, sds_rank
  
    if (istart(idim) < 0) istart(idim) = 0
    if (istride(idim) < 1) istride(idim) = 1
  
    if (iedge(idim) < 0) iedge(idim) = int(ceiling(real(sds_dims(idim))/real(istride(idim))))
    max_iedge(idim) = int(ceiling(real(sds_dims(idim) - istart(idim))/real(istride(idim))))
    iedge(idim) = min(max_iedge(idim),iedge(idim))
    
  end do
    
  if (allocated(buffer)) then
    if (size(buffer,1) < iedge(1) .or. &
        size(buffer,2) < iedge(2) .or. &
size(buffer,3) < iedge(3) .or. &
size(buffer,4) < iedge(4) .or. &
        size(buffer,5) < iedge(5)) then
      deallocate(buffer,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 4d HDF buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif
  
  if ((.not. allocated(buffer))) then
    allocate(buffer(iedge(1),iedge(2),iedge(3),iedge(4),iedge(5)),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 5d buffer.')",EXE_PROMPT
      stop
    endif
  endif
    
  istatus = sfrdata(sds_id, istart, istride, iedge, buffer)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  istatus = sfendacc(sds_id)
  
  return
    
END FUNCTION read_hdf_sds_int8_5d

!-------------------------------------------------------------------
! Subroutine to read 5D int16 hdf data.
!-------------------------------------------------------------------
  
FUNCTION read_hdf_sds_int16_5d(sd_id, sds_name, istart, istride, iedge, buffer, type) result(error_status)
  INTEGER(kind=int4), intent(in) :: sd_id
  CHARACTER(*), intent(in) :: sds_name
  INTEGER(kind=int4), dimension(5), intent(inout) :: istart, istride
  INTEGER(kind=int4), dimension(5), intent(inout) :: iedge
  INTEGER(kind=int2), dimension(:,:,:,:,:), allocatable, intent(inout) :: buffer
  INTEGER(kind=int4), intent(out), optional :: type
  INTEGER(kind=int4) :: astatus, istatus, sds_id, sds_rank, sds_type, sds_nattr, error_status, idim
  INTEGER(kind=int4), dimension(5) :: sds_dims, max_iedge
    
  integer :: sfselect, sfn2index, sfrdata, sfendacc, sfginfo
  
  error_status = sym%SUCCESS
  
  sds_id = sfselect(sd_id, sfn2index(sd_id,trim(sds_name)))
  if (sds_id == FAIL) then
    print "(a,'Error selecting ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    error_status = sym%FAILURE
    return
  endif
    
  istatus = sfginfo(sds_id, trim(sds_name), sds_rank, sds_dims, sds_type, sds_nattr)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (present(type)) then
    type = sds_type
    istatus = sfendacc(sds_id)
    error_status = sym%SUCCESS
    return
  endif
  
  if (sds_type /= DFNT_INT16) then
    print "(a,'Error reading (type mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (sds_rank /= 5) then
    print "(a,'Error reading (rank mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  do idim=1, sds_rank
  
    if (istart(idim) < 0) istart(idim) = 0
    if (istride(idim) < 1) istride(idim) = 1
  
    if (iedge(idim) < 0) iedge(idim) = int(ceiling(real(sds_dims(idim))/real(istride(idim))))
    max_iedge(idim) = int(ceiling(real(sds_dims(idim) - istart(idim))/real(istride(idim))))
    iedge(idim) = min(max_iedge(idim),iedge(idim))
    
  end do
    
  if (allocated(buffer)) then
    if (size(buffer,1) < iedge(1) .or. &
        size(buffer,2) < iedge(2) .or. &
size(buffer,3) < iedge(3) .or. &
size(buffer,4) < iedge(4) .or. &
        size(buffer,5) < iedge(5)) then
      deallocate(buffer,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 4d HDF buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif
  
  if ((.not. allocated(buffer))) then
    allocate(buffer(iedge(1),iedge(2),iedge(3),iedge(4),iedge(5)),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 5d buffer.')",EXE_PROMPT
      stop
    endif
  endif
    
  istatus = sfrdata(sds_id, istart, istride, iedge, buffer)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  istatus = sfendacc(sds_id)
  
  return
    
END FUNCTION read_hdf_sds_int16_5d

!-------------------------------------------------------------------
! Subroutine to read 5D int32 hdf data.
!-------------------------------------------------------------------
  
FUNCTION read_hdf_sds_int32_5d(sd_id, sds_name, istart, istride, iedge, buffer, type) result(error_status)
  INTEGER(kind=int4), intent(in) :: sd_id
  CHARACTER(*), intent(in) :: sds_name
  INTEGER(kind=int4), dimension(5), intent(inout) :: istart, istride
  INTEGER(kind=int4), dimension(5), intent(inout) :: iedge
  INTEGER(kind=int4), dimension(:,:,:,:,:), allocatable, intent(inout) :: buffer
  INTEGER(kind=int4), intent(out), optional :: type
  INTEGER(kind=int4) :: astatus, istatus, sds_id, sds_rank, sds_type, sds_nattr, error_status, idim
  INTEGER(kind=int4), dimension(5) :: sds_dims, max_iedge
    
  integer :: sfselect, sfn2index, sfrdata, sfendacc, sfginfo
  
  error_status = sym%SUCCESS
  
  sds_id = sfselect(sd_id, sfn2index(sd_id,trim(sds_name)))
  if (sds_id == FAIL) then
    print "(a,'Error selecting ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    error_status = sym%FAILURE
    return
  endif
    
  istatus = sfginfo(sds_id, trim(sds_name), sds_rank, sds_dims, sds_type, sds_nattr)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (present(type)) then
    type = sds_type
    istatus = sfendacc(sds_id)
    error_status = sym%SUCCESS
    return
  endif
  
  if (sds_type /= DFNT_INT32) then
    print "(a,'Error reading (type mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (sds_rank /= 5) then
    print "(a,'Error reading (rank mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  do idim=1, sds_rank
  
    if (istart(idim) < 0) istart(idim) = 0
    if (istride(idim) < 1) istride(idim) = 1
  
    if (iedge(idim) < 0) iedge(idim) = int(ceiling(real(sds_dims(idim))/real(istride(idim))))
    max_iedge(idim) = int(ceiling(real(sds_dims(idim) - istart(idim))/real(istride(idim))))
    iedge(idim) = min(max_iedge(idim),iedge(idim))
    
  end do
    
  if (allocated(buffer)) then
    if (size(buffer,1) < iedge(1) .or. &
        size(buffer,2) < iedge(2) .or. &
size(buffer,3) < iedge(3) .or. &
size(buffer,4) < iedge(4) .or. &
        size(buffer,5) < iedge(5)) then
      deallocate(buffer,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 5d HDF buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif
  
  if ((.not. allocated(buffer))) then
    allocate(buffer(iedge(1),iedge(2),iedge(3),iedge(4),iedge(5)),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 5d buffer.')",EXE_PROMPT
      stop
    endif
  endif
    
  istatus = sfrdata(sds_id, istart, istride, iedge, buffer)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  istatus = sfendacc(sds_id)
  
  return
    
END FUNCTION read_hdf_sds_int32_5d

!-------------------------------------------------------------------
! Subroutine to read 5D float32 hdf data.
!-------------------------------------------------------------------
  
FUNCTION read_hdf_sds_float32_5d(sd_id, sds_name, istart, istride, iedge, buffer, type) result(error_status)
  INTEGER(kind=int4), intent(in) :: sd_id
  CHARACTER(*), intent(in) :: sds_name
  INTEGER(kind=int4), dimension(5), intent(inout) :: istart, istride
  INTEGER(kind=int4), dimension(5), intent(inout) :: iedge
  REAL(kind=real4), dimension(:,:,:,:,:), allocatable, intent(inout) :: buffer
  INTEGER(kind=int4), intent(out), optional :: type
  INTEGER(kind=int4) :: astatus, istatus, sds_id, sds_rank, sds_type, sds_nattr, error_status, idim
  INTEGER(kind=int4), dimension(5) :: sds_dims, max_iedge
    
  integer :: sfselect, sfn2index, sfrdata, sfendacc, sfginfo
  
  error_status = sym%SUCCESS
  
  sds_id = sfselect(sd_id, sfn2index(sd_id,trim(sds_name)))
  if (sds_id == FAIL) then
    print "(a,'Error selecting ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    error_status = sym%FAILURE
    return
  endif
    
  istatus = sfginfo(sds_id, trim(sds_name), sds_rank, sds_dims, sds_type, sds_nattr)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (present(type)) then
    type = sds_type
    istatus = sfendacc(sds_id)
    error_status = sym%SUCCESS
    return
  endif
  
  if (sds_type /= DFNT_FLOAT32) then
    print "(a,'Error reading (type mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (sds_rank /= 5) then
    print "(a,'Error reading (rank mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  do idim=1, sds_rank
  
    if (istart(idim) < 0) istart(idim) = 0
    if (istride(idim) < 1) istride(idim) = 1
  
    if (iedge(idim) < 0) iedge(idim) = int(ceiling(real(sds_dims(idim))/real(istride(idim))))
    max_iedge(idim) = int(ceiling(real(sds_dims(idim) - istart(idim))/real(istride(idim))))
    iedge(idim) = min(max_iedge(idim),iedge(idim))
    
  end do
  
  if (allocated(buffer)) then
    if (size(buffer,1) < iedge(1) .or. &
        size(buffer,2) < iedge(2) .or. &
size(buffer,3) < iedge(3) .or. &
size(buffer,4) < iedge(4) .or. &
        size(buffer,5) < iedge(5)) then
      deallocate(buffer,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 5d HDF buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif
  
  if ((.not. allocated(buffer))) then
    allocate(buffer(iedge(1),iedge(2),iedge(3),iedge(4),iedge(5)),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 5d buffer.')",EXE_PROMPT
      stop
    endif
  endif
    
  istatus = sfrdata(sds_id, istart, istride, iedge, buffer)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  istatus = sfendacc(sds_id)
  
  return
    
END FUNCTION read_hdf_sds_float32_5d

!-------------------------------------------------------------------
! Subroutine to read 5D float64 hdf data.
!-------------------------------------------------------------------
  
FUNCTION read_hdf_sds_float64_5d(sd_id, sds_name, istart, istride, iedge, buffer, type) result(error_status)
  INTEGER(kind=int4), intent(in) :: sd_id
  CHARACTER(*), intent(in) :: sds_name
  INTEGER(kind=int4), dimension(5), intent(inout) :: istart, istride
  INTEGER(kind=int4), dimension(5), intent(inout) :: iedge
  REAL(kind=real8), dimension(:,:,:,:,:), allocatable, intent(inout) :: buffer
  INTEGER(kind=int4), intent(out), optional :: type
  INTEGER(kind=int4) :: astatus, istatus, sds_id, sds_rank, sds_type, sds_nattr, error_status, idim
  INTEGER(kind=int4), dimension(5) :: sds_dims, max_iedge
    
  integer :: sfselect, sfn2index, sfrdata, sfendacc, sfginfo
  
  error_status = sym%SUCCESS
  
  sds_id = sfselect(sd_id, sfn2index(sd_id,trim(sds_name)))
  if (sds_id == FAIL) then
    print "(a,'Error selecting ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    error_status = sym%FAILURE
    return
  endif
    
  istatus = sfginfo(sds_id, trim(sds_name), sds_rank, sds_dims, sds_type, sds_nattr)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (present(type)) then
    type = sds_type
    istatus = sfendacc(sds_id)
    error_status = sym%SUCCESS
    return
  endif
  
  if (sds_type /= DFNT_FLOAT64) then
    print "(a,'Error reading (type mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (sds_rank /= 5) then
    print "(a,'Error reading (rank mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  do idim=1, sds_rank
  
    if (istart(idim) < 0) istart(idim) = 0
    if (istride(idim) < 1) istride(idim) = 1
  
    if (iedge(idim) < 0) iedge(idim) = int(ceiling(real(sds_dims(idim))/real(istride(idim))))
    max_iedge(idim) = int(ceiling(real(sds_dims(idim) - istart(idim))/real(istride(idim))))
    iedge(idim) = min(max_iedge(idim),iedge(idim))
    
  end do
    
  if (allocated(buffer)) then
    if (size(buffer,1) < iedge(1) .or. &
        size(buffer,2) < iedge(2) .or. &
size(buffer,3) < iedge(3) .or. &
size(buffer,4) < iedge(4) .or. &
        size(buffer,5) < iedge(5)) then
      deallocate(buffer,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 5d HDF buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif
  
  if ((.not. allocated(buffer))) then
    allocate(buffer(iedge(1),iedge(2),iedge(3),iedge(4),iedge(5)),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 4d buffer.')",EXE_PROMPT
      stop
    endif
  endif
    
  istatus = sfrdata(sds_id, istart, istride, iedge, buffer)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  istatus = sfendacc(sds_id)
  
  return
    
END FUNCTION read_hdf_sds_float64_5d

!-------------------------------------------------------------------
! Subroutine to read 6D float64 hdf data.
!-------------------------------------------------------------------
  
FUNCTION read_hdf_sds_float64_6d(sd_id, sds_name, istart, istride, iedge, buffer, type) result(error_status)
  INTEGER(kind=int4), intent(in) :: sd_id
  CHARACTER(*), intent(in) :: sds_name
  INTEGER(kind=int4), dimension(6), intent(inout) :: istart, istride
  INTEGER(kind=int4), dimension(6), intent(inout) :: iedge
  REAL(kind=real8), dimension(:,:,:,:,:,:), allocatable, intent(inout) :: buffer
  INTEGER(kind=int4), intent(out), optional :: type
  INTEGER(kind=int4) :: astatus, istatus, sds_id, sds_rank, sds_type, sds_nattr, error_status, idim
  INTEGER(kind=int4), dimension(6) :: sds_dims, max_iedge
    
  integer :: sfselect, sfn2index, sfrdata, sfendacc, sfginfo
  
  error_status = sym%SUCCESS
  
  sds_id = sfselect(sd_id, sfn2index(sd_id,trim(sds_name)))
  if (sds_id == FAIL) then
    print "(a,'Error selecting ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    error_status = sym%FAILURE
    return
  endif
    
  istatus = sfginfo(sds_id, trim(sds_name), sds_rank, sds_dims, sds_type, sds_nattr)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (present(type)) then
    type = sds_type
    istatus = sfendacc(sds_id)
    error_status = sym%SUCCESS
    return
  endif
  
  if (sds_type /= DFNT_FLOAT64) then
    print "(a,'Error reading (type mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  if (sds_rank /= 6) then
    print "(a,'Error reading (rank mismatch) ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  
  do idim=1, sds_rank
  
    if (istart(idim) < 0) istart(idim) = 0
    if (istride(idim) < 1) istride(idim) = 1
  
    if (iedge(idim) < 0) iedge(idim) = int(ceiling(real(sds_dims(idim))/real(istride(idim))))
    max_iedge(idim) = int(ceiling(real(sds_dims(idim) - istart(idim))/real(istride(idim))))
    iedge(idim) = min(max_iedge(idim),iedge(idim))
    
  end do
    
  if (allocated(buffer)) then
    if (size(buffer,1) < iedge(1) .or. &
        size(buffer,2) < iedge(2) .or. &
size(buffer,3) < iedge(3) .or. &
size(buffer,4) < iedge(4) .or. &
        size(buffer,5) < iedge(5) .or. &
        size(buffer,6) < iedge(6)) then
      deallocate(buffer,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 6d HDF buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif
  
  if ((.not. allocated(buffer))) then
    allocate(buffer(iedge(1),iedge(2),iedge(3),iedge(4),iedge(5),iedge(6)),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 6d buffer.')",EXE_PROMPT
      stop
    endif
  endif
    
  istatus = sfrdata(sds_id, istart, istride, iedge, buffer)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    istatus = sfendacc(sds_id)
    error_status = sym%FAILURE
    return
  endif
  istatus = sfendacc(sds_id)
  
  return
    
END FUNCTION read_hdf_sds_float64_6d


!-------------------------------------------------------------------
! Subroutine to read in SDS dimensions.
!-------------------------------------------------------------------

FUNCTION hdf_sds_dimenions_reader(id, sds_name, rank, dims) result(error_status)
  INTEGER(kind=int4), intent(in) :: id
  CHARACTER(len=*), intent(in) :: sds_name
  INTEGER(kind=int4), intent(out) :: rank
  INTEGER(kind=int4), dimension(MAX_RANK_HDF), intent(out) :: dims
  
  INTEGER(kind=int4) :: error_status
  INTEGER(kind=int4) :: sds_id, istatus, sds_type, sds_nattr
    
  INTEGER :: sfselect, sfn2index, sfginfo, sfendacc
  
  error_status = sym%SUCCESS
  
  sds_id = sfselect(id, sfn2index(id,trim(sds_name)))
    
  istatus = sfginfo(sds_id, sds_name, rank, dims, sds_type, sds_nattr)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),id    
    error_status = sym%FAILURE
  endif
  
  istatus = sfendacc(sds_id)
  return
  
END FUNCTION hdf_sds_dimenions_reader






END MODULE hdf_read_module
