!$Id$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: SFC_EMISS.f90 (src)
!       SFC_EMISS (program)
!
! PURPOSE: Routines for opening, reading and closing the SEEBOR Emissivity database
!
! DESCRIPTION: 
!
! AUTHORS:
!  Andrew Heidinger, Andrew.Heidinger@noaa.gov
!  Andi Walther, CIMSS, andi.walther@ssec.wisc.edu
!  Denis Botambekov, CIMSS, denis.botambekov@ssec.wisc.edu
!  William Straka, CIMSS, wstraka@ssec.wisc.edu
!
! COPYRIGHT
! THIS SOFTWARE AND ITS DOCUMENTATION ARE CONSIDERED TO BE IN THE PUBLIC
! DOMAIN AND THUS ARE AVAILABLE FOR UNRESTRICTED PUBLIC USE. THEY ARE
! FURNISHED "AS IS." THE AUTHORS, THE UNITED STATES GOVERNMENT, ITS
! INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND AGENTS MAKE NO WARRANTY,
! EXPRESS OR IMPLIED, AS TO THE USEFULNESS OF THE SOFTWARE AND
! DOCUMENTATION FOR ANY PURPOSE. THEY ASSUME NO RESPONSIBILITY (1) FOR
! THE USE OF THE SOFTWARE AND DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL
! SUPPORT TO USERS.
!
!--------------------------------------------------------------------------------------
module SFC_EMISS
   use HDF
   use CONSTANTS
   use NUMERICAL_ROUTINES
  
   implicit none
  
   !--- routine access declaration
   private :: read_integrated_seebor_hdf
   public :: open_seebor_emiss, close_seebor_emiss, read_seebor_emiss

   !---------------------------------------------------------------------------------------
   INTEGER, parameter, private :: Num_Lat_Emiss = 3600
   INTEGER, parameter, private :: Num_Lon_Emiss = 7200
   REAL(kind=real4), parameter, private :: First_Lat_Emiss = 89.9750, last_lat_emiss = -89.9750
   REAL(kind=real4), parameter, private :: First_Lon_Emiss = -179.975, last_lon_emiss = 179.975
   REAL(kind=real4), parameter, private :: Del_Lat_Emiss = 0.05
   REAL(kind=real4), parameter, private :: Del_Lon_Emiss = 0.05
CONTAINS
  
   !====================================================================
   ! subroutine Name: open_seebor_emiss
   !
   ! Function:
   !   Determines which SEEBOR Emissivty file to use
   !
   ! Description:
   !   This subroutine, given the ancillary data directory and month
   !   outputs the HDF file id associated with the correct SEEBOR
   !   Emmissivity database file
   !====================================================================

   subroutine open_seebor_Emiss(data_dir, month, id)
      CHARACTER(len=*), intent(in) :: data_dir
      INTEGER(kind=int2), intent(in) :: month
      INTEGER(kind=int4), intent(out) :: id
  
      CHARACTER(len=256) :: filename
      CHARACTER(len=3) :: jday_str
      CHARACTER(len=4) :: year_str
  
  logical :: file_exists
  
  INTEGER :: sfstart
  
  year_str = "2005"
  
  select case (month)
  case (1)
    jday_str = "001"
  case (2)
    jday_str = "032"
  case (3)
    jday_str = "060"
  case (4)
    jday_str = "091"
  case (5)
    jday_str = "121"
  case (6)
    jday_str = "152"
  case (7)
    jday_str = "182"
  case (8)
    jday_str = "213"
  case (9)
    jday_str = "244"
  case (10)
    jday_str = "274"
  case (11)
    jday_str = "305"
  case (12)
    jday_str = "335"
  end select
  
  filename = trim(data_dir)//"/global_emiss_intABI_"//trim(year_str)//trim(jday_str)//".hdf"
  
  inquire(file = filename, exist = file_exists)
  if (.not. file_exists) then
    print "(/,a,'Surface emissivity file, ',a,' does not exist.')",EXE_PROMPT,trim(filename)
    stop
  endif
  
  id = sfstart(trim(filename), DFACC_READ)
  if (id == FAIL) then
    print "(/,a,'Failed to open, ',a)",EXE_PROMPT,trim(filename)
    stop
  endif

end subroutine open_seebor_emiss


!====================================================================
! subroutine Name: close_seebor_emiss
!
! Function:
!   Closes a given SEEBOR Emissivity file
!
!====================================================================

subroutine close_seebor_Emiss(id)
  INTEGER(kind=int4), intent(in) :: id
  
  INTEGER(kind=int4) :: istatus
  INTEGER :: sfend
  
  istatus = sfend(id)
  if (istatus /= 0) then
    print "(/,a,'Error closing surface emissivity hdf file.')",EXE_PROMPT
    stop
  endif

end subroutine close_seebor_emiss

   !====================================================================
   ! subroutine Name: read_seebor_emiss
   !
   ! Function:
   !  Read a given channel from the surface emissivity file for a segment of
   !  Data
   ! 
   ! Description:
   !  This subroutine, given the latitude, longitude, space mask and channel
   !   number reads in a segment of data from the appropriate SEEBOR Emissivity file
   !   and outputs to the appropriate global array
   !====================================================================
  
   subroutine read_seebor_Emiss(id, ichan, lat, lon, Space_Mask, emiss)
   
      INTEGER(kind=int4), intent(in) :: id, ichan
      REAL(kind=real4), dimension(:,:), intent(in) :: lat, lon
      INTEGER(kind=int1), dimension(:,:), intent(in) :: Space_Mask
      REAL(kind=real4), dimension(:,:), intent(out) :: emiss
  
      CHARACTER (len=100) :: sds_name  
      INTEGER :: astatus
      INTEGER :: ilat1, ilat2, ilon1, ilon2, ilat, ilon, ilat_ad, ilon_ad, &
             ilon1_2, ilon2_2
      INTEGER :: temp, nx, ny, i, j
      REAL(kind=real4), dimension(:,:), allocatable :: Emiss_Grid, Emiss_Grid_2
      REAL(kind=real4) :: wlon, elon, slat, nlat
      INTEGER(kind=int1) :: dateline_flg, space_check
  
      INTEGER, dimension(2) :: start_2d, stride_2d, edge_2d, &
                           start_2d_2, stride_2d_2, edge_2d_2
  
      ! if (ichan < 3 .or. ichan > 5) then
      if (ichan < 7) then
         print "(a,'Surface emissivity: invalid channel number ',i0,' - cannot read in surface emissivity')",EXE_PROMPT,ichan
         stop
      end if
  
      space_check = minval(Space_Mask)
      if (space_check == 1) then
         emiss = missing_value_real4
         return
      endif
  
      write(sds_name,'(a,i0)') 'emiss',ichan
      if (ichan == 3) then
         sds_name = trim(sds_name)//'b'
      endif 

      nx = size(emiss,1)
      ny = size(emiss,2)
  
      call find_bounds(lat,lon,wlon,elon,slat,nlat,dateline_flg)
     
      if (dateline_flg == 0) then
  
         ilat1 = max(1,min(Num_Lat_Emiss,int(abs(nlat - First_Lat_Emiss)/Del_Lat_Emiss) + 1))
         ilat2 = max(1,min(Num_Lat_Emiss,int(abs(slat - First_Lat_Emiss)/Del_Lat_Emiss) + 1))
  
         ilon1 = max(1,min(Num_Lon_Emiss,int(abs(wlon - First_Lon_Emiss)/Del_Lon_Emiss) + 1))
         ilon2 = max(1,min(Num_Lon_Emiss,int(abs(elon - First_Lon_Emiss)/Del_Lon_Emiss) + 1))
  
         if (ilat1 > ilat2) then
            temp = ilat1
            ilat1 = ilat2
            ilat2 = temp
         endif
  
         if (ilon1 > ilon2) then
            temp = ilon1
            ilon1 = ilon2
            ilon2 = temp
         endif
  
         start_2d = (/ilon1, ilat1/)
         stride_2d = (/1, 1/)
         edge_2d = (/(ilon2-ilon1)+1, (ilat2-ilat1)+1/)
  
         call read_integrated_seebor_hdf(id, trim(sds_name), start_2d, stride_2d, edge_2d, Emiss_Grid)
  
         do j = 1, ny
            do i = 1, nx
    
               if (Space_Mask(i,j) == sym%NO_SPACE) then

                  ilat = max(1,min(Num_Lat_Emiss,int(abs(lat(i,j) - First_Lat_Emiss)/Del_Lat_Emiss) + 1))
                  ilon = max(1,min(Num_Lon_Emiss,int(abs(lon(i,j) - First_Lon_Emiss)/Del_Lon_Emiss) + 1))
             
                  ilat_ad = max(1,min((ilat - start_2d(2)) + 1,size(Emiss_Grid,2)))
                  ilon_ad = max(1,min((ilon - start_2d(1)) + 1,size(Emiss_Grid,1)))
                  if (Emiss_Grid(ilon_ad,ilat_ad) < 0.0) then
                     Emiss(i,j) = 0.99
                  else
                     Emiss(i,j) = Emiss_Grid(ilon_ad,ilat_ad)
                  end if

               end if      
            end do
         end do
    
         deallocate(Emiss_Grid, stat=astatus)
         if (astatus /= 0) then
            print "(a,'Error deallocating surface emissivity grid.')",EXE_PROMPT
            stop
         end if
    
      else
  
         ilat1 = max(1,min(Num_Lat_Emiss,int(abs(nlat - First_Lat_Emiss)/Del_Lat_Emiss) + 1))
         ilat2 = max(1,min(Num_Lat_Emiss,int(abs(slat - First_Lat_Emiss)/Del_Lat_Emiss) + 1))
  
         ilon1 = max(1,min(Num_Lon_Emiss,int(abs(wlon - First_Lon_Emiss)/Del_Lon_Emiss) + 1))
         ilon2 = max(1,min(Num_Lon_Emiss,int(abs(180.0 - First_Lon_Emiss)/Del_Lon_Emiss) + 1))
    
         ilon1_2 = max(1,min(Num_Lon_Emiss,int(abs(-180.0 - First_Lon_Emiss)/Del_Lon_Emiss) + 1))
         ilon2_2 = max(1,min(Num_Lon_Emiss,int(abs((elon-360.0) - First_Lon_Emiss)/Del_Lon_Emiss) + 1))
  
         if (ilat1 > ilat2) then
            temp = ilat1
            ilat1 = ilat2
            ilat2 = temp
         endif
  
         if (ilon1 > ilon2) then
            temp = ilon1
            ilon1 = ilon2
            ilon2 = temp
         endif
    
         if (ilon1_2 > ilon2_2) then
            temp = ilon1_2
            ilon1_2 = ilon2_2
            ilon2_2 = temp
         endif

         !--- prevent one single value in a hemisphere - akh mod
         ilon1_2 = 1  !this must be for a dateline crossing segment
         ilon2 = Num_Lon_Emiss   !this must be for dateline crossing segment
         if (ilon1 == Num_Lon_Emiss) then
            ilon1 = Num_Lon_Emiss -1
         endif
         if (ilon2_2 == 1) then
            ilon2_2 =  2
         endif
  
         !--- read western segment
         start_2d = (/ilon1, ilat1/)
         stride_2d = (/1, 1/)
         edge_2d = (/(ilon2-ilon1)+1, (ilat2-ilat1)+1/)

         call READ_INTEGRATED_SEEBOR_HDF(id, trim(sds_name), start_2d, stride_2d, edge_2d, Emiss_Grid)
    
         start_2d_2 = (/ilon1_2, ilat1/)
         stride_2d_2 = (/1, 1/)
         edge_2d_2 = (/(ilon2_2-ilon1_2)+1, (ilat2-ilat1)+1/)
  
         call READ_INTEGRATED_SEEBOR_HDF(id, trim(sds_name), start_2d_2, stride_2d_2, edge_2d_2, Emiss_Grid_2)
    
         do j = 1, ny
            do i = 1, nx
    
               if (Space_Mask(i,j) == sym%NO_SPACE) then

               ilat = max(1,min(Num_Lat_Emiss,int(abs(lat(i,j) - First_Lat_Emiss)/Del_Lat_Emiss) + 1))
               ilon = max(1,min(Num_Lon_Emiss,int(abs(lon(i,j) - First_Lon_Emiss)/Del_Lon_Emiss) + 1))
               if (lon(i,j) >= 0.0) then
                  ilat_ad = max(1,min((ilat - start_2d(2)) + 1,size(Emiss_Grid,2)))
                  ilon_ad = max(1,min((ilon - start_2d(1)) + 1,size(Emiss_Grid,1)))
                  Emiss(i,j) = Emiss_Grid(ilon_ad,ilat_ad)
               else
                  ilat_ad = max(1,min((ilat - start_2d_2(2)) + 1,size(Emiss_Grid_2,2)))
                  ilon_ad = max(1,min((ilon - start_2d_2(1)) + 1,size(Emiss_Grid_2,1)))
                  Emiss(i,j) = Emiss_Grid_2(ilon_ad,ilat_ad)
               endif
  
               if (Emiss(i,j) < 0.0) then
                  Emiss(i,j) = 0.99
               endif

            endif
      
         enddo
      enddo
  
      deallocate(Emiss_Grid, Emiss_Grid_2, stat=astatus)
      if (astatus /= 0) then
         print "(a,'Error deallocating surface emissivity grid.')",EXE_PROMPT
         stop
      endif
    
   endif
  
   end subroutine read_seebor_emiss

!====================================================================
! subroutine Name: read_integrated_seebor_hdf
!
! Function:
!  Read a given SDS from the surface emissivity file 
! 
! Description:
!  This subroutine, given the SDS name and segment information
!   number reads in the data from the HDF file.
!====================================================================
  
subroutine read_integrated_seebor_hdf(sd_id, sds_name, istart, istride, iedge, buffer)
  INTEGER(kind=int4), intent(in) :: sd_id
  CHARACTER(*), intent(in) :: sds_name
  INTEGER(kind=int4), dimension(2), intent(in) :: istart, istride
  INTEGER(kind=int4), dimension(2), intent(inout) :: iedge
  REAL(kind=real4), dimension(:,:), allocatable, intent(out) :: buffer
  INTEGER(kind=int2), dimension(:,:), allocatable :: int_buffer
  INTEGER(kind=int4) :: astatus, istatus, sds_id, sds_rank, sds_type, sds_nattr, attr_index
  INTEGER(kind=int4), dimension(2) :: sds_dims
  REAL(kind=real4), dimension(1) :: scale_fac, offset
    
  integer :: sfselect, sfn2index, sfrdata, sfendacc, sfginfo, sffattr, sfrnatt
  
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
  
  if (iedge(1) < 0) iedge(1) = sds_dims(1)
  if (iedge(2) < 0) iedge(2) = sds_dims(2)
  
  iedge(1) = min((sds_dims(1) - istart(1)),iedge(1))
  iedge(2) = min((sds_dims(2) - istart(2)),iedge(2))
    
  allocate(buffer(iedge(1),iedge(2)),int_buffer(iedge(1),iedge(2)),stat=astatus)
  if (astatus /= 0) then
    print "(a,'Not enough memory to allocate 2d buffer.')",EXE_PROMPT
    stop
  endif
    
  istatus = sfrdata(sds_id, istart, istride, iedge, int_buffer)
  if (istatus /= 0) then
    print "(a,'Error reading ',a,' from sd_id: ',i0)",EXE_PROMPT,trim(sds_name),sd_id
    stop
  endif
  istatus = sfendacc(sds_id)
  
  buffer = int_buffer*scale_fac(1) + offset(1)
  
  deallocate(int_buffer, stat=astatus)
  if (astatus /= 0) then
    print "(a,'Error deallocating emissivity integer buffer.')",EXE_PROMPT
    stop
  endif
    
end subroutine read_integrated_seebor_hdf
end module SFC_EMISS
