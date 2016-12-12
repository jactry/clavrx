!$Id$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: modis_module.f90 (src)
!       MODIS_MODULE (program)
!
! PURPOSE: 
!
! DESCRIPTION: 
!
! AUTHORS:
!  Andrew Heidinger, Andrew.Heidinger@noaa.gov
!  Andi Walther, CIMSS, andi.walther@ssec.wisc.edu
!  Denis Botambekov, CIMSS, denis.botambekov@ssec.wisc.edu
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
module MODIS_MODULE

        use CX_CONSTANTS_MOD, only: &
                int4, real4 &
                , int2 &
                , int1 &
                , exe_prompt &
                , Missing_Value_Real4 &
                , sym &
                , missing_value_int1

 
        use PIXEL_COMMON, only: &
                sensor &
                , geo &
                , image &
                , nav &
                , ch &
                , scan_time &
                , scan_number &
                , Cloud_Mask_Aux_Read_Flag &
                , Cloud_Mask_Aux_Flag &
                , Cld_mask_Aux &
                , Cld_Phase_Aux &
                , Cld_Type_Aux &
                , Zc_Aux &
                , Tau_Aux &
                , Reff_Aux &
                , Line_Idx_Min_segment

        use PIXEL_ROUTINES,only: &
               qc_modis
         
        use date_tools_mod, only: &
            julian 
         
        use PLANCK , only: &
                convert_radiance &
                , compute_bt_array

        use VIEWING_GEOMETRY_MODULE, only: &
                glint_angle &
                , scattering_angle
 
        use FILE_TOOLS, only: &
                file_search &
                , get_lun
 
        use CALIBRATION_CONSTANTS, only: &
                sat_name &
                , solar_ch20 &
                , solar_ch20_nu &
                , ew_ch20 &
                , planck_a1 &
                , planck_a2 &
                , planck_nu 

        implicit none
        private

        public:: DETERMINE_MODIS_GEOLOCATION_FILE
        public:: DETERMINE_MODIS_CLOUD_MASK_FILE
        public:: READ_MODIS
        public:: READ_MODIS_INSTR_CONSTANTS
        public:: READ_MODIS_TIME_ATTR
        public:: READ_MODIS_SIZE_ATTR

        private:: READ_MODIS_THERMAL_BAND

        integer(kind=int2), parameter, private:: fill_value = 32767
        real(kind=real4), parameter, private:: missing_value = -999.0
        character(len=13), parameter:: MODULE_PROMPT="MODIS_MODULE:"
        

       include 'hdf.f90'
contains

!----------------------------------------------------------------
! read the modis constants into memory
!-----------------------------------------------------------------
subroutine READ_MODIS_INSTR_CONSTANTS(Instr_Const_File)
 character(len=*), intent(in):: Instr_Const_File
 integer:: ios0, erstat
 integer:: Instr_Const_lun

 Instr_Const_lun = GET_LUN()

 open(unit=Instr_Const_lun,file=trim(Instr_Const_File),status="old",position="rewind",action="read",iostat=ios0)

 print *, EXE_PROMPT, MODULE_PROMPT, " Opening ", trim(Instr_Const_File)
 erstat = 0
 if (ios0 /= 0) then
    erstat = 19
    print *, EXE_PROMPT, MODULE_PROMPT, "Error opening MODIS constants file, ios0 = ", ios0
    stop 19
 endif

   read(unit=Instr_Const_lun,fmt="(a3)") sat_name
   read(unit=Instr_Const_lun,fmt=*) Solar_Ch20
   read(unit=Instr_Const_lun,fmt=*) Ew_Ch20
   read(unit=Instr_Const_lun,fmt=*) planck_a1(20), planck_a2(20), planck_nu(20)
   read(unit=Instr_Const_lun,fmt=*) planck_a1(21), planck_a2(21), planck_nu(21)
   read(unit=Instr_Const_lun,fmt=*) planck_a1(22), planck_a2(22), planck_nu(22)
   read(unit=Instr_Const_lun,fmt=*) planck_a1(23), planck_a2(23), planck_nu(23)
   read(unit=Instr_Const_lun,fmt=*) planck_a1(24), planck_a2(24), planck_nu(24)
   read(unit=Instr_Const_lun,fmt=*) planck_a1(25), planck_a2(25), planck_nu(25)
   read(unit=Instr_Const_lun,fmt=*) planck_a1(27), planck_a2(27), planck_nu(27)
   read(unit=Instr_Const_lun,fmt=*) planck_a1(28), planck_a2(28), planck_nu(28)
   read(unit=Instr_Const_lun,fmt=*) planck_a1(29), planck_a2(29), planck_nu(29)
   read(unit=Instr_Const_lun,fmt=*) planck_a1(30), planck_a2(30), planck_nu(30)
   read(unit=Instr_Const_lun,fmt=*) planck_a1(31), planck_a2(31), planck_nu(31)
   read(unit=Instr_Const_lun,fmt=*) planck_a1(32), planck_a2(32), planck_nu(32)
   read(unit=Instr_Const_lun,fmt=*) planck_a1(33), planck_a2(33), planck_nu(33)
   read(unit=Instr_Const_lun,fmt=*) planck_a1(34), planck_a2(34), planck_nu(34)
   read(unit=Instr_Const_lun,fmt=*) planck_a1(35), planck_a2(35), planck_nu(35)
   read(unit=Instr_Const_lun,fmt=*) planck_a1(36), planck_a2(36), planck_nu(36)
   close(unit=Instr_Const_lun)

   !-- convert solar flux in channel 20 to mean with units mW/m^2/cm^-1
   Solar_Ch20_Nu = 1000.0 * Solar_Ch20 / Ew_Ch20

end subroutine READ_MODIS_INSTR_CONSTANTS

!----------------------------------------------------------------------
subroutine DETERMINE_MODIS_GEOLOCATION_FILE(Modis_1b_Name &
            , Dir_Modis_Geo &
            , Auxiliary_Geolocation_File_Name)

    character(len=*), intent(in):: Modis_1b_Name
    character(len=*), intent(in):: Dir_Modis_Geo
    character(len=*), intent(out):: Auxiliary_Geolocation_File_Name
    integer, parameter:: nc = 24
    character(len=nc):: Search_String
    integer:: ilen
    integer(kind=int4):: Num_Files
    character(len=1020), dimension(:), pointer:: Files

    Search_String = trim(Modis_1b_Name(1:nc))
   
    if (trim(Sensor%Sensor_Name) == 'MODIS' .and. trim(Sensor%Platform_Name)  == 'AQUA') then
      Search_String = "MYD03"//trim(Search_String(9:nc) )//"*"
    else if (trim(Sensor%Sensor_Name) == 'MODIS-MAC') then
       Search_String = "MAC03S0"//trim(Search_String(9:nc))//"*"   
    else if (trim(Sensor%Sensor_Name) == 'MODIS-CSPP') then
       Search_String = trim(Modis_1b_Name(1:13)) //".geo.hdf"    
    else
       Search_String = "MOD03"//trim(Search_String(9:nc)//"*" )
    end if
    
    Files => FILE_SEARCH(Dir_Modis_Geo,Search_String)

    Num_Files = size(Files)

    if (Num_Files == 0) then 
       print *, EXE_PROMPT, MODULE_PROMPT, "No MODIS Geolocation File Found"
       Auxiliary_Geolocation_File_Name = "no_file"
       return
    end	if

    if (Num_Files > 1) then
       print *, EXE_PROMPT, MODULE_PROMPT, "Multiple MODIS Geolocation Files Found"
       Auxiliary_Geolocation_File_Name = "no_file"
    endif

    Auxiliary_Geolocation_File_Name = Files(1)

    Files => null()
  
    !--- if an AQUA Mac Name, trim something off
    ilen = len_trim(Auxiliary_Geolocation_File_Name)
    if (trim(Sensor%Sensor_Name) == 'MODIS-MAC') then
         Auxiliary_Geolocation_File_Name = Auxiliary_Geolocation_File_Name(ilen-42:ilen)
    else if (trim(Sensor%Sensor_Name) == 'MODIS-CSPP') then
         Auxiliary_Geolocation_File_Name = TRIM(Auxiliary_Geolocation_File_Name)
    else  
         Auxiliary_Geolocation_File_Name = Auxiliary_Geolocation_File_Name(ilen-40:ilen)
    endif

    print *, EXE_PROMPT, MODULE_PROMPT, "Will use MODIS Geolocation File = ", trim(Auxiliary_Geolocation_File_Name)
   
end subroutine DETERMINE_MODIS_GEOLOCATION_FILE

!----------------------------------------------------------------------
subroutine DETERMINE_MODIS_CLOUD_MASK_FILE(Modis_1b_Name, Dir_Modis_Cloud_Mask, Auxiliary_Cloud_Mask_File_Name)
    character(len=*), intent(in):: Modis_1b_Name
    character(len=*), intent(in):: Dir_Modis_Cloud_Mask
    character(len=*), intent(out):: Auxiliary_Cloud_Mask_File_Name
    integer, parameter:: nc = 25
    character(len=nc):: Search_String
    integer:: ilen
    integer(kind=int4):: Num_Files
    character(len=1020), dimension(:), pointer:: Files

    Search_String = trim(Modis_1b_Name(1:nc-1))

    if (Sensor%Spatial_Resolution_Meters == 1000) then
      if (trim(Sensor%Sensor_Name) == 'MODIS' .and. trim(Sensor%Platform_Name)  == 'AQUA') then
         Search_String = "MYD35_L2"//trim(Search_String(9:nc))//"*"
       else if (trim(Sensor%Sensor_Name) == 'MODIS-MAC') then
         Search_String = "MAC35S0"//trim(Search_String(9:nc))//"*"
      else
         Search_String = "MOD35"//trim(Search_String(9:nc))//"*"
      endif
    else
      if (trim(Sensor%Sensor_Name) == 'MODIS' .and. trim(Sensor%Platform_Name)  == 'AQUA') then
         Search_String = "MYDATML2"//trim(Search_String(9:nc))//"*"
      else
         Search_String = "MODATML2"//trim(Search_String(9:nc))//"*"
      endif
    endif

    Files => FILE_SEARCH(trim(Dir_Modis_Cloud_Mask),trim(Search_String),count=Num_Files)

    if (Num_Files == 0) then 
       print *, EXE_PROMPT, MODULE_PROMPT, "No MODIS Cloud Mask File Found"
       Auxiliary_Cloud_Mask_File_Name = "no_file"
       return
    endif

    if (Num_Files > 1) then
       print *, EXE_PROMPT, MODULE_PROMPT, "Multiple MODIS Cloud Mask Files Found"
       Auxiliary_Cloud_Mask_File_Name = "no_file"
    endif

    Auxiliary_Cloud_Mask_File_Name = Files(1)

    Files => null()
   
    ilen = len_trim(Auxiliary_Cloud_Mask_File_Name)
    if (trim(Sensor%Sensor_Name) == 'MODIS-MAC') then
      Auxiliary_Cloud_Mask_File_Name = Auxiliary_Cloud_Mask_File_Name(ilen-45:ilen)
    else  
      Auxiliary_Cloud_Mask_File_Name = Auxiliary_Cloud_Mask_File_Name(ilen-43:ilen)
    endif

end subroutine DETERMINE_MODIS_CLOUD_MASK_FILE

!----------------------------------------------------------------------
! Read MODIS data. 
!
! For reflectance channels, the output is scaled radiance
!
! For emissive channels, the output is brightness temperature
!
! nx, ny = size of input, output arrays
! nx_local, ny_local = size of arrays read in from hdf files
!----------------------------------------------------------------------
subroutine READ_MODIS_LEVEL1B(path,file_name,iband, &
                              calibrated_data_out, &
                              uncert_data_out, &
                              nx,ny,Seg_Idx,ny_total,ny_local_temp, &
                              Error_Status) 
      ! use CONSTANTS
      ! use HDF

      character(len=*), intent(in):: path
      character(len=*), intent(in):: file_name
      integer(kind=int4), intent(in):: iband
      integer(kind=int4), intent(in):: nx
      integer(kind=int4), intent(in):: ny
      integer(kind=int4), intent(in):: Seg_Idx
      integer(kind=int4), intent(in):: ny_total
      integer(kind=int4), intent(out):: ny_local_temp
      integer(kind=int4), intent(out):: Error_Status
      real(kind=int4), dimension(:,:), intent(out):: calibrated_data_out
      integer(kind=int1), dimension(:,:), intent(out):: uncert_data_out
      integer(kind=int4):: nx_local
      integer(kind=int4):: ny_local
      integer(kind=int4):: ny_start
      integer(kind=int4):: ny_end
      character(len=120):: sds_name
      character(len=120):: sds_name_uncert
      character(len=120):: sds_name_temp
      character(len=120):: scale_name
      character(len=120):: offset_name
      integer(kind=int4):: Sd_Id
      integer(kind=int4):: Sds_Id
      integer(kind=int4):: sds_uncert_id

      integer(kind=int4):: Status_Flag
      integer(kind=int4):: iband_sds
      integer(kind=int4):: nbands
      integer(kind=int4):: sfend, sfstart, sfselect, &
                           sfginfo, sfn2index
      integer(kind=int4):: num_attrs, sds_data_type, sds_rank
      integer(kind=int4), dimension(3):: sds_dims
      integer(kind=int4), dimension(:,:), allocatable:: i4_buffer
      integer(kind=int2), dimension(:,:,:), allocatable:: i2_buffer
      integer(kind=int2), dimension(:,:,:), allocatable:: i1_buffer
      real(kind=int4), dimension(:), allocatable:: scales
      real(kind=int4), dimension(:), allocatable:: offsets
      real(kind=int4), dimension(:,:), allocatable:: calibrated_data
      integer(kind=int4), dimension(3):: sds_start
      integer(kind=int4), dimension(3):: sds_stride
      integer(kind=int4), dimension(3):: sds_edges
      integer:: nx_min
      integer:: ny_min
      integer:: therm_flag
      integer:: iend
      integer(kind=int4):: sfrdata
      integer(kind=int4):: sfrnatt
      integer(kind=int4):: sfrcatt
      integer(kind=int4):: sffattr
      integer(kind=int4):: sfgainfo
      integer(kind=int4):: attr_index
      integer(kind=int4):: attr_data_type
      integer(kind=int4):: attr_count
      character(len=120):: attr_name
      character(len=1), dimension(:), allocatable:: band_names
      integer, dimension(:), allocatable:: band_int_names
      character(len=120):: band_name
      integer:: num_char_band_names
      integer::ii

      Error_Status = 0
      iend = 0
      Status_Flag = 0
      Sd_Id = 0
      therm_flag = sym%NO

error_check: do while (Status_Flag == 0 .and. iend == 0)

      !---- open file
      Sd_Id = sfstart(trim(path)//trim(file_name), DFACC_read)

      !-- if file is unreadable, exist this loop
      if (Sd_Id <= 0) then
         Status_Flag = 1
         exit
      endif

      !--- based on channel, set appropriate names and sds band index

      if (iband == 1 .or. iband == 2) then
        sds_name = "EV_250_Aggr1km_RefSB"
        sds_name_uncert = "EV_250_Aggr1km_RefSB_Uncert_Indexes"
        scale_name = "reflectance_scales"
        offset_name = "reflectance_offsets"
        band_name = "Band_250M"
      endif

      if (iband >= 3 .and. iband <= 7) then
        sds_name = "EV_500_Aggr1km_RefSB"
        sds_name_uncert = "EV_500_Aggr1km_RefSB_Uncert_Indexes"
        scale_name = "reflectance_scales"
        offset_name = "reflectance_offsets"
        band_name = "Band_500M"
      endif

      if ((iband >=8 .and. iband <= 19) .or. iband == 26) then
        sds_name = "EV_1KM_RefSB"
        sds_name_uncert = "EV_1KM_RefSB_Uncert_Indexes"
        scale_name = "reflectance_scales"
        offset_name = "reflectance_offsets"
        band_name = "Band_1km_RefSB"
      endif

      if (iband >=20 .and. iband <= 36 .and. iband /= 26) then
        sds_name = "EV_1KM_Emissive"
        sds_name_uncert = "EV_1KM_Emissive_Uncert_Indexes"
        scale_name = "radiance_scales"
        offset_name = "radiance_offsets"
        band_name = "Band_1KM_Emissive"
        therm_flag = sym%YES
      endif

      !--- open sds
      Sds_Id = sfselect(Sd_Id, sfn2index(Sd_Id,trim(sds_name)))
      sds_uncert_id = sfselect(Sd_Id, sfn2index(Sd_Id,trim(sds_name_uncert)))

      !--- get information
      Status_Flag = sfginfo(Sds_Id, sds_name_temp, sds_rank, sds_dims,  &
                        sds_data_type, num_attrs) + Status_Flag

!     print *, '================================'
!     print *, 'iband = ', iband
!     print *, 'Sds_Id = ', Sds_Id
!     print *, 'sds_name_temp = ', sds_name_temp
!     print *, 'sds_rank = ', sds_rank
!     print *, 'sds_dims = ', sds_dims
!     print *, 'num_attrs = ', num_attrs

      !--- define size of data
      nx_local = sds_dims(1)
      ny_local = sds_dims(2)
      nbands = sds_dims(3)

      !---  assign field band number from hdf file

      !--- convert these to a vector of integers
      attr_index = sffattr(Sds_Id,'band_names')
      Status_Flag = sfgainfo(Sds_Id,attr_index,attr_name, attr_data_type, attr_count)
      num_char_band_names = attr_count

      allocate(band_names(num_char_band_names))
      allocate(band_int_names(nbands))

      Status_Flag = sfrcatt(Sds_Id,attr_index,band_names)
      where(band_names == "o" .or. band_names =="h" .or.  &
            band_names == "l" .or.  band_names=="i" .or. &
            ichar(band_names) == 0)
              band_names = " "
      endwhere

      band_name = ""
      do ii = 1,num_char_band_names
            band_name = trim(band_name) // band_names(ii)
      enddo
    
      read(band_name,fmt=*) band_int_names

      !--- look for the band in this list of bands
      iband_sds = -1
      DO ii = 1 , nbands
         if (band_int_names(ii) == iband) iband_sds = ii 
      END DO
      
      if (iband_sds == -1) then
              print *, EXE_PROMPT, "Modis Band not Available "
              Status_Flag = 1
      endif

      !--- compute number of lines to read for this segment
      ny_start = (Seg_Idx-1)*ny + 1
      ny_end = min(ny_start+ny-1,ny_total)
      ny_local_temp = ny_end - ny_start + 1

      nx_min = min(nx,nx_local)
      ny_min = min(ny,ny_local_temp)

      allocate(i4_buffer(nx_local,ny_local_temp))
      allocate(i2_buffer(nx_local,ny_local_temp,1))
      allocate(i1_buffer(nx_local,ny_local_temp,1))

      allocate(scales(nbands))
      allocate(offsets(nbands))
      allocate(calibrated_data(nx_local,ny_local_temp))

      !-- read channel data for specified channel
      sds_stride = (/1, 1, 1/)
      sds_start = (/0, ny_start-1,iband_sds-1/)
      sds_edges = (/nx_local, ny_local_temp, 1/)

      !-- read scaled data
      Status_Flag = sfrdata(Sds_Id, sds_start, sds_stride, sds_edges, i2_buffer) + Status_Flag

      !-- read uncertainty indexes
      Status_Flag = sfrdata(sds_uncert_id, sds_start, sds_stride, sds_edges, i1_buffer) + Status_Flag

      !-- read calibration information
      Status_Flag = sfrnatt(Sds_Id, sffattr(Sds_Id,scale_name), scales) + Status_Flag
      Status_Flag = sfrnatt(Sds_Id, sffattr(Sds_Id,offset_name), offsets) + Status_Flag

     !Status_Flag =  HDF_SDS_readER(Sd_Id, sds_name, sds_start, sds_stride, sds_edges, i2_buffer) + Status_Flag
     !Status_Flag =  HDF_SDS_ATTRIBUTE_readER(Sd_Id, sds_name,scale_name,scales ) + Status_Flag
     !Status_Flag =  HDF_SDS_ATTRIBUTE_readER(Sd_Id, sds_name,offset_name,offsets ) + Status_Flag

      !--- close file
      Status_Flag = sfend(Sd_Id)

      !--- convert from i2 to i4 because it's unsigned integer (Denis B.)
      i4_buffer=i2_buffer(:,:,1)
      where (i4_buffer < 0 .and. i4_buffer .ne. fill_value) i4_buffer=i4_buffer+2**16

      !----- calibrate counts
      !--radiances have units of  Watts/m^2/micrometer/steradian
      !--reflectances range from 0.0 to 1.0

      calibrated_data = scales(iband_sds)*(i4_buffer - offsets(iband_sds))

      !--- scale reflectaces to 0 to 100%
      if (therm_flag == sym%NO) then
              calibrated_data = 100.0*calibrated_data
                   where(calibrated_data < 0.0)
                          calibrated_data = missing_value
                   endwhere
      endif

      where(i4_buffer == fill_value) 
              calibrated_data = missing_value
      endwhere

      ! --- specific of aqua 1.6um fix (Denis B. - 06.05.2015)
      if (iband == 6) then
         where(i4_buffer == 65528)
            calibrated_data = missing_value                                                                                                    
         endwhere
      endif

      calibrated_data_out = missing_value
      calibrated_data_out(1:nx_min,1:ny_min) = calibrated_data(1:nx_min,1:ny_min)

      uncert_data_out = missing_value_int1
      uncert_data_out(1:nx_min,1:ny_min) = i1_buffer(1:nx_min,1:ny_min,1)

      iend = 1

      enddo error_check  !end of while loop

      !--- deallocate memory
      if (allocated(band_names)) deallocate(band_names)
      if (allocated(band_int_names)) deallocate(band_int_names)
      if (allocated(i4_buffer)) deallocate(i4_buffer)
      if (allocated(i2_buffer)) deallocate(i2_buffer)
      if (allocated(i1_buffer)) deallocate(i1_buffer)
      if (allocated(scales)) deallocate(scales)
      if (allocated(offsets)) deallocate(offsets)
      if (allocated(calibrated_data)) deallocate(calibrated_data)

      if (Status_Flag /= 0) then

              !--- set output error status
              Error_Status = 1

              print *, EXE_PROMPT, "Reading of MODIS Level1b failed, skipping this file"

              !--- close hdf file, if loop exiting after opening it
              if (Sd_Id > 0)  then 
                      Status_Flag = sfend(Sd_Id)
              endif

      endif

end subroutine READ_MODIS_LEVEL1B

!----------------------------------------------------------------------
! Read the geolocation information
!
! nx = number of elements (x-dir)
! ny = numer of lines (y-dir)
! Seg_Idx = the segment number (starts with 1)
! ny_total = total number of lines in file
! nx_local = actual number of elements in this file
! ny_local = actual number of lines in this file
! ny_local_temp = actual number of lines from this file for this segment
! time_ms_start = millisecond time at start of file (line=1)
! time_ms_end = millisecond time at end of file (line=Image%Number_Of_Lines)
! scan_number = number of line within this file
! asc_des_flag = 0 for ascending, 1 for descednding
!----------------------------------------------------------------------
subroutine READ_MODIS_LEVEL1B_GEOLOCATION(path,file_name,  &
                              lon_out,  &
                              lat_out, &
                              sensor_zenith_out,  &
                              sensor_azimuth_out, &
                              solar_zenith_out, &
                              solar_azimuth_out, &
                              relative_azimuth_out, &
                              nx,ny,Seg_Idx,ny_total,ny_local_temp, &
                              time_ms_start, time_ms_end, time_ms, &
                              scan_number, asc_des_flag, Error_Status) 

      character(len=*), intent(in):: path
      character(len=*), intent(in):: file_name
      integer(kind=int4), intent(in):: nx
      integer(kind=int4), intent(in):: ny
      integer(kind=int4), intent(in):: Seg_Idx
      integer(kind=int4), intent(in):: ny_total
      integer(kind=int4), intent(in):: time_ms_start
      integer(kind=int4), intent(in):: time_ms_end
      integer(kind=int4), intent(out), dimension(:):: time_ms
      integer(kind=int4), intent(out), dimension(:):: scan_number
      integer(kind=int1), intent(out), dimension(:):: asc_des_flag
      integer(kind=int4), intent(out):: ny_local_temp
      integer(kind=int4), intent(out):: Error_Status
      real(kind=real4), dimension(:,:), intent(out):: lon_out
      real(kind=real4), dimension(:,:), intent(out):: lat_out
      real(kind=real4), dimension(:,:), intent(out):: sensor_zenith_out
      real(kind=real4), dimension(:,:), intent(out):: sensor_azimuth_out
      real(kind=real4), dimension(:,:), intent(out):: solar_zenith_out
      real(kind=real4), dimension(:,:), intent(out):: solar_azimuth_out
      real(kind=real4), dimension(:,:), intent(out):: relative_azimuth_out

      real(kind=real4), allocatable, dimension(:,:):: r4_buffer
      integer(kind=int2), allocatable, dimension(:,:):: i2_buffer
      character(len=120):: sds_name
      integer(kind=int4), dimension(2):: sds_dims
      integer(kind=int4), dimension(2):: sds_start
      integer(kind=int4), dimension(2):: sds_stride
      integer(kind=int4), dimension(2):: sds_edges
      integer(kind=int4):: nx_local
      integer(kind=int4):: ny_local
      integer(kind=int4):: ny_start
      integer(kind=int4):: ny_end
      integer:: nx_min
      integer:: ny_min
      integer:: Status_Flag
      integer(kind=int4):: Sd_Id
      integer(kind=int4):: Sds_Id
      integer(kind=int4):: num_attrs, sds_data_type, sds_rank

      integer(kind=int4):: sfend, sfstart, sfselect, sfrdata, sfendacc, &
                           sfginfo, sfn2index
      integer(kind=int4):: iline
      real(kind=real4):: dtime_dline
      integer(kind=int4):: iend

      integer(kind=int4), parameter:: time_last_gran_start = 86100000
      integer(kind=int4), parameter:: time_last_gran_end   = 86400000 

      Status_Flag = 0
      Error_Status = 0
      iend = 0

error_check: do while (Status_Flag == 0 .and. iend == 0)

      !---- open file
      Sd_Id = sfstart(trim(path)//trim(file_name), DFACC_read)

      !-- if file is unreadable, exist this loop
      if (Sd_Id <= 0) then
         Status_Flag = 1
         exit
      endif
 
      !--- Open Latitude and extract data sizes
      Sds_Id = sfselect(Sd_Id, sfn2index(Sd_Id,trim('Latitude')))
      !-- if latitude is missing, exit this loop  (note,the do while is working)
      if (Sds_Id <= 0) then
         Status_Flag = Status_Flag + 1
         exit
      endif

      Status_Flag = sfginfo(Sds_Id, sds_name, sds_rank, sds_dims, sds_data_type, num_attrs) + Status_Flag
      Status_Flag = sfendacc(Sds_Id) + Status_Flag
      nx_local = sds_dims(1)
      ny_local = sds_dims(2)
      sds_stride = (/1, 1/)
 
      !--- compute number of lines to read for this segment
      ny_start = (Seg_Idx-1)*ny + 1
      ny_end = min(ny_start+ny-1,ny_total)
      ny_local_temp = ny_end - ny_start + 1

      sds_start = (/0, ny_start-1/)
      sds_edges = (/nx_local, ny_local_temp/)

      nx_min = min(nx,nx_local)
      ny_min = min(ny,ny_local_temp)

      !--- allocate space for data to be read in
      allocate(r4_buffer(nx_local,ny_local_temp))
      allocate(i2_buffer(nx_local,ny_local_temp))

      !--- Latitude
      Sds_Id = sfselect(Sd_Id, sfn2index(Sd_Id,trim('Latitude')))
      Status_Flag = sfrdata(Sds_Id, sds_start, sds_stride, sds_edges, r4_buffer) + Status_Flag
      lat_out(1:nx_min,1:ny_min) = r4_buffer(1:nx_min,1:ny_min)
      Status_Flag = sfendacc(Sds_Id) + Status_Flag
 
      !--Longitude
      Sds_Id = sfselect(Sd_Id, sfn2index(Sd_Id,trim('Longitude')))
      Status_Flag = sfrdata(Sds_Id, sds_start, sds_stride, sds_edges, r4_buffer) + Status_Flag
      lon_out(1:nx_min,1:ny_min) = r4_buffer(1:nx_min,1:ny_min)
      Status_Flag = sfendacc(Sds_Id) + Status_Flag

      !--Sensor Zenith
      Sds_Id = sfselect(Sd_Id, sfn2index(Sd_Id,trim('SensorZenith')))
      Status_Flag = sfrdata(Sds_Id, sds_start, sds_stride, sds_edges, i2_buffer) + Status_Flag
      sensor_zenith_out(1:nx_min,1:ny_min) = i2_buffer(1:nx_min,1:ny_min)/100.0
      Status_Flag = sfendacc(Sds_Id) + Status_Flag
      
      !--Solar Zenith
      Sds_Id = sfselect(Sd_Id, sfn2index(Sd_Id,trim('SolarZenith')))
      Status_Flag = sfrdata(Sds_Id, sds_start, sds_stride, sds_edges, i2_buffer) + Status_Flag
      solar_zenith_out(1:nx_min,1:ny_min) = i2_buffer(1:nx_min,1:ny_min)/100.0
      Status_Flag = sfendacc(Sds_Id) + Status_Flag

      !--Sensor Azimuth 
      Sds_Id = sfselect(Sd_Id, sfn2index(Sd_Id,trim('SensorAzimuth')))
      Status_Flag = sfrdata(Sds_Id, sds_start, sds_stride, sds_edges, i2_buffer) + Status_Flag
      sensor_azimuth_out(1:nx_min,1:ny_min) = i2_buffer(1:nx_min,1:ny_min)/100.0
      Status_Flag = sfendacc(Sds_Id) + Status_Flag

      !--Solar Azimuth
      Sds_Id = sfselect(Sd_Id, sfn2index(Sd_Id,trim('SolarAzimuth')))
      Status_Flag = sfrdata(Sds_Id, sds_start, sds_stride, sds_edges, i2_buffer) + Status_Flag
      solar_azimuth_out(1:nx_min,1:ny_min) = i2_buffer(1:nx_min,1:ny_min)/100.0
      Status_Flag = sfendacc(Sds_Id) + Status_Flag
 
      !--- close file
      Status_Flag = sfend(Sd_Id) + Status_Flag

      !--- make relative azimuth
      relative_azimuth_out = abs(solar_azimuth_out - sensor_azimuth_out)
      where (relative_azimuth_out > 180.0)
                relative_azimuth_out = 360 - relative_azimuth_out
      end where
      relative_azimuth_out = 180 - relative_azimuth_out

      !---- compute additional angles

      !--- compute the cosine of the glint zenith angle
      Geo%Glintzen = glint_angle ( Geo%Solzen , Geo%Satzen , Geo%Relaz  )

      !--- compute the cosine of the scattering angle
      Geo%Scatangle = scattering_angle( Geo%Solzen , Geo%Satzen , Geo%Relaz )

      !--- redefine solar azimuth to be consistent with avhrr
      where (Geo%Solaz /= Missing_Value_Real4)
       Geo%Solaz = 180.0 - abs(Geo%Solaz)
      end where

      !---- compute time
      if (time_ms_end .eq. 0 .and. time_ms_start .eq. time_last_gran_start) then
      ! end time of last granule is adjusted from 0 so time increment is positive
         dtime_dline = (time_ms_end + time_last_gran_end - time_ms_start) / (Image%Number_Of_Lines-1) 
      else
         dtime_dline = (time_ms_end - time_ms_start) / (Image%Number_Of_Lines-1)
      endif

      do iline = 1,ny_local_temp
        scan_number(iline) = ny_start + iline - 1
        time_ms(iline) = int(float(iline - 1 + ny_start-1)*dtime_dline) + time_ms_start 
      enddo

      !--- compute ascending descending flag
      asc_des_flag = 0
      do iline = 1,ny_local_temp-1
          if (lat_out(nx_local/2,iline+1) > lat_out(nx_local/2,iline)) then
                  asc_des_flag(iline) = 0
          else
                  asc_des_flag(iline) = 1
          endif
      enddo
      asc_des_flag(ny_local_temp) = asc_des_flag(ny_local_temp-1)

      iend = 1

      enddo  error_check ! end of while loop

      !--- deallocate memory
      !--- clean up memory
      if (allocated(r4_buffer)) deallocate(r4_buffer)
      if (allocated(i2_buffer)) deallocate(i2_buffer)

      if (Status_Flag /= 0) then

              !--- set output error status
              Error_Status = 1

              print *, EXE_PROMPT, "Reading of MODIS Geolocation failed, skipping this file"

              !--- close hdf file, if loop exiting after opening it
              if (Sd_Id > 0)  then
                      Status_Flag = sfend(Sd_Id)
              endif

      endif

end subroutine READ_MODIS_LEVEL1B_GEOLOCATION

!----------------------------------------------------------------------
! Read the cloud mask
!
! nx = number of elements (x-dir)
! ny = numer of lines (y-dir)
! Seg_Idx = the segment number (starts with 1)
! ny_total = total number of lines in file
! nx_local = actual number of elements in this file
! ny_local = actual number of lines in this file
! ny_local_temp = actual number of lines from this file for this segment
! scan_number = number of line within this file
!----------------------------------------------------------------------
subroutine READ_MODIS_LEVEL1B_CLOUD_MASK(path,file_name,  &
                              Cloud_Mask_Out,  &
                              Cloud_Phase_Out,  &
                              Cloud_Type_Out,  &
                              Cloud_Height_Out,  &
                              Cloud_Opd_Out, &
                              Cloud_Reff_Out, &
                              nx,ny,Seg_Idx,ny_total,ny_local_temp, &
                              Error_Status) 

      character(len=*), intent(in):: path
      character(len=*), intent(in):: file_name
      integer(kind=int4), intent(in):: nx
      integer(kind=int4), intent(in):: ny
      integer(kind=int4), intent(in):: Seg_Idx
      integer(kind=int4), intent(in):: ny_total
      integer(kind=int4), intent(out):: ny_local_temp
      integer(kind=int4), intent(out):: Error_Status
      integer(kind=int1), dimension(:,:), intent(out):: Cloud_Mask_Out
      integer(kind=int1), dimension(:,:), intent(out):: Cloud_Phase_Out
      integer(kind=int1), dimension(:,:), intent(out):: Cloud_Type_Out
      real(kind=real4), dimension(:,:), intent(out):: Cloud_Height_Out
      real(kind=real4), dimension(:,:), intent(out):: Cloud_Opd_Out
      real(kind=real4), dimension(:,:), intent(out):: Cloud_Reff_Out

      integer(kind=int1), allocatable, dimension(:,:,:):: i1_buffer
      integer(kind=int2), allocatable, dimension(:,:):: i2_buffer
      character(len=120):: sds_name
      integer(kind=int4), dimension(3):: sds_dims
      integer(kind=int4), dimension(3):: sds_start
      integer(kind=int4), dimension(3):: sds_stride
      integer(kind=int4), dimension(3):: sds_edges
      integer(kind=int4):: nx_local
      integer(kind=int4):: ny_local
      integer(kind=int4):: ny_start
      integer(kind=int4):: ny_end
      integer:: nx_min
      integer:: ny_min
      integer:: Status_Flag
      integer(kind=int4):: Sd_Id
      integer(kind=int4):: Sds_Id
      integer(kind=int4):: num_attrs, sds_data_type, sds_rank

      integer(kind=int4):: sfend, sfstart, sfselect, sfrdata, sfendacc, &
                           sfginfo, sfn2index
      integer(kind=int4):: iend
      logical:: Atml2_File



      Status_Flag = 0
      Error_Status = 0
      iend = 0

      Atml2_File = .false.
      if (index(file_name, "ATML2") > 0) Atml2_File = .true.

error_check: do while (Status_Flag == 0 .and. iend == 0)

      !---- open file
      Sd_Id = sfstart(trim(path)//trim(file_name), DFACC_read)

      !-- if file is unreadable, exist this loop
      if (Sd_Id <= 0) then
         Status_Flag = 1
         exit
      endif
 
      !--- Open cloud mask and extract data sizes
      Sds_Id = sfselect(Sd_Id, sfn2index(Sd_Id,trim('Cloud_Mask')))
      if (Sds_Id <=0) Status_Flag = Status_Flag + 1
      Status_Flag = sfginfo(Sds_Id, sds_name, sds_rank, sds_dims, sds_data_type, num_attrs) + Status_Flag
      Status_Flag = sfendacc(Sds_Id) + Status_Flag
      nx_local = sds_dims(1)
      ny_local = sds_dims(2)
      sds_stride = (/1, 1, 1/)

      !--- compute number of lines to read for this segment
      ny_start = (Seg_Idx-1)*ny + 1
      ny_end = min(ny_start+ny-1,ny_total)
      ny_local_temp = ny_end - ny_start + 1

      nx_min = min(nx,nx_local)
      ny_min = min(ny,ny_local_temp)

      !--- allocate space for data to be read in
      if (Sensor%Spatial_Resolution_Meters == 1000) then
        allocate(i1_buffer(nx_local,ny_local_temp,6))
        allocate(i2_buffer(nx_local,ny_local_temp))
        sds_start = (/0, ny_start-1, 0/)
        sds_edges = (/nx_local, ny_local_temp,6/)
      else
        allocate(i1_buffer(nx_local,ny_local_temp,1))
        allocate(i2_buffer(nx_local,ny_local_temp))
        sds_start = (/0, ny_start-1, 0/)
        sds_edges = (/nx_local, ny_local_temp,1/)
      endif

      !--- Read Cloud Mask Bytes (packed)
      Sds_Id = sfselect(Sd_Id, sfn2index(Sd_Id,trim('Cloud_Mask')))
      Status_Flag = sfrdata(Sds_Id, sds_start, sds_stride, sds_edges, i1_buffer) + Status_Flag
      Cloud_Mask_Out(1:nx_min,1:ny_min) = i1_buffer(1:nx_min,1:ny_min,1)
      Status_Flag = sfendacc(Sds_Id) + Status_Flag

      if (Atml2_File) then
       !--- Read IR Cloud Phase
       !--- 0 = cloud free, 1 = water, 2 = ice, 3 = mixed, 6 = unknown)
       Sds_Id = sfselect(Sd_Id, sfn2index(Sd_Id,trim('Cloud_Phase_Infrared_1km')))
       Status_Flag = sfrdata(Sds_Id, sds_start, sds_stride, sds_edges, i1_buffer) + Status_Flag
       Cloud_Phase_Out(1:nx_min,1:ny_min) = i1_buffer(1:nx_min,1:ny_min,1)
       Status_Flag = sfendacc(Sds_Id) + Status_Flag

       !--- Read Cloud Top Height (i2, scale = 1, add_offset = 0, fill=-999)
       Sds_Id = sfselect(Sd_Id, sfn2index(Sd_Id,trim('Cloud_Top_Height_1km')))
       Status_Flag = sfrdata(Sds_Id, sds_start, sds_stride, sds_edges, i2_buffer) + Status_Flag
       Cloud_Height_Out(1:nx_min,1:ny_min) = real(i2_buffer(1:nx_min,1:ny_min))
       where(i2_buffer == -999)
          Cloud_Height_Out = Missing_Value_Real4
       endwhere
       Status_Flag = sfendacc(Sds_Id) + Status_Flag

       !--- Read Cloud Optical Depth (i2, scale = 0.01, add_offset = 0, fill=-9999)
       Sds_Id = sfselect(Sd_Id, sfn2index(Sd_Id,trim('Cloud_Optical_Thickness')))
       Status_Flag = sfrdata(Sds_Id, sds_start, sds_stride, sds_edges, i2_buffer) + Status_Flag
       Cloud_Opd_Out(1:nx_min,1:ny_min) = 0.0 + 0.01*real(i2_buffer(1:nx_min,1:ny_min))
       where(i2_buffer == -9999)
          Cloud_Opd_Out = Missing_Value_Real4
       endwhere
       Status_Flag = sfendacc(Sds_Id) + Status_Flag

       !--- Read Cloud Particle Size( i2, scale = 0.01, add_offset = 0, fill=-9999)
       Sds_Id = sfselect(Sd_Id, sfn2index(Sd_Id,trim('Cloud_Effective_Radius')))
       Status_Flag = sfrdata(Sds_Id, sds_start, sds_stride, sds_edges, i2_buffer) + Status_Flag
       Cloud_Reff_Out(1:nx_min,1:ny_min) = 0.0 + 0.01*real(i2_buffer(1:nx_min,1:ny_min))
       where(i2_buffer == -9999)
          Cloud_Reff_Out = Missing_Value_Real4
       endwhere
       Status_Flag = sfendacc(Sds_Id) + Status_Flag

      endif

      !--- close file
      Status_Flag = sfend(Sd_Id) + Status_Flag

      iend = 1

      enddo  error_check ! end of while loop

      !--- unpacked needed information for 4-level cloud mask
      Cloud_Mask_Out = ishft(ishft(Cloud_Mask_Out,5),-6)

      !--- switch CLAVR-x convection for mask
      Cloud_Mask_Out = 3-Cloud_Mask_Out

      !--- switch CLAVR-x convection for phase
      where(Cloud_Phase_Out == 0)
          Cloud_Type_Out = sym%CLEAR_TYPE
      endwhere
      where(Cloud_Phase_Out == 1)
          Cloud_Type_Out = sym%WATER_TYPE
      endwhere
      where(Cloud_Phase_Out == 2)
          Cloud_Type_Out = sym%OPAQUE_ICE_TYPE
      endwhere
      where(Cloud_Phase_Out == 3)
          Cloud_Type_Out = sym%MIXED_TYPE
      endwhere
      where(Cloud_Phase_Out == 6)
          Cloud_Type_Out = sym%UNKNOWN_TYPE
      endwhere

      !--- modify phase (0,1,3 have same meaning as CLAVR-x)
      where(Cloud_Phase_Out == 2)
          Cloud_Phase_Out = sym%ICE_PHASE
      endwhere
      where(Cloud_Phase_Out == 5)
          Cloud_PHASE_Out = sym%UNKNOWN_PHASE
      endwhere


      !--- deallocate memory
      !--- clean up memory
      if (allocated(i1_buffer)) deallocate(i1_buffer)
      if (allocated(i2_buffer)) deallocate(i2_buffer)

      if (Status_Flag /= 0) then

              !--- set output error status
              Error_Status = 1

              print *, EXE_PROMPT, "Reading of MODIS Cloud Mask failed, skipping this file"

              !--- close hdf file, if loop exiting after opening it
              if (Sd_Id > 0)  then
                      Status_Flag = sfend(Sd_Id)
              endif

      endif

end subroutine READ_MODIS_LEVEL1B_CLOUD_MASK


!======================================================================
subroutine READ_MODIS(Seg_Idx,Error_Status)

    integer, intent(in):: Seg_Idx
    integer, intent(out):: Error_Status

    integer:: End_Flag
    integer:: Chan_Idx
    real, dimension(20:36):: Eff_Wavenumber
    character (len=1020) :: File_Name_Tmp

    Error_Status = 0
    End_Flag = 0

    Eff_Wavenumber = (/Planck_Nu(20), Planck_Nu(21), Planck_Nu(22), Planck_Nu(23),  &
                       Planck_Nu(24), Planck_Nu(25), Missing_Value_Real4,  &
                       Planck_Nu(27), Planck_Nu(28), Planck_Nu(29), Planck_Nu(30),  &
                       Planck_Nu(31), Planck_Nu(32), Planck_Nu(33), Planck_Nu(34),  &
                       Planck_Nu(35), Planck_Nu(36)/)

error_check: do while (Error_Status == 0 .and. End_Flag == 0)

    if (trim(Sensor%Sensor_Name) == 'MODIS' .or. trim(Sensor%Sensor_Name) == 'MODIS-CSPP') then
  
       !-- read geolocation (and compute scan time)
       if (Sensor%Spatial_Resolution_Meters == 1000) then
          File_Name_Tmp = trim(Image%Auxiliary_Geolocation_File_Name)
       else
          File_Name_Tmp = trim(Image%Level1b_Name)
       endif

       call READ_MODIS_LEVEL1B_GEOLOCATION(trim(Image%Level1b_Path),  &
                                           trim(File_Name_Tmp), &
                                           Nav%Lon_1b, Nav%Lat_1b, & 
                                           Geo%Satzen, Geo%Sataz, Geo%Solzen, Geo%Solaz, Geo%Relaz, &
                                           Image%Number_Of_Elements,Image%Number_Of_Lines_Per_Segment, &
                                           Seg_Idx,Image%Number_Of_Lines,Image%Number_Of_Lines_Read_This_Segment, &
                                           Image%start_Time,Image%End_Time,scan_time, &
                                           Scan_Number,Nav%Ascend,Error_Status)


       if (Error_Status /= 0) exit

       !--- read cloud mask - assume path is same as level-1b
       Cloud_Mask_Aux_Read_Flag = sym%NO
 
       if (Cloud_Mask_Aux_Flag /= sym%NO_AUX_CLOUD_MASK) then

          call READ_MODIS_LEVEL1B_CLOUD_MASK(trim(Image%Level1b_Path),  &
                                             trim(Image%Auxiliary_Cloud_Mask_File_Name), &
                                             Cld_Mask_Aux, & 
                                             Cld_Phase_Aux, & 
                                             Cld_Type_Aux, & 
                                             Zc_Aux, & 
                                             Tau_Aux, &
                                             Reff_Aux, &
                                             Image%Number_Of_Elements,Image%Number_Of_Lines_Per_Segment, &
                                             Seg_Idx,Image%Number_Of_Lines,Image%Number_Of_Lines_Read_This_Segment, &
                                             Error_Status)


          if (Error_Status == 0) then
             Cloud_Mask_Aux_Read_Flag = sym%YES
          endif

       endif

       !--- read channel data
       do Chan_Idx = 1,36

         if (Sensor%Chan_On_Flag_Default (Chan_Idx) == sym%NO) cycle

         if (Chan_Idx < 20 .or. Chan_Idx == 26) then
            call READ_MODIS_LEVEL1B(trim(Image%Level1b_Path),trim(Image%Level1b_Name), &
                               Chan_Idx, &
                               ch(Chan_Idx)%Ref_Toa,  &
                               ch(Chan_Idx)%Unc, &
                               Image%Number_Of_Elements,Image%Number_Of_Lines_Per_Segment, &
                               Seg_Idx,Image%Number_Of_Lines,Image%Number_Of_Lines_Read_This_Segment,Error_Status)
         endif

         if (Chan_Idx >= 20 .and. Chan_Idx /= 26) then
           call READ_MODIS_THERMAL_BAND(Chan_Idx, &
                                        Seg_Idx, &
                                        Eff_Wavenumber(Chan_Idx), &
                                        Image%Number_Of_Elements, &
                                        Image%Number_Of_Lines_Per_Segment, &
                                        trim(Image%Level1b_Path), &
                                        trim(Image%Level1b_Name), &
                                        ch(Chan_Idx)%Rad_Toa, &
                                        ch(Chan_Idx)%Bt_Toa, &
                                        ch(Chan_Idx)%Unc,Error_Status)
         endif

       enddo

       !--- Quality check MODIS data
       call QC_MODIS(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)

    endif

    End_Flag = 1

    enddo error_check

    if (Error_Status /= 0) then
            print *, EXE_PROMPT, "Error reading MODIS data, skipping file"
    endif

    end subroutine READ_MODIS

!======================================================================
! READ_MODIS_TIME_ATTR
! Gets the start time/date and end time/date from metadata
!======================================================================

    subroutine READ_MODIS_TIME_ATTR(Path, File_Name, Start_Year, Start_Day, &
                                    Start_Time, End_Year, End_Day, End_Time)

      character(len=*), intent(in):: Path
      character(len=*), intent(in):: File_Name
      integer(kind=int2), intent(out):: Start_Year
      integer(kind=int2), intent(out):: Start_Day
      integer(kind=int2), intent(out):: End_Year
      integer(kind=int2), intent(out):: End_Day      
      integer(kind=int4), intent(out):: Start_Time
      integer(kind=int4), intent(out):: End_Time 

      integer(kind=int4):: Status_Flag
      integer(kind=int4):: Sd_Id
      integer(kind=int4):: Sds_Id
      integer(kind=int4):: sfend
      integer(kind=int4):: sfstart
      integer(kind=int4):: sfrcatt
      integer(kind=int4):: sffattr
      integer(kind=int4):: sfgainfo

      integer(kind=int4), parameter:: string_length_assumed = 20000
      !string to get metadata
      character(len=string_length_assumed):: Metadata_Temp

      character(len=200):: attr_name
      integer:: data_type
      integer:: count
      
      !char arrays from mega array
      character(len=15):: start_time_eos
      character(len=15):: end_time_eos
      character(len=10):: start_date
      character(len=10):: end_date
      
      !string index in mega array
      integer(kind=int4) :: String_Index

      !integers to get info from subarrays
      integer(kind=int4):: hour
      integer(kind=int4):: minute
      integer(kind=int4):: seconds
 
      integer(kind=int4):: day
      integer(kind=int4):: month
      integer(kind=int4):: year
      integer(kind=int4):: jday

      Status_Flag = 0
     
      !---- open file
      Sd_Id = sfstart(trim(Path)//trim(File_Name), DFACC_read)

      !--- determine attribute index
      Sds_Id = sffattr(Sd_Id, "CoreMetadata.0")

      !--- get information on attribute
      Status_Flag = sfgainfo(Sd_Id, Sds_Id, attr_name, data_type, count)

      !--- check if attribute's size exceeds expectation, if so warn
      if (count > string_length_assumed) then
        print *, EXE_PROMPT, MODULE_PROMPT, " Warning: inconsistent attribute length in READ_MODIS_TIME_ATTR"
      endif
      
      !--- read attribute
      Status_Flag = sfrcatt(Sd_Id, Sds_Id, Metadata_Temp) + Status_Flag  

      !--- close file
      Status_Flag = sfend(Sd_Id) + Status_Flag
           
      !Lets do the begining time now. First we get the index in the string
      !where the string appears
      String_Index=index(Metadata_Temp, 'RANGEBEGINNINGTIME')
      start_time_eos = Metadata_Temp(String_Index+80:String_Index+80+14)

      !parse and make start/end time in msec
      
      read(start_time_eos(1:2), fmt="(I2)") hour     
      read(start_time_eos(4:5), fmt="(I2)") minute
      read(start_time_eos(7:8), fmt="(I2)") seconds
      
      Start_Time =  (hour + (float(minute)/60.0) + ((float(seconds)/60.0/60.0))) * &
      60.0 * 60.0 * 1000.0
      
      !now do the end time
      
      String_Index = index(Metadata_Temp, 'RANGEENDINGTIME')
      end_time_eos = Metadata_Temp(String_Index+77:String_Index+77+14)

      read(end_time_eos(1:2), fmt="(I2)") hour     
      read(end_time_eos(4:5), fmt="(I2)") minute
      read(end_time_eos(7:8), fmt="(I2)") seconds
      
      End_Time =  (hour + (float(minute)/60.0) + ((float(seconds)/60.0/60.0))) * &
      60.0 * 60.0 * 1000.0
             
      !now we have to get the date (which is in yyyy-mm-dd)
      String_Index=index(Metadata_Temp, 'RANGEBEGINNINGDATE')
      start_date = Metadata_Temp(String_Index+80:String_Index+80+9)
     
      read(start_date(1:4), fmt="(I4)") year     
      read(start_date(6:7), fmt="(I2)") month
      read(start_date(9:10), fmt="(I2)") day
      
      CALL JULIAN(day, month, year, jday)
     
      Start_Year= year
      Start_Day = jday

      !Now we do it for the end day
      String_Index = index(Metadata_Temp, 'RANGEENDINGDATE')
      end_date = Metadata_Temp(String_Index+77:String_Index+77+9)

      read(end_date(1:4), fmt="(I4)") End_Year     
      read(end_date(6:7), fmt="(I2)") month
      read(end_date(9:10), fmt="(I2)") day
      
      CALL JULIAN(day, month, year, jday)
      
      End_Year = year
      End_Day = jday
    
    end subroutine READ_MODIS_TIME_ATTR
    !======================================================================
    ! READ_MODIS_SIZE_ATTR
    ! Gets the size of the array from metadata
    !======================================================================

    subroutine READ_MODIS_SIZE_ATTR(File_Name_Full, Num_Elements, Num_Lines)

      character(len=*), intent(in):: File_Name_Full
      integer(kind=int4), intent(out):: Num_Elements
      integer(kind=int4), intent(out):: Num_Lines

      integer(kind=int4):: Status_Flag
      integer(kind=int4):: Sd_Id
      integer(kind=int4):: Sds_Id
      integer(kind=int4):: sfend
      integer(kind=int4):: sfstart
      integer(kind=int4):: sfrcatt
      integer(kind=int4):: sffattr
      integer(kind=int4):: sfgainfo

      integer(kind=int4), parameter:: string_length_assumed = 32000
      
      !Mega character array to get metadata
      character(len=string_length_assumed):: Metadata_Temp
      character(len=200):: attr_name
      integer:: data_type, count

      character(len=15):: Temp_Char
      character(len=15):: Temp_5km_Char
      integer(kind=int4) :: String_Index
      integer(kind=int4) :: Num_Char, i
      logical :: File_5km_Test

      Status_Flag = 0
      File_5km_Test = .FALSE.
      
      !---- open file
      Sd_Id = sfstart(trim(File_Name_Full), DFACC_read)

      !--- determine attribute index
      Sds_Id = sffattr(Sd_Id, "StructMetadata.0")

      !--- get information on attribute
      Status_Flag = sfgainfo(Sd_Id, Sds_Id, attr_name, data_type, count)

      !--- check if assumed size is wrong, if so warn
      if (count > string_length_assumed) then
        print *, EXE_PROMPT, MODULE_PROMPT, " Warning: inconsistent attribute length in READ_MODIS_SIZE_ATTR"
      endif
      
      !--- read attribute
      Status_Flag = sfrcatt(Sd_Id, Sds_Id, Metadata_Temp) + Status_Flag  

      !--- close file
      Status_Flag = sfend(Sd_Id) + Status_Flag

      !--- In order to test for the 5km file, we know that the 2*nscans will 
      !    be different than 10*scans 
      String_Index=index(Metadata_Temp, '2*nscans')
      Temp_5km_Char = trim(Metadata_Temp(String_Index+19:String_Index+19+5))
           
      !--- the total number of scans
      String_Index=index(Metadata_Temp, '10*nscans')
      Temp_Char = trim(Metadata_Temp(String_Index+20:String_Index+20+5))
     
      
      !Test to see if '10*nscans' .eqv. '2*nscans'. This is only true
      ! for 5km MODIS files. This test does NOT key off the file name
      ! just in case the file name is different than 02SSH.
      if (trim(Temp_Char) .EQ. trim(Temp_5km_Char)) THEN
           File_5km_Test = .TRUE. 
      endif
      
      !usage of the I3 and I4 read have to be done due to the restrictions
      ! of the READ function. I3 means to take the first 3 characters, I4 takes
      ! the first 4 characters of a string. If the 4th character is a blank
      ! space, then you get a "Bad value during integer read" error
      
      ! The issue with this method is if there are less than 1000 Scans for 
      ! a 1km pass, in which case this will fail. This is a possiblity (though
      ! has never been seen) for DB data if less than 5min of data is downloaded
      Num_Char = 0
      do i = 1, len(Temp_Char)
        if (Temp_Char(i:i) .ge. '0' .and. Temp_Char(i:i) .le. '9') then
          Num_Char = Num_Char + 1
        endif
      enddo
      
      if (Num_Char .eq. 1) THEN
        read(Temp_Char, fmt="(I1)") Num_Lines
      elseif (Num_Char .eq. 2) THEN
        read(Temp_Char, fmt="(I2)") Num_Lines   
      elseif (Num_Char .eq. 3) THEN
        read(Temp_Char, fmt="(I3)") Num_Lines
      elseif (Num_Char .eq. 4) THEN
        read(Temp_Char, fmt="(I4)") Num_Lines
      elseif (Num_Char .eq. 5) THEN
        read(Temp_Char, fmt="(I5)") Num_Lines
      endif
      
      
      !--- number of elements
      String_Index=index(Metadata_Temp, 'Max_EV_frames')
      Temp_Char = Metadata_Temp(String_Index+24:String_Index+24+5)

      Num_Char = 0
      do i = 1, len(Temp_Char) 
        if (Temp_Char(i:i) .ge. '0' .and. Temp_Char(i:i) .le. '9') then
          Num_Char = Num_Char + 1
        endif
      enddo

      if (Num_Char .eq. 1) THEN
        read(Temp_Char, fmt="(I1)") Num_Elements
      elseif (Num_Char .eq. 2) THEN
        read(Temp_Char, fmt="(I2)") Num_Elements
      elseif (Num_Char .eq. 3) THEN
        read(Temp_Char, fmt="(I3)") Num_Elements
      elseif (Num_Char .eq. 4) THEN
        read(Temp_Char, fmt="(I4)") Num_Elements
      elseif (Num_Char .eq. 5) THEN
        read(Temp_Char, fmt="(I5)") Num_Elements
      endif
		
		
		! exceptiopm atml
		if ( index(file_name_full, 'MYDATML') > 0 ) then
			print*,'WARNING: MYDATML file: number of lines and elements hardcoded in modis_module.f90 line 1330 '
			print*,'WARNING:              This should be read in from file  AW 04/24/2016'
			Num_elements = 271
			Num_lines=408
		
		end if

    end subroutine READ_MODIS_SIZE_ATTR
!==============================================================================
!
!==============================================================================
subroutine READ_MODIS_THERMAL_BAND(Chan_Idx, &
                                   Seg_Idx, &
                                   Eff_Wavenumber, &
                                   Number_Elements, &
                                   Number_Lines, &
                                   File_Path, &
                                   File_Name, &
                                   Radiance, &
                                   Brightness_Temp, &
                                   Uncertainty, &
                                   Error_Status)
   integer, intent(in):: Chan_Idx
   integer, intent(in):: Seg_Idx
   integer, intent(in):: Number_Elements
   integer, intent(in):: Number_Lines
   real, intent(in):: Eff_Wavenumber
   character(len=*), intent(in):: File_Path
   character(len=*), intent(in):: File_Name
   real, intent(out), dimension(:,:):: Radiance
   real, intent(out), dimension(:,:):: Brightness_Temp
   integer(kind=int1), intent(out), dimension(:,:):: Uncertainty
   integer(kind=int4), intent(out):: Error_Status

   call READ_MODIS_LEVEL1B(File_Path,File_Name, &
                           Chan_Idx, Radiance, Uncertainty, Number_Elements,Number_Lines, &
                           Seg_Idx,Image%Number_Of_Lines,Image%Number_Of_Lines_Read_This_Segment,Error_Status)
   call CONVERT_RADIANCE(Radiance,Eff_Wavenumber,Missing_Value_Real4)
   call COMPUTE_BT_ARRAY(Brightness_Temp,Radiance,Chan_Idx,Missing_Value_Real4)

end subroutine  READ_MODIS_THERMAL_BAND

end module MODIS_MODULE
