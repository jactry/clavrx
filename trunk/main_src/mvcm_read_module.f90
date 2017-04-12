!$Id:$
!--------------------------------------------------------------------------------------
! This module reads the MVCM cloud mask to allow CLAVR-x to simulate the MODAWG chain
!
! Previously, these files were read in via the IFF module but we need to be able these
! files from the original NASA Level1b data.
!
!
!VGEOM_snpp_d20130426_t101800_c20161107162802.nc
!IFFCMO_npp_d20130426_t101800_c20161029052326_ssec_dev.hdf
!
!
!IFFCMO_aqua_d20121229_t051500_c20170114104543_ssec_dev.hdf 
! 
!-------------------------------------------------------------------------------------
module MVCM_READ_MODULE

    use PIXEL_COMMON, only: Image, Cldmask, Cloud_Mask_Aux_Flag, Cloud_Mask_Aux_Read_Flag
    use NUMERICAL_ROUTINES, only: COMPUTE_MONTH, COMPUTE_DAY, LEAP_YEAR_FCT
    use FILE_TOOLS, only: FILE_SEARCH
    use CONSTANTS, only: &
             Real4 &
           , Int4 &
           , Int2 &
           , Int1 &
           , Exe_Prompt &
           , Missing_Value_Real4 &
           , Sym &
           , Missing_Value_Int1

    implicit none

    public:: DETERMINE_MVCM_NAME
    public:: READ_MVCM_DATA

    character(len=13), parameter:: MVCM_PROMPT="MVCM_MODULE:"

    contains

!------------------------------------------------------------------------------------
! Determine the name of the MVCM file for the Level1b file
!------------------------------------------------------------------------------------
 subroutine DETERMINE_MVCM_NAME()

  character(len=100):: Search_String
  character(len=1020), dimension(:), pointer:: Files
  integer:: Num_Files
  character(len=4):: Year_String, Time_String
  character(len=2):: Month_String, Dom_String
  character(len=3):: Doy_String
  integer:: Year, Month, Dom, Doy,Ileap

  Image%Auxiliary_Cloud_Mask_File_Name = 'no_file'

  !--- NASA VIIRS Level1b
  if (index(Image%Level1b_Name,'VGEOM') == 1) then 

    !--- search should be the date and start time (ie. d20130426_t083000
    !Search_String = 'IFFCMO_npp_'//Image%Level1b_Name(12:28)//'*.hdf'
    Search_String = 'VCLDMK_snpp_'//Image%Level1b_Name(12:28)//'*.nc'

    Files => FILE_SEARCH(trim(Image%Level1b_Path),trim(Search_String),count=Num_Files)

    if (Num_Files == 0 .or. Num_Files > 1) then
        print *, EXE_PROMPT, MVCM_PROMPT, "MVCM File Not Found, "
        return
    endif

    Image%Auxiliary_Cloud_Mask_File_Name = Files(1)

    print *, EXE_PROMPT, MVCM_PROMPT, "NASA VIIRS Level1b MVCM File Found, ",trim(Image%Auxiliary_Cloud_Mask_File_Name)

  endif

  !--- NASA MODIS Level1b
  if (index(Image%Level1b_Name,'MYD021KM') == 1) then 

    !--- search should be the date and start time (ie. d20130426_t083000
    Search_String = 'IFFCMO_aqua_'//Image%Level1b_Name(11:27)//'*.hdf'
    Year_String = Image%Level1b_Name(11:14)
    Time_String = Image%Level1b_Name(19:22)
    Doy_String = Image%Level1b_Name(15:17)

    read(Doy_String,*) Doy
    read(Year_String,*) Year
    Ileap = LEAP_YEAR_FCT(Year)
    Month = COMPUTE_MONTH(Doy, ileap)
    Dom = COMPUTE_Day(Doy, ileap)
    write(Dom_String,fmt="(I2.2)") Dom
    write(Month_String,fmt="(I2.2)") Month

    Search_String = 'IFFCMO_aqua_d'//Year_String//Month_String//Dom_String//"_t"//Time_String//'*.hdf'

    Files => FILE_SEARCH(trim(Image%Level1b_Path),trim(Search_String),count=Num_Files)

    if (Num_Files == 0 .or. Num_Files > 1) then
        print *, EXE_PROMPT, MVCM_PROMPT, "Multiple NASA MODIS Level1b MVCM File Found, ",trim(Image%Auxiliary_Cloud_Mask_File_Name)
        return
    endif

    Image%Auxiliary_Cloud_Mask_File_Name = Files(1)

    print *, EXE_PROMPT, MVCM_PROMPT, "NASA MODIS Level1b MVCM File Found, ",trim(Image%Auxiliary_Cloud_Mask_File_Name)

  endif

  !--- SIPS IFF VIIRS Level1b
  !--- SIPS IFF MODIS Level1b

  
  Files => null() 

 end subroutine DETERMINE_MVCM_NAME
    
!------------------------------------------------------------------------------------------
! open, read a slab from the MVCM file and close it
!
! slab refers to the data in the file
! seg refers to the data segment needed for clavr-x
!------------------------------------------------------------------------------------------
 subroutine READ_MVCM_DATA(Seg_Idx)

    use NETCDF_READ_MODULE

    integer, intent(in):: Seg_Idx
    integer:: Sd_Id, Group_Id, Sds_Id, Status_Flag
    integer, dimension(2):: DimID
    integer(kind=int4), dimension(2):: Sds_Dims, Sds_Stride, Sds_Start, Sds_Edges
    character(len=120):: Sds_Name, Sds_Name_Temp
    integer(kind=int4), dimension(:,:), allocatable:: I4_Buffer
    integer:: Nx_Slab_Read_Start, Nx_Seg, Nx_Slab_Count
    integer:: Ny_Slab_Read_Start, Ny_Seg, Ny_Seg_Max, Ny_Slab_Count

    Status_Flag = 0
    Sds_Name = "Integer_Cloud_Mask"
    Cloud_Mask_Aux_Read_Flag = sym%NO   !will be set to yes if successful

    if (trim(Image%Auxiliary_Cloud_Mask_File_Name) == "no_file") then
         print *, EXE_PROMPT, MVCM_PROMPT, "MVCM File Could Not Be Found: "
         return
    endif

    !--- open file for read
    Status_Flag = nf90_open(trim(Image%Level1b_Path)//trim(Image%Auxiliary_Cloud_Mask_File_Name), &
                   nf90_nowrite, Sd_Id) + Status_Flag

    !--- if file is unreadable, exit
    if (Sd_Id <= 0) then
         Status_Flag = 1
         print *, EXE_PROMPT, MVCM_PROMPT, "MVCM File Could Not Be Opened: "
         return
     endif

    !--- get information
    Status_Flag = nf90_inq_ncid(Sd_Id, "geophysical_data", Group_Id) + Status_Flag
    Status_Flag = nf90_inq_varid(Group_Id, trim(Sds_Name), Sds_Id) + Status_Flag
    ! --- get dimention of file
    Status_Flag = nf90_inquire_variable(Group_Id, Sds_Id, dimids = DimID) + Status_Flag
    Status_Flag = nf90_inquire_dimension(Group_Id, DimID(1), len = Sds_Dims(1)) + Status_Flag
    Status_Flag = nf90_inquire_dimension(Group_Id, DimID(2), len = Sds_Dims(2)) + Status_Flag

    if (Status_Flag /= 0) then
         print *, EXE_PROMPT, MVCM_PROMPT, "MVCM Sds information could not be read: "
         return
    endif

    !determine expected size of data slab
    Nx_Seg = Image%Number_of_Elements
    Ny_Seg_Max = Image%Number_of_Lines_Per_Segment
    Ny_Seg = Image%Number_of_Lines_Read_This_Segment

    !--- define size of data
    Nx_Slab_Read_Start = 1
    Ny_Slab_Read_Start = (Seg_Idx-1) * Ny_Seg_Max + 1
    Nx_Slab_Count = Sds_Dims(1)
    Ny_Slab_Count = Ny_Seg

    !--- constrain to size of data
    if ((Seg_Idx * Ny_Seg) .gt. Sds_Dims(2)) Ny_Slab_Count = Sds_Dims(2) - Ny_Slab_Read_Start

    Sds_Stride = (/1, 1/)
    Sds_Start = (/Nx_Slab_Read_Start, Ny_Slab_Read_Start/)
    Sds_Edges = (/Nx_Slab_Count, Ny_Slab_Count/)

    if (.not. allocated(I4_Buffer)) allocate(I4_Buffer(Nx_Slab_Count, Ny_Slab_Count),stat=Status_Flag)
    call read_netcdf_2d_int(Group_Id, Sds_Start, Sds_Edges, trim(Sds_Name), I4_Buffer)

    if (Status_Flag /= 0) then
         print *, EXE_PROMPT, MVCM_PROMPT, "MVCM Sds could not be read: "
         return
    endif

    !--- close file
    Status_Flag = nf90_close(Sd_Id) + Status_Flag
    if (Status_Flag /= 0) then
         print *, EXE_PROMPT, MVCM_PROMPT, "MVCM file could not be closed: "
         return
    endif

    !--- switch CLAVR-x convection for mask
    Cldmask%Cld_Mask_Aux(:,1:Ny_Slab_Count) = 3-I4_Buffer(:,1:Ny_Slab_Count)

    ! --- set missing to CLAVR-x missing
    where (Cldmask%Cld_Mask_Aux .lt. 0 .or. Cldmask%Cld_Mask_Aux .gt. 3)
       Cldmask%Cld_Mask_Aux = Missing_Value_Int1
    endwhere

    if (allocated(I4_Buffer)) deallocate(I4_Buffer,stat=Status_Flag)

    Cloud_Mask_Aux_Read_Flag = sym%YES

 end subroutine READ_MVCM_DATA

!-------------------------------------------
!
!-------------------------------------------
end module MVCM_READ_MODULE
