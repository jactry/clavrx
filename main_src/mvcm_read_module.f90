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

    use HDF, only: DFACC_read

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
    Search_String = 'IFFCMO_npp_'//Image%Level1b_Name(12:28)//'*.hdf'

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
    integer, intent(in):: Seg_Idx
    integer:: Sd_Id, Sds_Id, Status_Flag, Sds_Rank, Sds_Data_Type, Num_Attrs
    integer(kind=int4), dimension(3):: Sds_Dims, Sds_Stride, Sds_Start, Sds_Edges
    character(len=120):: Sds_Name, Sds_Name_Temp
    integer:: sfstart, sfginfo, sfrdata, sfend, sfendacc, sfn2index, sfselect
    integer(kind=int1), dimension(:,:,:), allocatable:: I1_Buffer
    integer:: Nx_Slab_Data, Nx_Slab_Read_Start, Nx_Slab_Read_End, Nx_Seg, Nx_Slab_Count
    integer:: Ny_Slab_Data, Ny_Slab_Read_Start, Ny_Slab_Read_End, Ny_Seg, Ny_Seg_Max, Ny_Slab_Count

    Status_Flag = 0
    Sds_Name = "Cloud_Mask"
    Cloud_Mask_Aux_Read_Flag = sym%NO   !will be set to yes if successful

    if (trim(Image%Auxiliary_Cloud_Mask_File_Name) == "no_file") then
         print *, EXE_PROMPT, MVCM_PROMPT, "MVCM File Could Not Be Found: "
         return
    endif

    !--- open file for read
    Sd_Id = sfstart(trim(Image%Level1b_Path)//trim(Image%Auxiliary_Cloud_Mask_File_Name), DFACC_read)

    !--- if file is unreadable, exit
    if (Sd_Id <= 0) then
         Status_Flag = 1
         print *, EXE_PROMPT, MVCM_PROMPT, "MVCM File Could Not Be Opened: "
         return
     endif

    !determine expected size of data slab
    Nx_Seg = Image%Number_of_Elements
    Ny_Seg_Max = Image%Number_of_Lines_Per_Segment
    Ny_Seg = Image%Number_of_Lines_Read_This_Segment

    !--- open sds
    Sds_Id = sfselect(Sd_Id, sfn2index(Sd_Id,trim(Sds_Name)))
    if (Sds_Id <= 0) then
         Status_Flag = 1
         Status_Flag = sfend(Sd_Id) + Status_Flag
         print *, EXE_PROMPT, MVCM_PROMPT, "MVCM Sds Could Not Be Opened: "
         return
     endif

    !--- get information
    Status_Flag = sfginfo(Sds_Id, Sds_Name_Temp, Sds_Rank, Sds_Dims,  &
                          Sds_Data_Type, Num_Attrs) + Status_Flag
    if (Status_Flag /= 0) then
         print *, EXE_PROMPT, MVCM_PROMPT, "MVCM Sds information could not be read: "
         Status_Flag = sfend(Sd_Id) + Status_Flag
         return
    endif

    !--- define size of data
    Nx_Slab_Data = Sds_Dims(1)
    Ny_Slab_Data = Sds_Dims(2)

    Nx_Slab_Read_Start = 0
    Nx_Slab_Read_End = Nx_Slab_Data-1
    Ny_Slab_Read_Start = (Seg_Idx-1)*Ny_Seg_Max
    Ny_Slab_Read_End = Ny_Slab_Read_Start + Ny_Seg -1
    !--- constrain to size of data
    Ny_Slab_Read_End = min(Ny_Slab_Read_End,Ny_Slab_Data)

    Nx_Slab_Count = Nx_Slab_Read_End - Nx_Slab_Read_Start + 1 
    Ny_Slab_Count = Ny_Slab_Read_End - Ny_Slab_Read_Start + 1 

    Sds_Stride = (/1, 1, 1/)
    Sds_Start = (/Nx_Slab_Read_Start, Ny_Slab_Read_Start, 0/)
    Sds_Edges = (/Nx_Slab_Count, Ny_Slab_Count,1/) 

    !--- compute number of lines to read for this segment
    if (.not. allocated(I1_Buffer)) allocate(I1_Buffer(Nx_Slab_Count, Ny_Slab_Count, 1))

    Status_Flag = sfrdata(Sds_Id, Sds_Start, Sds_Stride, Sds_Edges, I1_Buffer) + Status_Flag
    if (Status_Flag /= 0) then
         print *, EXE_PROMPT, MVCM_PROMPT, "MVCM Sds could not be read: "
         Status_Flag = sfend(Sd_Id) + Status_Flag
         return
    endif

    Status_Flag = sfendacc(Sds_Id) + Status_Flag
    if (Status_Flag /= 0) then
         print *, EXE_PROMPT, MVCM_PROMPT, "MVCM Sds could not be closed: "
         Status_Flag = sfend(Sd_Id) + Status_Flag
         return
    endif

    !--- close file
    Status_Flag = sfend(Sd_Id) + Status_Flag
    if (Status_Flag /= 0) then
         print *, EXE_PROMPT, MVCM_PROMPT, "MVCM file could not be closed: "
         return
    endif

    Cldmask%Cld_Mask_Aux(:,1:Ny_Slab_Count) = ishft(ishft(I1_Buffer(:,:,1),5),-6) 

    !--- switch CLAVR-x convection for mask
    Cldmask%Cld_Mask_Aux(:,1:Ny_Slab_Count) = 3-Cldmask%Cld_Mask_Aux(:,1:Ny_Slab_Count)

    if (allocated(I1_Buffer)) deallocate(I1_Buffer)

    Cloud_Mask_Aux_Read_Flag = sym%YES

 end subroutine READ_MVCM_DATA

!-------------------------------------------
!
!-------------------------------------------
end module MVCM_READ_MODULE
