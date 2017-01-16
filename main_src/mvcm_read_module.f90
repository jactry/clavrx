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
!-------------------------------------------------------------------------------------
module MVCM_READ_MODULE

    use PIXEL_COMMON, only: Image, Cldmask
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

  Image%Auxiliary_Cloud_Mask_File_Name = 'no_file'

  !--- NASA VIIRS Level1b
  if (index(Image%Level1b_Name,'VGEOM') == 1) then 

    Search_String = 'IFFCMO_npp_'//Image%Level1b_Name(8:28)//'*.hdf'

    Files => FILE_SEARCH(trim(Image%Level1b_Path),trim(Search_String),count=Num_Files)

    if (Num_Files == 0 .or. Num_Files > 1) then
        print *, EXE_PROMPT, MVCM_PROMPT, "MVCM File Not Found, "
        return
    endif

    Image%Auxiliary_Cloud_Mask_File_Name = Files(1)

    print *, EXE_PROMPT, MVCM_PROMPT, "MVCM File Found, ",trim(Image%Auxiliary_Cloud_Mask_File_Name)

  endif
  !--- NASA MODIS Level1b
  !--- SIPS IFF VIIRS Level1b
  !--- SIPS IFF MODIS Level1b

  
  Files => null() 

 end subroutine DETERMINE_MVCM_NAME
    
!------------------------------------------------------------------------------------------
! open, read a slab from the MVCM file and close it
!------------------------------------------------------------------------------------------
 subroutine READ_MVCM_DATA(Seg_Idx)
    integer, intent(in):: Seg_Idx
    integer:: Sd_Id, Sds_Id, Status_Flag, Sds_Rank, Sds_Data_Type, Num_Attrs
    integer(kind=int4), dimension(2):: Sds_Dims, Sds_Stride, Sds_Start, Sds_Edges
    character(len=120):: Sds_Name, Sds_Name_Temp
    integer:: sfstart, sfginfo, sfrdata, sfend, sfendacc
    integer(kind=int2), dimension(:,:), allocatable:: I2_Buffer

    Status_Flag = 0
    Sds_Name = "Cloud_Mask"

    if (trim(Image%Auxiliary_Cloud_Mask_File_Name) == "no_file") then
         print *, EXE_PROMPT, MVCM_PROMPT, "MVCM File Could Not Be Found: "
         return
    endif

    Sd_Id = sfstart(trim(Image%Level1b_Path)//trim(Image%Auxiliary_Cloud_Mask_File_Name), DFACC_read)

    !--- if file is unreadable, exit
    if (Sd_Id <= 0) then
         Status_Flag = 1
         print *, EXE_PROMPT, MVCM_PROMPT, "MVCM File Could Not Be Opened: "
         return
     endif

    !determine expected size of data slab
    Nx = Image%Number_of_Elements

    !--- open sds
    Sds_Id = sfselect(Sd_Id, sfn2index(Sd_Id,trim(Sds_Name)))
    if (Sds_Id <= 0) then
         Status_Flag = 1
         print *, EXE_PROMPT, MVCM_PROMPT, "MVCM Sds Could Not Be Opened: "
         return
     endif

    !--- get information
    Status_Flag = sfginfo(Sds_Id, Sds_Name_Temp, Sds_Rank, Sds_Dims,  &
                          Sds_Data_Type, Num_Attrs) + Status_Flag
    if (Status_Flag /= 0) then
         print *, EXE_PROMPT, MVCM_PROMPT, "MVCM Sds information could not be read: "
         return
    endif

    !--- define size of data
    Nx_Local = Sds_Dims(1)
    Ny_Local = Sds_Dims(2)

    Sds_Stride = (/1, 1/)
    Sds_Start = (/0, Ny_Start-1/)
    Sds_Edges = (/Nx_Local, Ny_Local_Temp/) 

    !--- compute number of lines to read for this segment
    Ny_Start = (Seg_Idx-1)*Ny + 1
    Ny_End = min(Ny_Start+Ny-1,Ny_Total)
    Ny_Local_Temp = Ny_End - Ny_Start + 1

    Nx_Min = min(Nx,Nx_Local)
    Ny_Min = min(Ny,Ny_Local_Temp)

    if (.not. allocated(I2_Buffer)) allocate(I2_Buffer(Nx_Local, Ny_Local_Temp))

    Status_Flag = sfrdata(Sds_Id, Sds_Start, Sds_Stride, Sds_Edges, I2_Buffer) + Status_Flag
    if (Status_Flag /= 0) then
         print *, EXE_PROMPT, MVCM_PROMPT, "MVCM Sds could not be read: "
         return
    endif

    Cldmask%Cld_Mask_Aux(1:Nx_Min,1:Ny_Min) = I2_Buffer(1:Nx_Min,1:Ny_Min) 

    Status_Flag = sfendacc(Sds_Id) + Status_Flag
    if (Status_Flag /= 0) then
         print *, EXE_PROMPT, MVCM_PROMPT, "MVCM Sds could not be closed: "
         return
    endif

    !--- close file
    Status_Flag = sfend(Sd_Id) + Status_Flag
    if (Status_Flag /= 0) then
         print *, EXE_PROMPT, MVCM_PROMPT, "MVCM file could not be closed: "
         return
    endif

    if (allocated(I2_Buffer)) deallocate(I2_Buffer)

 end subroutine READ_MVCM_DATA

!-------------------------------------------
!
!-------------------------------------------
end module MVCM_READ_MODULE
