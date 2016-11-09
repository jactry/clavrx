!$Id:
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE
!
! NAME: oca_module.f90 (src)
!       OCA_MODULE (program)
!
! PURPOSE: This code finds and reads level2 cloud EUMETSAT product (OCA)
!          and saves them to CLAVR-x level2 output files.
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
! HISTORY:
!   Developed - Denis B. (03/18/2016)
!
!--------------------------------------------------------------------------------------
module OCA_MODULE

   use FILE_TOOLS, only: &
             FILE_SEARCH &
              , get_lun 

   use PIXEL_COMMON, only: &
             Sensor &
           , Cloud_Mask_Aux_Read_Flag

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
   private
   include 'hdf.f90'
   public :: READ_OCA
   

   integer(kind=int2), parameter, private:: fill_value = 32767
   real(kind=real4), parameter, private:: missing_value = -999.0
   character(len=13), parameter:: OCA_PROMPT="OCA_MODULE:"

contains

!--------------------------------------------------------------------
! read oca into memory
!--------------------------------------------------------------------
subroutine READ_OCA (L1b_Path, L1b_File, Nx, Ny, Seg_Idx, &
                       Ny_Total, Ny_Local_Temp, Cost_Out, &
                       Cloud_Press1_Out,Cloud_Press2_Out, &
                       Cloud_Press_Uncer1_Out,Cloud_Press_Uncer2_Out, &
                       Cod_Out,Phase_Out)

      character(len=*), intent(in):: L1b_Path
      character(len=*), intent(in):: L1b_File
      integer(kind=int4), intent(in):: Nx
      integer(kind=int4), intent(in):: Ny
      integer(kind=int4), intent(in):: Seg_Idx
      integer(kind=int4), intent(in):: Ny_Total
      integer(kind=int4), intent(out):: Ny_Local_Temp
      real(kind=real4), dimension(:,:), intent(out):: Cost_Out
      real(kind=real4), dimension(:,:), intent(out):: Cloud_Press1_Out
      real(kind=real4), dimension(:,:), intent(out):: Cloud_Press2_Out
      real(kind=real4), dimension(:,:), intent(out):: Cloud_Press_Uncer1_Out
      real(kind=real4), dimension(:,:), intent(out):: Cloud_Press_Uncer2_Out
      real(kind=real4), dimension(:,:), intent(out):: Cod_Out
      integer(kind=int1), dimension(:,:), intent(out):: Phase_Out

      character(len=1020):: Oca_File
      character(len=120):: Sds_Name
      character(len=120):: Sds_Name_Temp
      integer(kind=int4):: Iend
      integer(kind=int4):: Status_Flag
      integer(kind=int4):: Sd_Id
      integer(kind=int4):: Sds_Id
      integer(kind=int4):: sfend, sfstart, sfselect, &
                           sfginfo, sfn2index, sfrdata, sfendacc
      integer(kind=int4):: Num_Attrs, Sds_Data_Type, Sds_Rank
      integer(kind=int4), dimension(2):: Sds_Dims
      integer(kind=int4), dimension(2):: Sds_Start
      integer(kind=int4), dimension(2):: Sds_Stride
      integer(kind=int4), dimension(2):: Sds_Edges
      integer(kind=int4):: Nx_Local
      integer(kind=int4):: Ny_Local
      integer(kind=int4):: Ny_Start
      integer(kind=int4):: Ny_End
      integer:: Nx_Min
      integer:: Ny_Min
      real(kind=real4), allocatable, dimension(:,:):: R4_Buffer
      real(kind=real4), parameter:: Missing = -999.0
      integer(kind=int1), parameter:: Missing_I1 = -128

      ! --- Find OCA file name
      call DETERMINE_OCA_FILE(trim(L1b_Path),trim(L1b_File),Oca_File)

      !--- if file not found exit
      if (trim(Oca_File) == "no_file") RETURN

      Iend = 0
      Status_Flag = 0
      Sd_Id = 0

error_check: do while (Status_Flag == 0 .and. Iend == 0)

      !---- open file
      Sd_Id = sfstart(trim(L1b_Path)//trim(Oca_File), DFACC_read)

      !--- if file is unreadable, exit
      if (Sd_Id <= 0) then
         Status_Flag = 1
         exit
      endif

      !--- open sds
      Sds_Name = "Cost"
      Sds_Id = sfselect(Sd_Id, sfn2index(Sd_Id,trim(Sds_Name)))

      !--- get information
      Status_Flag = sfginfo(Sds_Id, Sds_Name_Temp, Sds_Rank, Sds_Dims,  &
                        Sds_Data_Type, Num_Attrs) + Status_Flag

!     print *, '================================'
!     print *, 'Sds_Id = ', Sds_Id
!     print *, 'sds_name_temp = ', trim(Sds_Name_Temp)
!     print *, 'sds_rank = ', Sds_Rank
!     print *, 'sds_dims = ', Sds_Dims
!     print *, 'num_attrs = ', Num_Attrs

      !--- compute number of lines to read for this segment
      Nx_Local = Sds_Dims(1)
      Ny_Local = Sds_Dims(2)
      Ny_Start = (Seg_Idx-1)*Ny + 1
      Ny_End = min(Ny_Start+Ny-1,Ny_Total)
      Ny_Local_Temp = Ny_End - Ny_Start + 1

      Nx_Min = min(Nx,Nx_Local)
      Ny_Min = min(Ny,Ny_Local_Temp)

      allocate(R4_Buffer(Nx_Local,Ny_Local_Temp))

      !-- read channel data for specified channel
      Sds_Stride = (/1, 1/)
      Sds_Start = (/0, Ny_Start-1/)
      Sds_Edges = (/Nx_Local, Ny_Local_Temp/)

      !--- read data
      Status_Flag = sfrdata(Sds_Id, Sds_Start, Sds_Stride, Sds_Edges, R4_Buffer) + Status_Flag

      !--- reset missing values to clavr-x standard
      where (R4_Buffer .ge. 2000)
         R4_Buffer = Missing
      endwhere

      !--- save data to output
      Cost_Out(1:Nx_Min,1:ny_Min) = R4_Buffer(1:Nx_Min,1:Ny_Min)
      R4_Buffer = Missing
      Status_Flag = sfendacc(Sds_Id) + Status_Flag

      !--- open sds Layer1_CTP
      Sds_Name = "Layer1_CTP"
      Sds_Id = sfselect(Sd_Id, sfn2index(Sd_Id,trim(Sds_Name)))

      !--- read data
      Status_Flag = sfrdata(Sds_Id, Sds_Start, Sds_Stride, Sds_Edges, R4_Buffer) + Status_Flag

      !--- reset missing values to clavr-x standard
      where (R4_Buffer .gt. 2000)
         R4_Buffer = Missing
      endwhere

      !--- save data to output
      Cloud_Press1_Out(1:Nx_Min,1:ny_Min) = R4_Buffer(1:Nx_Min,1:Ny_Min)
      R4_Buffer = Missing
      Status_Flag = sfendacc(Sds_Id) + Status_Flag

      !--- open sds Layer2_CTP
      Sds_Name = "Layer2_CTP"
      Sds_Id = sfselect(Sd_Id, sfn2index(Sd_Id,trim(Sds_Name)))

      !--- read data
      Status_Flag = sfrdata(Sds_Id, Sds_Start, Sds_Stride, Sds_Edges, R4_Buffer) + Status_Flag

      !--- reset missing values to clavr-x standard
      where (R4_Buffer .gt. 2000)
         R4_Buffer = Missing
      endwhere

      !--- save data to output
      Cloud_Press2_Out(1:Nx_Min,1:ny_Min) = R4_Buffer(1:Nx_Min,1:Ny_Min)
      R4_Buffer = Missing
      Status_Flag = sfendacc(Sds_Id) + Status_Flag

      !--- open sds Layer1_CTP_Err
      Sds_Name = "Layer1_CTP1_Err"
      Sds_Id = sfselect(Sd_Id, sfn2index(Sd_Id,trim(Sds_Name)))

      !--- read data
      Status_Flag = sfrdata(Sds_Id, Sds_Start, Sds_Stride, Sds_Edges, R4_Buffer) + Status_Flag

      !--- reset missing values to clavr-x standard
      where (R4_Buffer .ge. 2000)
         R4_Buffer = Missing
      endwhere

      !--- save data to output
      Cloud_Press_Uncer1_Out(1:Nx_Min,1:ny_Min) = R4_Buffer(1:Nx_Min,1:Ny_Min)
      R4_Buffer = Missing
      Status_Flag = sfendacc(Sds_Id) + Status_Flag

      !--- open sds Layer2_CTP_Err
      Sds_Name = "Layer2_CTP_Err"
      Sds_Id = sfselect(Sd_Id, sfn2index(Sd_Id,trim(Sds_Name)))

      !--- read data
      Status_Flag = sfrdata(Sds_Id, Sds_Start, Sds_Stride, Sds_Edges, R4_Buffer) + Status_Flag

      !--- reset missing values to clavr-x standard
      where (R4_Buffer .ge. 2000)
         R4_Buffer = Missing
      endwhere

      !--- save data to output
      Cloud_Press_Uncer2_Out(1:Nx_Min,1:ny_Min) = R4_Buffer(1:Nx_Min,1:Ny_Min)
      R4_Buffer = Missing
      Status_Flag = sfendacc(Sds_Id) + Status_Flag

      !--- open sds Layer1_COT
      Sds_Name = "Layer1_COT"
      Sds_Id = sfselect(Sd_Id, sfn2index(Sd_Id,trim(Sds_Name)))

      !--- read data
      Status_Flag = sfrdata(Sds_Id, Sds_Start, Sds_Stride, Sds_Edges, R4_Buffer) + Status_Flag

      !--- reset missing values to clavr-x standard
      where (R4_Buffer .ge. 2000)
         R4_Buffer = Missing
      endwhere

      !--- save data to output
      Cod_Out(1:Nx_Min,1:ny_Min) = R4_Buffer(1:Nx_Min,1:Ny_Min)
      R4_Buffer = Missing
      Status_Flag = sfendacc(Sds_Id) + Status_Flag

      !--- open sds Phase
      Sds_Name = "Phase"
      Sds_Id = sfselect(Sd_Id, sfn2index(Sd_Id,trim(Sds_Name)))

      !--- read data
      Status_Flag = sfrdata(Sds_Id, Sds_Start, Sds_Stride, Sds_Edges, R4_Buffer) + Status_Flag
!print *,R4_Buffer(2120:2130,420:440)
!print *,'---------------'

      !--- reset missing values to clavr-x standard
      where (R4_Buffer .ge. 2000)
         R4_Buffer = float(Missing_I1)
      endwhere

      !--- re-map OCA phase to CLAVR-x
      !Phase" Range 110 to 113
      !110 - no cloud -> 0 - CLEAR_PHASE
      !111 single layer water -> 1 - WATER_PHASE
      !112 single layer ice -> 4 - ICE_PHASE
      !113 two-layer (ice over ?) -> 4 - ICE_PHASE
      where (R4_Buffer .eq. 110.)
         R4_Buffer = 0
      endwhere
      where (R4_Buffer .eq. 111.)
         R4_Buffer = 1
      endwhere
      where (R4_Buffer .eq. 112. .or. R4_Buffer .eq. 113.)
         R4_Buffer = 4
      endwhere

      !--- save data to output
      Phase_Out(1:Nx_Min,1:ny_Min) = int(R4_Buffer(1:Nx_Min,1:Ny_Min))
      R4_Buffer = Missing
!print *,Phase_Out(2120:2130,420:440)
!print *,'---------------'
!stop
      Status_Flag = sfendacc(Sds_Id) + Status_Flag


      !--- close file
      Status_Flag = sfend(Sd_Id) + Status_Flag

      iend = 1

      enddo error_check  !end of while loop

      !--- deallocate memory
      !--- clean up memory
      if (allocated(R4_Buffer)) deallocate(R4_Buffer)

      if (Status_Flag /= 0) then
         print *, EXE_PROMPT, "Reading of OCA Level2 failed, skipping this file"

         !--- close hdf file, if exiting after opening it
         if (Sd_Id > 0)  then
            Status_Flag = sfend(Sd_Id)
         endif
      endif

end subroutine READ_OCA

!--------------------------------------------------------------------
! find oca filename and check if it exists
!--------------------------------------------------------------------
subroutine DETERMINE_OCA_FILE(Path_In,File_In,File_Out)

      character(len=*), intent(in):: Path_In
      character(len=*), intent(in):: File_In
      character(len=*), intent(out):: File_Out

      character(len=500):: Search_String
      character(len=1020), dimension(:), pointer:: Files
      integer(kind=int4):: Num_Files


      !0        1         2         3         4         5
      !12345678901234567890123456789012345678901234567890
      !met10_1_2016_053_1300.area
      !met10_2016_053_1300.oca.hdf

      Search_String = trim(File_In(1:6))//trim(File_In(9:22))//'oca.hdf'

      Files => FILE_SEARCH(trim(Path_In),trim(Search_String),count=Num_Files)

      if (Num_Files == 0) then
         print *, EXE_PROMPT, OCA_PROMPT, "No OCA File Found"
         File_Out = "no_file"
         return
      endif

      if (Num_Files > 1) then
         print *, EXE_PROMPT, OCA_PROMPT, "Multiple OCA Files Found"
!         File_Out = "no_file"
      endif

      File_Out = Files(1)
      print *, EXE_PROMPT, OCA_PROMPT, "OCA File Found, ",trim(File_Out)

      Files => null()

end subroutine DETERMINE_OCA_FILE

!=====================================================================

end module OCA_MODULE

