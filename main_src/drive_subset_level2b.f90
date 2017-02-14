!$id:$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: drive_subset_level2b.f90 (src)
!       DRIVE_SUBSET_LEVEL2B (program)
!
! PURPOSE: Fortran code to subset level2b over the spatial dimensions and 
!          parameter contents
!
! DESCRIPTION: This code reads the spatial and parameter requirements from the 
!              drive_subset_level2b_input file
!
!              It will subset all level2b files in the path
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
program DRIVE_SUBSET_LEVEL2B
 use CX_CONSTANTS_MOD
 use HDF_PARAMS
 use NUMERICAL_TOOLS_MOD
 use SCALING_PARAMETERS,only:

 use FILE_TOOLS,only: &
  get_lun
 

 implicit none

 real(kind=real4):: Lon_West_Out
 real(kind=real4):: Lon_East_Out
 real(kind=real4):: Dlon_Out
 real(kind=real4):: Lat_North_Out
 real(kind=real4):: Lat_South_Out
 real(kind=real4):: Dlat_Out
 character(len=1020):: Path_Input
 character(len=1020):: File_Input
 character(len=1020):: Path_Output
 character(len=1020):: File_Output
 character(len=1020):: Data_Description_String
 character(len=128),dimension(:), allocatable:: Sds_Output_Names
 integer(kind=int4):: Nsds_Out
 integer(kind=int4):: Num_Files
 character(len=1020), dimension(:), allocatable:: Files
 integer(kind=int4):: Idx
 integer(kind=int4):: Num_Command_Line_Arguments
 integer(kind=int4):: Lun
 character(len=22), parameter:: LOCAL_EXE_PROMPT = "DRIVE_SUBSET_LEVEL2B: "


 !----------------------------------------------------------------------
 ! Begin Executable Code
 !----------------------------------------------------------------------
 Num_Command_Line_Arguments = iargc()
 if (Num_Command_Line_Arguments == 2) then
     call getarg(1, Path_Input)
     call getarg(2, Path_Output)
 else
     print *, LOCAL_EXE_PROMPT, "ERROR:: Unexpected Number of Command Line Arguments, stopping"
     stop
 endif

 if (Path_Input == Path_Output) then
     print *, LOCAL_EXE_PROMPT, "ERROR:: Input and Output paths are the same, stopping"
     stop
 endif

 Lun = GET_LUN()
 open(unit=Lun,file='drive_subset_level2b_input',status='old',action='read')
 read(unit=Lun,fmt=*) Lon_West_Out, Lon_East_Out, Dlon_Out
 read(unit=Lun,fmt=*) Lat_South_Out, Lat_North_Out, Dlat_Out
 read(unit=Lun,fmt=*) Data_Description_String
 read(unit=Lun,fmt=*) Nsds_Out
 allocate(Sds_Output_Names(Nsds_Out))
 do Idx = 1, Nsds_Out
   read(unit=Lun,fmt=*) Sds_Output_Names(Idx)
 enddo

!-----------------------------------------------------------------------------
! Search for level2b in the path
!-----------------------------------------------------------------------------
!call FILE_SEARCH(Path_Input,'level2b.hdf',Num_Files,Files)

File_Loop: do Idx = 1, Num_Files

        File_Input = Files(Idx)
        File_Output = File_Input
        print *, "processing ", trim(File_Input)
        call SUBSET_LEVEL2b(File_Input, &
                            Path_Input, &
                            File_Output, &
                            Path_Output, &
                            Lon_West_Out, &
                            Lon_East_Out, &
                            Dlon_Out, &
                            Lat_South_Out, &
                            Lat_North_Out, &
                            Dlat_Out, &
                            Data_Description_String, &
                            Sds_Output_Names)

enddo File_Loop
   
if (allocated(Files)) deallocate(Files)

deallocate(Sds_Output_Names)

!----------------------------------------------------------------------
! End Executable Code
!----------------------------------------------------------------------
end program DRIVE_SUBSET_LEVEL2B
