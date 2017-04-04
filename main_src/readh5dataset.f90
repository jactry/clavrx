!$Id$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: readh5dataset.f90 (src)
!       ReadH5Dataset (program)
!
! PURPOSE: interface to HDF-5 reading routines
!
! DESCRIPTION:  see below
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
! REVISION HISTORY:
!  A.Walther 2013: added all substring tools
!   (AW) 17 Oct 2013: fixed bug in H5ReadInteger2DSub ( use of counts instead of dims)
!
!
!  H5ReadDataSet
!  =============
! 
!  This module contains generic subroutines for reading the most common
!  data and attribute types in an hdf-5 type file:
!  a) 0-, 1-, 2- or 3-dimensional datasets or attributes, which can be
!     of integer, real or double precision data type (0-dim. means here
!     a scalar).
!  b) 0- and 1-dimensional character strings in attributes or datasets
!     (0-dim. means here a single string, 1-dim. a 1-D array of strings).
!  c) 0- and 1-dimensional integer, real or double parts of compound
!     attributes or datasets (0-dim. means here a scalar); note that
!     the different parts of a compound set have to be read separately.
! 
!  Note: Depending on your machine, specify the single and double
!        precision kind parameters (SP and DP) below.
! 
!  The routines handle datasets of groups, attributes of groups, as well
!  as attributes of datasets of groups:
!      /<N_groups>/dataset
!      /<N_groups>/attribute
!      /<N_groups>/dataset/attribute
!  where 'N' is 0, 1, 2, ...
! 
!  To use these routines, add a "USE ReadH5dataset" statement to the
!  calling routine. The calling routine must contain a variable of the
!  same type and dimension as the dataset or attribute. The exact sizes
!  of arrays are determined automatically.
! 
!  When reading a dataset or attribute from an HDF 5 file, call the
!  generic routine H5ReadDataset:
! 
!    CALL H5ReadDataset ( filename, dsetname, dataset )
!    CALL H5ReadDataset ( filename, dsetname, compoundsetname, dataset )
! 
!  respectively for normal datasets and attributes, and for parts of
!  compound datasets and attributes
! 
!  Alternatively, one can make the following generic calls (which sound
!  better when reading an attribute) with exactly the same functionality: 
! 
!    CALL H5ReadAttribute ( filename, attribname, attribute )
!    CALL H5ReadAttribute ( filename, attribname, compoundsetname, attribute )
! 
!  If an error occurs in either of the steps of the access to the
!  data or reading the data, the variable ErrorFlag is set to -1
!  and the ErrorMessage string is filled. These two can be used
!  for error tracing in the calling routine.
!  Note that some routines return a possitive ErrorFlag number,
!  which does _not_ indicate an error (very confusing).
! 
! 
!  In the above calls are:
!    'filename'   = a character type variable containing the name of the
!                   HDF-5 file
!    'dsetname'   = the character type name of the dataset or attribute,
!                   including the full HDF-5 group hierarchy, using slashes
!                   as the divider symbol (e.g. "/DATA/GEOMETRY/radius");
!                   be sure to use the right case on the hierarchy names
!    'attribname' = same as 'dsetname' in functionality
!    'compoundsetname' = name of the part of a compound set to read
!    'dataset'    = a pointer type variable of the right rank and type,
!                   that will contain the dataset or attribute
!    'attribute'  = same as 'dataset' in functionality
! 
!  Additionally available outside the module is the LengString function,
!  which determines the length of a string, that is: the last non-space
!  character in a string.
!
!   HISTORY: 
!      2014:Added sub array routines (2104) AW
!      04/08/2015: close dsapce and memspace in H5Read2Integer2DSub, which caused memory leak.
!                 may be also have to close in other routines (AW)
!--------------------------------------------------------------------------------------
MODULE ReadH5Dataset
! *******************************************************************
! *** Interface and Module Declaration Section                    ***
! *******************************************************************

  !<<<<<<<<<<<<<<<<<<<<<<<<< Used Modules >>>>>>>>>>>>>>>>>>>>>>>>>>>
  USE HDF5

  !<<<<<<<<<<<<<<<<<<<<<<<<< Interfaces >>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  ! Close all for the outside world
  PRIVATE

  ! Make only this available to the outside world
  PUBLIC :: H5ReadDataset, H5ReadAttribute
  PUBLIC :: ErrorFlag, ErrorMessage, LengString
  PUBLIC :: H5_DATASET_DIMENSIONS

  ! The interfaces themselves
  INTERFACE H5ReadDataset
     MODULE PROCEDURE  &
      & H5ReadInteger0D, H5ReadInteger1D, H5ReadInteger2D, H5ReadInteger3D, &
      & H5ReadReal0D,    H5ReadReal1D,    H5ReadReal2D,    H5ReadReal3D,    &
      & H5ReadDouble0D,  H5ReadDouble1D,  H5ReadDouble2D,  H5ReadDouble3D,  &
      & H5ReadString0D,  H5ReadString1D,  H5ReadI8nteger1D    ,                              &
      & H5ReadCompInteger0D, H5ReadCompReal0D, H5ReadCompDouble0D,          &
      & H5ReadCompInteger1D, H5ReadCompReal1D, H5ReadCompDouble1D,   &
      & H5ReadReal2DSub, H5ReadInteger2DSub, H5Read2Integer2D, H5Read2Integer2DSub
  END INTERFACE

  INTERFACE H5ReadAttribute
     MODULE PROCEDURE  &
      & H5ReadInteger0D, H5ReadInteger1D, H5ReadInteger2D, H5ReadInteger3D, &
      & H5ReadReal0D,    H5ReadReal1D,    H5ReadReal2D,    H5ReadReal3D,    &
      & H5ReadDouble0D,  H5ReadDouble1D,  H5ReadDouble2D,  H5ReadDouble3D,  &
      & H5ReadString0D,  H5ReadString1D,  H5ReadI2nteger0D , &                                 
      & H5ReadCompInteger0D, H5ReadCompReal0D, H5ReadCompDouble0D,          &
      & H5ReadCompInteger1D, H5ReadCompReal1D, H5ReadCompDouble1D
  END INTERFACE

  !<<<<<<<<<<<<<<<<<<<<<<<<< Module variables >>>>>>>>>>>>>>>>>>>>>>>

  INTEGER, PARAMETER  :: DP = 8  ! byte precision of double variables
  INTEGER, PARAMETER  :: SP = 4  ! byte precision of real variables

  INTEGER(HID_T)      :: file_id
  INTEGER(HID_T)      :: root_id
  CHARACTER(LEN = 1)  :: rootname
  INTEGER(HID_T)      :: id
  INTEGER(HID_T)      :: p_id
  INTEGER(HID_T)      :: d_id
  INTEGER(HID_T)      :: type_id
  INTEGER             :: dspace
  integer, parameter :: int8 = selected_int_kind(10)
  INTEGER(HID_T), DIMENSION(:), ALLOCATABLE :: OpenLevels_id
  INTEGER(HID_T), DIMENSION(:), ALLOCATABLE :: OpenLevels_type
  INTEGER                                   :: NLevels

  ! Error handling -- made available for use in the calling routines
  INTEGER             :: ErrorFlag
  CHARACTER(LEN=100)  :: ErrorMessage
  INTEGER             :: AllocStat

  ! For debugging the routines.
  ! If set true here all routines show debugging messages.
  ! If set false here only routines that have their own debug
  ! activated show debugging messages.
  !LOGICAL             :: DebugMsg = .true.
  LOGICAL             :: DebugMsg = .false.


! *******************************************************************


CONTAINS


   ! **************************** 
   !    Returns dataset dimensions
   ! *******************************

   subroutine H5_DATASET_DIMENSIONS  &
        ( &
        h5filename &
        , datasetname &
        , dims )
      implicit none
      
      character ( len = *) :: h5filename
      character ( len = *) :: datasetname
      integer, pointer :: dims(:)
      
      INTEGER   :: ltype
      
      
      INTEGER, PARAMETER                                   :: maxdims = 4
      INTEGER                                              :: ndims
      INTEGER(hsize_t), DIMENSION(maxdims) , target        :: datadims
      INTEGER(hsize_t), DIMENSION(maxdims)                 :: maxdatadims
      
     
      CALL H5Read_init ( H5filename, datasetname )
      IF (ErrorFlag.lt.0) return
      
      ltype=OpenLevels_type(NLevels)

      IF ( ltype == H5G_DATASET_F) THEN
         CALL h5dget_space_f(d_id,dspace,ErrorFlag)
      ELSE  !  ltype.eqv.3
         CALL h5aget_space_f(d_id,dspace,ErrorFlag)
      ENDIF
      
      
      CALL h5sget_simple_extent_ndims_f(dspace,ndims,ErrorFlag)
     
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace dimensionality"
       return
    ENDIF
    CALL DebugMessage("     > No. of dim.: ",ndims) 

    ALLOCATE(dims(ndims),stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error allocating dims"
       return
    ENDIF

    CALL h5sget_simple_extent_dims_f(dspace,datadims,maxdatadims,ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace size"
       return
    ENDIF
      
      
      dims = int(datadims(1:ndims))
      CALL H5Read_close
      
   
   end subroutine H5_DATASET_DIMENSIONS

! *******************************************************************
! *** General Routines                                            ***
! *******************************************************************


  ! Initialise the HDF file and locate the dataset or attribute
  ! to read, opening all levels along the way.

  SUBROUTINE H5Read_init          &
       (                          &
       H5filename,                & !Coming in
       setname                    & !Coming in
       )

    !<<<<<<<<<<<<<<<<<<<<<<< Implicit statement >>>>>>>>>>>>>>>>>>>>>
    IMPLICIT NONE

    !<<<<<<<<<<<<<<<<<<<<<<< Arguments >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !Coming in
    CHARACTER (len = *)              :: H5filename
    CHARACTER (len = *)              :: setname

    !Going out:
    !None

    !<<<<<<<<<<<<<<<<<<<<<<< Local variables >>>>>>>>>>>>>>>>>>>>>>>>
    CHARACTER(LEN = LEN(trim(setname))+1)  :: local_setname
    CHARACTER(LEN = LEN(trim(setname))+1)  :: subsetname
    INTEGER                          :: ilevel
    INTEGER                          :: ibeg,iend,iendnext
    INTEGER                          :: ltype,ltypenext
    LOGICAL                          :: LocalDebug


    ! For debugging, in addition to the debugging of the calling routine
    LocalDebug = .false.
   ! LocalDebug = .true.

    !<<<<<<<<<<<<<<<<<<<<<<< Start of routine code >>>>>>>>>>>>>>>>>>

    ! No errors yet

    ErrorFlag=0
    ErrorMessage=""
    
    

    ! Copy the setname to a local variable, making sure it start with a '/'

    local_setname = trim(setname)
    IF ( local_setname(1:1) .NE. "/") local_setname = "/"//setname

    ! Initialise the HDF routines and open the file

    CALL DebugMessage(" >>> Showing debug messages of HDF5 file reading <<<")
    CALL h5open_f(ErrorFlag)

    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error initialising HDF routines"
       return
    ENDIF

    CALL DebugMessage(" --- Opening HDF file")

    CALL h5fopen_f(H5filename, H5F_ACC_RDONLY_F, file_id, ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error opening HDF file"
       return
    ENDIF

    ! Open the root group.

    CALL DebugMessage(" --- Opening root group")
    rootname = "/"
    CALL h5gopen_f(file_id,rootname,root_id,ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error opening root group"
       return
    ENDIF

    p_id = root_id  ! the rootgroup was just opened.
    ibeg=2          ! the slash before each level name is to be ignored

    ! First determine how many levels there are to be opened.

    NLevels=0
    do ilevel=1,LengString(local_setname)
       if (local_setname(ilevel:ilevel) == "/") NLevels=NLevels+1
    enddo
    CALL DebugMessage(" --- No. of levels: ",NLevels)

    ! Then open all levels, storing what is opened so that we can
    ! close things at the end properly.
 
  
    ALLOCATE(OpenLevels_id(NLevels),stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error allocating OpenLevels_id"
      
       return
    ENDIF
 
    ALLOCATE(OpenLevels_type(NLevels),stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error allocating OpenLevels_type"
       return
    ENDIF

    ! The very last level is the one to be read: dataset or attribute
    ! The level before that is either a group or a dataset.

    ! The very first level (after the root group) could be a group,
    ! a dataset or an attribute.

    ! Locate the end of the first element name
    iend=INDEX(local_setname(ibeg:),"/")
   
    IF (iend == 0) iend=LengString(subsetname(ibeg:))+ ibeg - 1
 
    call CheckMembers("/", &
             local_setname(ibeg:iend),ltypenext)
    IF (ErrorFlag.lt.0) return
    ltype=ltypenext  !  type of the level to be opened

    ! If there is only one level, there is no group after the root
    ! group, in which case the following do-loop is skipped.
 
    DO ilevel=1,NLevels-1
     
       ! Locate the end of the current element name
       iend=INDEX(local_setname(ibeg:),"/")

       IF (LocalDebug) THEN
          print*,""
          print*,local_setname(1:LengString(local_setname))
          print*,'il=',ilevel,'->',ibeg,iend,local_setname(ibeg:iend)
       ENDIF

       ! Name of the element to be opened in the next round
       subsetname=local_setname(iend+1:)
       iendnext=INDEX(subsetname(ibeg:),"/")
       IF (iendnext == 0) iendnext=LengString(subsetname(ibeg:))+ibeg-1

       IF (LocalDebug) THEN
          print*,subsetname
          print*,'  =',ibeg,iendnext,subsetname(ibeg:iendnext)
       ENDIF

       IF (ltype == H5G_GROUP_F) THEN

          CALL CheckMembers(local_setname(ibeg:iend), &
                subsetname(ibeg:iendnext),ltypenext)
          IF (ErrorFlag.lt.0) return

       ELSE  !  ltype.eqv.H5G_DATASET_F

          ltypenext=3 ! attribute = only possibility here

       ENDIF

       ! Open the current level

       IF (ltype == H5G_GROUP_F) THEN

          CALL DebugMessage( " --- Opening group for passing through: " &
                          // local_setname(ibeg:iend) )
          CALL h5gopen_f(p_id,local_setname(ibeg:iend),id,ErrorFlag)
          IF (ErrorFlag.lt.0) THEN
             ErrorMessage=" *** Error opening group " &
                        // local_setname(ibeg:iend)
             return
          ENDIF

       ELSE

          CALL DebugMessage( " --- Opening dataset for passing through: " &
                          // local_setname(ibeg:iend) )
          CALL h5dopen_f(p_id,local_setname(ibeg:iend),id,ErrorFlag)
          IF (ErrorFlag.lt.0) THEN
             ErrorMessage=" *** Error opening dataset " &
                        // local_setname(ibeg:iend)
             return
          ENDIF

       ENDIF

       ! Store details, so that we can close things properly at the end.
       OpenLevels_id(ilevel)=id
       OpenLevels_type(ilevel)=ltype

       local_setname=local_setname(iend+1:)
       p_id=id
       ltype=ltypenext

    ENDDO

    ! We have arrived at the level (dataset or atrribute) to open
    ! for reading.

    subsetname=local_setname(ibeg:LengString(local_setname))
    ltype=ltypenext

    IF (LocalDebug) THEN
       print*,"last level reached: ",subsetname(1:LengString(subsetname))
       print*,"ltype = ",ltype
    ENDIF

    IF ( ltype == H5G_DATASET_F ) THEN

       CALL DebugMessage( " --- Opening dataset for reading: " &
                       // subsetname(1:LengString(subsetname)) )
       CALL h5dopen_f(p_id, subsetname, d_id,ErrorFlag)
       IF (ErrorFlag.lt.0) THEN
          ErrorMessage=" *** Error opening dataset " &
                     // local_setname(ibeg:iend)
          return
       ENDIF

    ELSE  !  ltype.eqv.3

       CALL DebugMessage( " --- Opening attribute for reading: " &
                       // subsetname(1:LengString(subsetname)) )
                       
       CALL h5aopen_name_f(p_id, subsetname, d_id,ErrorFlag)
       IF (ErrorFlag.lt.0) THEN
          ErrorMessage=" *** Error opening attribute " &
                     // local_setname(ibeg:iend)
          return
       ENDIF

    ENDIF

    OpenLevels_id(NLevels)=d_id
    OpenLevels_type(NLevels)=ltype

  END SUBROUTINE H5Read_init

! *******************************************************************

  ! Close all groups, opened when searching for the dataset or
  ! attribute to read

  SUBROUTINE H5Read_close

    !<<<<<<<<<<<<<<<<<<<<<<< Implicit statement >>>>>>>>>>>>>>>>>>>>>
    IMPLICIT NONE

    !<<<<<<<<<<<<<<<<<<<<<<< Arguments >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !Coming in
    !None

    !Going out:
    !None

    !<<<<<<<<<<<<<<<<<<<<<<< Local variables >>>>>>>>>>>>>>>>>>>>>>>>
    INTEGER            :: n, ltype

    !<<<<<<<<<<<<<<<<<<<<<<< Start of routine code >>>>>>>>>>>>>>>>>>

    CALL DebugMessage(" --- Closing open levels")

    IF ( ALLOCATED(OpenLevels_id) ) THEN

       ! Let us close things in reversed order: seems nicest.

       DO n = NLevels,1,-1

          id = OpenLevels_id(n)
          ltype=OpenLevels_type(n)

          if ( ltype == H5G_GROUP_F ) then
             CALL h5gclose_f(id,ErrorFlag)
          else if ( ltype == H5G_DATASET_F ) then
             CALL h5dclose_f(id,ErrorFlag)
          else  !  ( ltype.eqv.3 )
             CALL h5aclose_f(id,ErrorFlag)
          endif
          IF (ErrorFlag.lt.0) THEN
             ErrorMessage=" *** Error closing level "
             return
          ENDIF

       END DO
  
       DEALLOCATE(OpenLevels_type,stat=AllocStat)
       IF ( AllocStat.ne.0 ) THEN
          ErrorFlag=-1
          ErrorMessage=" *** Error deallocating OpenLevels_type"
          return
       ENDIF

       DEALLOCATE(OpenLevels_id,stat=AllocStat)
       IF ( AllocStat.ne.0 ) THEN
          ErrorFlag=-1
          ErrorMessage=" *** Error deallocating OpenLevels_id"
          return
       ENDIF

    END IF

    ! Close the root group

    CALL DebugMessage(" --- Closing root group")
    CALL h5gclose_f(root_id,ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error closing root group"
       return
    ENDIF

    ! Close the HDF file and HDF interface

    CALL DebugMessage(" --- Closing HDF file")
    CALL h5fclose_f(file_id,ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error closing HDF file"
       return
    ENDIF

    CALL h5close_f(ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error closing HDF routines"
       return
    ENDIF

  END SUBROUTINE H5Read_close

! *******************************************************************

  ! Determine a length of a text in a given string by going backwards.
  ! The 'len' command gives the defined length of the string.

  SUBROUTINE CheckMembers(levelstr,sublevelstr,ltype)

    !<<<<<<<<<<<<<<<<<<<<<<< Implicit statement >>>>>>>>>>>>>>>>>>>>>
    IMPLICIT NONE

    !<<<<<<<<<<<<<<<<<<<<<<< Arguments >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !Coming in
    CHARACTER (len = *), INTENT(in) :: levelstr,sublevelstr

    !Going out:
    INTEGER, INTENT(out)            :: ltype

    !<<<<<<<<<<<<<<<<<<<<<<< Local variables >>>>>>>>>>>>>>>>>>>>>>>>
    INTEGER            :: nmembers, n, obj_type
    LOGICAL            :: FOUND
    CHARACTER(LEN=100) :: obj_name

    !<<<<<<<<<<<<<<<<<<<<<<< Start of routine code >>>>>>>>>>>>>>>>>>

    CALL DebugMessage( " --- At level: " &
                    // levelstr(1:LengString(levelstr)) )

    CALL h5gn_members_f(p_id,levelstr,nmembers,ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining level members"
       return
    ENDIF
    CALL DebugMessage("     > No. members: ",nmembers)

    ltype=-1
    FOUND=.false.
    
  
    CALL DebugMessage( " --- Determining type of next level: " &
                    // sublevelstr(1:LengString(sublevelstr)) )

    ! According to the manual, these can be the object types:
    !   H5G_LINK_F
    !   H5G_GROUP_F   --> appears to have number 1
    !   H5G_DATASET_F --> appears to have number 2
    !   H5G_TYPE_F
    ! Update 2011/01/24: found that above values depend on system or version. So use
    ! variable names rather than absolute numbers.
    ! Do not know what the first and last are ... do not support them.
    ! Use 'ltype' to indicate the level type:
    !   ltype = 1 : group
    !           2 : dataset
    !           3 : attribute (in fact: not 'group' or 'dataset')
    ! and use these to open the level for the final reading and
    ! for closing the levels.
    !
    ! Note that once a level is open, one could use
    !    CALL h5iget_type_f(d_id,ltype,ErrorFlag)
    ! to find the type "hanging" on d_id, but that routine gives
    ! other numbers than the types stored here, coming from the
    ! h5gget_obj_info_idx_f routine; e.g.: 'dataset' is 2 from
    ! h5gge... and 5 from h5iget...
    !
    ! So let us stick to the numbers defined here throughout
    ! this set of routines.

    DO n=1,nmembers
      
       CALL h5gget_obj_info_idx_f(p_id,levelstr,n-1, &
                  obj_name,obj_type,ErrorFlag)
       IF (ErrorFlag.lt.0) THEN
          ErrorMessage=" *** Error reading member list"
          return
       ENDIF
       CALL h5iget_type_f(p_id,ltype,ErrorFlag)
       
       IF (obj_name == sublevelstr) THEN
      
          if ((obj_type.lt.H5G_GROUP_F).or.(obj_type.gt.H5G_DATASET_F)) then
             ErrorMessage=" *** Unknown level type found for " &
                         // sublevelstr(1:LengString(sublevelstr))
             ErrorFlag=-1
             return
          endif
          ltype=obj_type
          FOUND=.true.
          goto 10  !  do not need to go through the rest of the list
       ENDIF
    ENDDO

 10 if (.not.FOUND) ltype=3      !  attribute = only possibility now

    CALL DebugMessage("     > Type number: ",ltype)

  END SUBROUTINE CheckMembers

! *******************************************************************

  ! Determine a length of a text in a given string by going backwards
  ! through the defined length searching for the last character that
  ! is not a space or a NULL, tab, newline, etc. character.
  !
  ! The 'len' command gives the defined length of the string.

  INTEGER FUNCTION LengString(str)
    !
    !<<<<<<<<<<<<<<<<<<<<<<< Implicit statement >>>>>>>>>>>>>>>>>>>>>
    IMPLICIT NONE

    !<<<<<<<<<<<<<<<<<<<<<<< Arguments >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !Coming in
    CHARACTER (len = *), INTENT(in) :: str

    !Going out:
    !result of the function

    !<<<<<<<<<<<<<<<<<<<<<<< Local variables >>>>>>>>>>>>>>>>>>>>>>>>
    INTEGER :: n

    !<<<<<<<<<<<<<<<<<<<<<<< Start of routine code >>>>>>>>>>>>>>>>>>

    do n=len(str),1,-1
       if ( (str(n:n).ne." ")            &
            .and. (str(n:n).ne.char(0))  &    ! NULL character
            .and. (str(n:n).ne.char(9))  &    ! horizontal tab 
            .and. (str(n:n).ne.char(10)) &    ! newline
            .and. (str(n:n).ne.char(12)) &    ! formfeed
            .and. (str(n:n).ne.char(13)) &    ! carriage return
          ) then
          LengString=n
          return
       endif
    enddo

  END FUNCTION LengString

! *******************************************************************

  ! Print a debug message, if needed.

  SUBROUTINE DebugMessage (DebugStr,DebugValue)

    !<<<<<<<<<<<<<<<<<<<<<<< Implicit statement >>>>>>>>>>>>>>>>>>>>>
    IMPLICIT NONE

    !<<<<<<<<<<<<<<<<<<<<<<< Arguments >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !Coming in
    CHARACTER (len = *), INTENT(in) :: DebugStr
    INTEGER, OPTIONAL,   INTENT(in) :: DebugValue

    !Going out:
    !none

    !<<<<<<<<<<<<<<<<<<<<<<< Local variables >>>>>>>>>>>>>>>>>>>>>>>>
    !none

    !<<<<<<<<<<<<<<<<<<<<<<< Start of routine code >>>>>>>>>>>>>>>>>>

    IF (DebugMsg) THEN
       IF (PRESENT(DebugValue)) THEN
          write(6,'(a,1x,i15)')DebugStr,DebugValue
       ELSE
          write(6,'(a)')DebugStr
       ENDIF
    ENDIF

  END SUBROUTINE DebugMessage




  SUBROUTINE H5ReadI2nteger0D     &
       (                         &
       H5filename,               & !Coming in
       datasetname,              & !Coming in
       dataset                   & !Going out
       )
    !
    !<<<<<<<<<<<<<<<<<<<<<<< Implicit statement >>>>>>>>>>>>>>>>>>>>>
    IMPLICIT NONE

    !<<<<<<<<<<<<<<<<<<<<<<< Arguments >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !Coming in
    CHARACTER (len = *), INTENT(in)                      :: H5filename
    CHARACTER (len = *), INTENT(in)                      :: datasetname

    !Going out:
    INTEGER  ( kind = 2)                                            :: dataset

    !<<<<<<<<<<<<<<<<<<<<<<< Local variables >>>>>>>>>>>>>>>>>>>>>>>>
    INTEGER, PARAMETER                                   :: maxdims = 4
  
    INTEGER(hsize_t), DIMENSION(1)                       :: dims
    INTEGER                                              :: H5dataset
    INTEGER                                              :: ltype

    !<<<<<<<<<<<<<<<<<<<<<<< Start of routine code >>>>>>>>>>>>>>>>>>

    CALL DebugMessage(" === in H5ReadInteger0D ===")

    ! Initialise the HDF5 file and open all levels up to
    ! the one that needs to be read.

    CALL H5Read_init ( H5filename, datasetname )
    IF (ErrorFlag.lt.0) return

    ! Get details of the dataset

    CALL DebugMessage(" --- Determining details of the data")

    ltype=OpenLevels_type(NLevels)

    dims=1

    ! Read the dataset and transfer the result

    CALL DebugMessage(" --- Reading the data")
    IF ( ltype == H5G_DATASET_F) THEN
       CALL h5dread_f(d_id, H5T_NATIVE_INTEGER, H5dataset, dims, ErrorFlag)
    ELSE  !  ltype.eqv.3
       CALL h5aread_f(d_id, H5T_NATIVE_INTEGER, H5dataset, dims, ErrorFlag)
    ENDIF
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error reading data"
       return
    ENDIF

    dataset = H5dataset

    ! Close all that is open

    CALL H5Read_close

  END SUBROUTINE H5ReadI2nteger0D

! *******************************************************************
! *** Reading 0-, 1-, 2- and 3-D Integer Datasets or Attributes   ***
! *******************************************************************


  SUBROUTINE H5ReadInteger0D     &
       (                         &
       H5filename,               & !Coming in
       datasetname,              & !Coming in
       dataset                   & !Going out
       )
    !
    !<<<<<<<<<<<<<<<<<<<<<<< Implicit statement >>>>>>>>>>>>>>>>>>>>>
    IMPLICIT NONE

    !<<<<<<<<<<<<<<<<<<<<<<< Arguments >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !Coming in
    CHARACTER (len = *), INTENT(in)                      :: H5filename
    CHARACTER (len = *), INTENT(in)                      :: datasetname

    !Going out:
    INTEGER                                              :: dataset

    !<<<<<<<<<<<<<<<<<<<<<<< Local variables >>>>>>>>>>>>>>>>>>>>>>>>
    INTEGER, PARAMETER                                   :: maxdims = 4
  
    INTEGER(hsize_t), DIMENSION(1)                       :: dims
    INTEGER                                              :: H5dataset
    INTEGER                                              :: ltype

    !<<<<<<<<<<<<<<<<<<<<<<< Start of routine code >>>>>>>>>>>>>>>>>>

    CALL DebugMessage(" === in H5ReadInteger0D ===")

    ! Initialise the HDF5 file and open all levels up to
    ! the one that needs to be read.

    CALL H5Read_init ( H5filename, datasetname )
    IF (ErrorFlag.lt.0) return

    ! Get details of the dataset

    CALL DebugMessage(" --- Determining details of the data")

    ltype=OpenLevels_type(NLevels)

    dims=1

    ! Read the dataset and transfer the result

    CALL DebugMessage(" --- Reading the data")
    IF ( ltype == H5G_DATASET_F) THEN
       CALL h5dread_f(d_id, H5T_NATIVE_INTEGER, H5dataset, dims, ErrorFlag)
    ELSE  !  ltype.eqv.3
       CALL h5aread_f(d_id, H5T_NATIVE_INTEGER, H5dataset, dims, ErrorFlag)
    ENDIF
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error reading data"
       return
    ENDIF

    dataset = H5dataset

    ! Close all that is open

    CALL H5Read_close

  END SUBROUTINE H5ReadInteger0D

! *******************************************************************

  SUBROUTINE H5ReadInteger1D     &
       (                         &
       H5filename,               & !Coming in
       datasetname,              & !Coming in
       dataset                   & !Going out
       )
    !
    !<<<<<<<<<<<<<<<<<<<<<<< Implicit statement >>>>>>>>>>>>>>>>>>>>>
    IMPLICIT NONE

    !<<<<<<<<<<<<<<<<<<<<<<< Arguments >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !Coming in
    CHARACTER (len = *), INTENT(in)                      :: H5filename
    CHARACTER (len = *), INTENT(in)                      :: datasetname

    !Going out:
    INTEGER, DIMENSION(:), POINTER                       :: dataset

    !<<<<<<<<<<<<<<<<<<<<<<< Local variables >>>>>>>>>>>>>>>>>>>>>>>>
    INTEGER, PARAMETER                                   :: maxdims = 4
    INTEGER                                              :: ndims
    INTEGER(hsize_t), DIMENSION(1)                       :: datadims
    INTEGER(hsize_t), DIMENSION(1)                       :: maxdatadims
    INTEGER(hsize_t), DIMENSION(1)                       :: dims
    INTEGER, DIMENSION(:), ALLOCATABLE, TARGET           :: H5dataset
    INTEGER                                              :: ltype

    !<<<<<<<<<<<<<<<<<<<<<<< Start of routine code >>>>>>>>>>>>>>>>>>

    CALL DebugMessage(" === in H5ReadInteger1D ===")

    ! Initialise the HDF5 file and open all levels up to
    ! the one that needs to be read.

    CALL H5Read_init ( H5filename, datasetname )
    IF (ErrorFlag.lt.0) return

    ! Get details of the dataset

    CALL DebugMessage(" --- Determining details of the data")

    ltype=OpenLevels_type(NLevels)

    IF ( ltype == H5G_DATASET_F) THEN
       CALL h5dget_space_f(d_id,dspace,ErrorFlag)
    ELSE  !  ltype.eqv.3
       CALL h5aget_space_f(d_id,dspace,ErrorFlag)
    ENDIF
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace"
       return
    ENDIF

    CALL h5sget_simple_extent_ndims_f(dspace,ndims,ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace dimensionality"
       return
    ENDIF
    CALL DebugMessage("     > No. of dim.: ",ndims) 

    CALL h5sget_simple_extent_dims_f(dspace,datadims,maxdatadims,ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace size"
       return
    ENDIF

    dims(1) = datadims(1)
    CALL DebugMessage("     > Size dim. 1: ",int(dims(1)))

    ! Read the dataset and transfer the result

    ALLOCATE(H5dataset(dims(1)),stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error allocating H5dataset"
       return
    ENDIF

    CALL DebugMessage(" --- Reading the data")
    IF ( ltype == H5G_DATASET_F) THEN
       CALL h5dread_f(d_id, H5T_NATIVE_INTEGER, H5dataset, dims, ErrorFlag)
    ELSE  !  ltype.eqv.3
       CALL h5aread_f(d_id, H5T_NATIVE_INTEGER, H5dataset, dims, ErrorFlag)
    ENDIF
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error reading data"
       return
    ENDIF

    ALLOCATE(dataset(dims(1)),stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error allocating dataset"
       return
    ENDIF

    dataset = H5dataset

    DEALLOCATE(H5dataset,stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error deallocating H5dataset"
       return
    ENDIF

    ! Close all that is open

    CALL H5Read_close

  END SUBROUTINE H5ReadInteger1D
! *******************************************************************

  SUBROUTINE H5ReadI8nteger1D     &
       (                         &
       H5filename,               & !Coming in
       datasetname,              & !Coming in
       dataset                   & !Going out
       )
    
    use ISO_C_BINDING
    !
    !<<<<<<<<<<<<<<<<<<<<<<< Implicit statement >>>>>>>>>>>>>>>>>>>>>
    IMPLICIT NONE

    !<<<<<<<<<<<<<<<<<<<<<<< Arguments >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !Coming in
    CHARACTER (len = *), INTENT(in)                      :: H5filename
    CHARACTER (len = *), INTENT(in)                      :: datasetname

    !Going out:
    INTEGER ( kind = int8), DIMENSION(:), POINTER                       :: dataset

    !<<<<<<<<<<<<<<<<<<<<<<< Local variables >>>>>>>>>>>>>>>>>>>>>>>>
    INTEGER, PARAMETER                                   :: maxdims = 4
    INTEGER                                              :: ndims
    INTEGER(hsize_t), DIMENSION(1)                       :: datadims
    INTEGER(hsize_t), DIMENSION(1)                       :: maxdatadims
    INTEGER(hsize_t), DIMENSION(1)                       :: dims
    INTEGER(kind=int8), DIMENSION(:), ALLOCATABLE, TARGET           :: H5dataset
    INTEGER                                              :: ltype
    INTEGER(HID_T) :: Nat_Type_Id
    INTEGER(HID_T) :: Dtype_Id
    TYPE(C_PTR) :: f_ptr
    !<<<<<<<<<<<<<<<<<<<<<<< Start of routine code >>>>>>>>>>>>>>>>>>

    CALL DebugMessage(" === in H5ReadInteger1D int8 ===")
      print*,'could be wrong in readhdf5 check.... stop here'
      !stop 
    ! Initialise the HDF5 file and open all levels up to
    ! the one that needs to be read.

    CALL H5Read_init ( H5filename, datasetname )
    IF (ErrorFlag.lt.0) return

    ! Get details of the dataset

    CALL DebugMessage(" --- Determining details of the data")
    
     CALL H5Dget_type_f(D_Id,Dtype_Id,ErrorFlag)
     CALL H5Tget_native_type_f(Dtype_Id,1,Nat_Type_Id,ErrorFlag)
    
    ltype=OpenLevels_type(NLevels)
     
    IF ( ltype == H5G_DATASET_F) THEN
       CALL h5dget_space_f(d_id,dspace,ErrorFlag)
    ELSE  !  ltype.eqv.3
       CALL h5aget_space_f(d_id,dspace,ErrorFlag)
    ENDIF
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace"
       return
    ENDIF

    CALL h5sget_simple_extent_ndims_f(dspace,ndims,ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace dimensionality"
       return
    ENDIF
    CALL DebugMessage("     > No. of dim.: ",ndims) 

    CALL h5sget_simple_extent_dims_f(dspace,datadims,maxdatadims,ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace size"
       return
    ENDIF

    dims(1) = datadims(1)
    CALL DebugMessage("     > Size dim. 1: ",int(dims(1)))

    ! Read the dataset and transfer the result

    ALLOCATE(H5dataset(dims(1)),stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error allocating H5dataset"
       return
    ENDIF

    CALL DebugMessage(" --- Reading the data")
    f_ptr = C_LOC(H5dataset(1))
    IF ( ltype == H5G_DATASET_F) THEN
      ! CALL h5dread_f(d_id, NAT_TYPE_ID, H5dataset , dims, ErrorFlag)
    ELSE  !  ltype.eqv.3
     !  CALL h5aread_f(d_id, H5T_NATIVE_INTEGER, H5dataset, dims, ErrorFlag)
    ENDIF
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error reading data"
       return
    ENDIF

    ALLOCATE(dataset(dims(1)),stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error allocating dataset"
       return
    ENDIF

    dataset = H5dataset

    DEALLOCATE(H5dataset,stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error deallocating H5dataset"
       return
    ENDIF

    ! Close all that is open

    CALL H5Read_close

  END SUBROUTINE H5ReadI8nteger1D
! *******************************************************************

  SUBROUTINE H5ReadInteger2D     &
       (                         &
       H5filename,               & !Coming in
       datasetname,              & !Coming in
       dataset                   & !Going out
       )
    !
    !<<<<<<<<<<<<<<<<<<<<<<< Implicit statement >>>>>>>>>>>>>>>>>>>>>
    IMPLICIT NONE

    !<<<<<<<<<<<<<<<<<<<<<<< Arguments >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !Coming in
    CHARACTER (len = *), INTENT(in)                      :: H5filename
    CHARACTER (len = *), INTENT(in)                      :: datasetname

    !Going out:
    INTEGER,  DIMENSION(:,:), POINTER                    :: dataset

    !<<<<<<<<<<<<<<<<<<<<<<< Local variables >>>>>>>>>>>>>>>>>>>>>>>>
    INTEGER, PARAMETER                                   :: maxdims = 4
    INTEGER                                              :: ndims
    INTEGER(hsize_t), DIMENSION(maxdims)                 :: datadims
    INTEGER(hsize_t), DIMENSION(maxdims)                 :: maxdatadims
    INTEGER(hsize_t), DIMENSION(:), ALLOCATABLE          :: dims
    INTEGER, DIMENSION(:,:), ALLOCATABLE, TARGET         :: H5dataset
    INTEGER                                              :: ltype

    !<<<<<<<<<<<<<<<<<<<<<<< Start of routine code >>>>>>>>>>>>>>>>>>

    CALL DebugMessage(" === in H5ReadInteger2D ===")

    ! Initialise the HDF5 file and open all levels up to
    ! the one that needs to be read.

    CALL H5Read_init ( H5filename, datasetname )
    IF (ErrorFlag.lt.0) return

    ! Get details of the dataset

    CALL DebugMessage(" --- Determining details of the data")

    ltype=OpenLevels_type(NLevels)

    IF ( ltype == H5G_DATASET_F) THEN
       CALL h5dget_space_f(d_id,dspace,ErrorFlag)
    ELSE  !  ltype.eqv.3
       CALL h5aget_space_f(d_id,dspace,ErrorFlag)
    ENDIF
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace"
       return
    ENDIF

    CALL h5sget_simple_extent_ndims_f(dspace,ndims,ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace dimensionality"
       return
    ENDIF
    CALL DebugMessage("     > No. of dim.: ",ndims) 

    ALLOCATE(dims(ndims),stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error allocating dims"
       return
    ENDIF

    CALL h5sget_simple_extent_dims_f(dspace,datadims,maxdatadims,ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace size"
       return
    ENDIF

    dims(1) = datadims(1)
    dims(2) = datadims(2)
    CALL DebugMessage("     > Size dim. 1: ",int(dims(1)))
    CALL DebugMessage("     > Size dim. 2: ",int(dims(2)))

    ! Read the dataset and transfer the result

    ALLOCATE(H5dataset(dims(1),dims(2)),stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error allocating H5dataset"
       return
    ENDIF

    CALL DebugMessage(" --- Reading the data")
    IF ( ltype == H5G_DATASET_F) THEN
       CALL h5dread_f(d_id, H5T_NATIVE_INTEGER, H5dataset, dims, ErrorFlag)
    ELSE  !  ltype.eqv.3
       CALL h5aread_f(d_id, H5T_NATIVE_INTEGER, H5dataset, dims, ErrorFlag)
    ENDIF
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error reading data"
       return
    ENDIF

    ALLOCATE(dataset(dims(1),dims(2)),stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error allocating dataset"
       return
    ENDIF

    dataset = H5dataset

    DEALLOCATE(H5dataset)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error deallocating H5dataset"
       return
    ENDIF

    DEALLOCATE(dims)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error deallocating dims"
       return
    ENDIF

    ! Close all that is open

    CALL H5Read_close

  END SUBROUTINE H5ReadInteger2D


! *******************************************************************

  SUBROUTINE H5Read2Integer2D     &
       (                         &
       H5filename,               & !Coming in
       datasetname,              & !Coming in
       dataset                   & !Going out
       )
    !
    !<<<<<<<<<<<<<<<<<<<<<<< Implicit statement >>>>>>>>>>>>>>>>>>>>>
    IMPLICIT NONE

    !<<<<<<<<<<<<<<<<<<<<<<< Arguments >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !Coming in
    CHARACTER (len = *), INTENT(in)                      :: H5filename
    CHARACTER (len = *), INTENT(in)                      :: datasetname

    !Going out:
    INTEGER(kind = 2),  DIMENSION(:,:), POINTER                    :: dataset

    !<<<<<<<<<<<<<<<<<<<<<<< Local variables >>>>>>>>>>>>>>>>>>>>>>>>
    INTEGER, PARAMETER                                   :: maxdims = 4
    INTEGER                                              :: ndims
    INTEGER(hsize_t), DIMENSION(maxdims)                 :: datadims
    INTEGER(hsize_t), DIMENSION(maxdims)                 :: maxdatadims
    INTEGER(hsize_t), DIMENSION(:), ALLOCATABLE          :: dims
    INTEGER, DIMENSION(:,:), ALLOCATABLE, TARGET         :: H5dataset
    INTEGER                                              :: ltype

    !<<<<<<<<<<<<<<<<<<<<<<< Start of routine code >>>>>>>>>>>>>>>>>>

    CALL DebugMessage(" === in H5Read2Integer2D ===")

    ! Initialise the HDF5 file and open all levels up to
    ! the one that needs to be read.

    CALL H5Read_init ( H5filename, datasetname )
    IF (ErrorFlag.lt.0) return

    ! Get details of the dataset

    CALL DebugMessage(" --- Determining details of the data")

    ltype=OpenLevels_type(NLevels)

    IF ( ltype == H5G_DATASET_F) THEN
       CALL h5dget_space_f(d_id,dspace,ErrorFlag)
    ELSE  !  ltype.eqv.3
       CALL h5aget_space_f(d_id,dspace,ErrorFlag)
    ENDIF
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace"
       return
    ENDIF

    CALL h5sget_simple_extent_ndims_f(dspace,ndims,ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace dimensionality"
       return
    ENDIF
    CALL DebugMessage("     > No. of dim.: ",ndims) 

    ALLOCATE(dims(ndims),stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error allocating dims"
       return
    ENDIF

    CALL h5sget_simple_extent_dims_f(dspace,datadims,maxdatadims,ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace size"
       return
    ENDIF

    dims(1) = datadims(1)
    dims(2) = datadims(2)
    CALL DebugMessage("     > Size dim. 1: ",int(dims(1)))
    CALL DebugMessage("     > Size dim. 2: ",int(dims(2)))

    ! Read the dataset and transfer the result

    ALLOCATE(H5dataset(dims(1),dims(2)),stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error allocating H5dataset"
       return
    ENDIF

    CALL DebugMessage(" --- Reading the data")
    IF ( ltype == H5G_DATASET_F) THEN
       CALL h5dread_f(d_id, H5T_NATIVE_INTEGER, H5dataset, dims, ErrorFlag)
    ELSE  !  ltype.eqv.3
       CALL h5aread_f(d_id, H5T_NATIVE_INTEGER, H5dataset, dims, ErrorFlag)
    ENDIF
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error reading data"
       return
    ENDIF

    ALLOCATE(dataset(dims(1),dims(2)),stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error allocating dataset"
       return
    ENDIF

    dataset = H5dataset

    DEALLOCATE(H5dataset)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error deallocating H5dataset"
       return
    ENDIF

    DEALLOCATE(dims)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error deallocating dims"
       return
    ENDIF

    ! Close all that is open

    CALL H5Read_close

  END SUBROUTINE H5Read2Integer2D


! *******************************************************************

  SUBROUTINE H5ReadInteger2DSub     &
       (                         &
       H5filename,               & !Coming in
       datasetname,              & !Coming in
       offset_in, &
       count_in, &
       dataset                   & !Going out
       )
    !
    !<<<<<<<<<<<<<<<<<<<<<<< Implicit statement >>>>>>>>>>>>>>>>>>>>>
    IMPLICIT NONE

    !<<<<<<<<<<<<<<<<<<<<<<< Arguments >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !Coming in
    CHARACTER (len = *), INTENT(in)                      :: H5filename
    CHARACTER (len = *), INTENT(in)                      :: datasetname
    
    INTEGER , DIMENSION(2) , INTENT(in) :: offset_in, count_in   

    !Going out:
    INTEGER,  DIMENSION(:,:), POINTER                    :: dataset

    !<<<<<<<<<<<<<<<<<<<<<<< Local variables >>>>>>>>>>>>>>>>>>>>>>>>
    INTEGER, PARAMETER                                   :: maxdims = 4
    INTEGER                                              :: ndims
    INTEGER(hsize_t), DIMENSION(maxdims)                 :: datadims
    INTEGER(hsize_t), DIMENSION(maxdims)                 :: maxdatadims
    INTEGER(hsize_t), DIMENSION(:), ALLOCATABLE          :: dims
    INTEGER, DIMENSION(:,:), ALLOCATABLE, TARGET         :: H5dataset
   
    
    integer ( hsize_t ), dimension(2) :: offset , count
    integer ( hid_t ) :: memspace
    integer , parameter :: memrank = 2
    INTEGER(HSIZE_T), DIMENSION(2) :: offset_out       !offset writing
     INTEGER(HSIZE_T), DIMENSION(2) :: count_out       !offset writing 

    !<<<<<<<<<<<<<<<<<<<<<<< Start of routine code >>>>>>>>>>>>>>>>>>

    CALL DebugMessage(" === in H5ReadInteger2DSub ===")
    
    offset = offset_in
    count = count_in
    
    offset_out = [0 , 0]
    count_out = count
    

    ! Initialise the HDF5 file and open all levels up to
    ! the one that needs to be read.
   
    CALL H5Read_init ( H5filename, datasetname )
    IF (ErrorFlag.lt.0) return

    ! Get details of the dataset

    CALL DebugMessage(" --- Determining details of the data")
  
    CALL h5dget_space_f(d_id,dspace,ErrorFlag)
    
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace"
       return
    ENDIF

    CALL h5sget_simple_extent_ndims_f(dspace,ndims,ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace dimensionality"
       return
    ENDIF
    CALL DebugMessage("     > No. of dim.: ",ndims) 

    ALLOCATE(dims(ndims),stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error allocating dims"
       return
    ENDIF

    CALL h5sget_simple_extent_dims_f(dspace,datadims,maxdatadims,ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace size"
       return
    ENDIF

    dims(1) = count(1)
    dims(2) = count(2)
    CALL DebugMessage("     > Size dim. 1: ",int(dims(1)))
    CALL DebugMessage("     > Size dim. 2: ",int(dims(2)))

    ! Read the dataset and transfer the result

    ALLOCATE(H5dataset(count(1),count(2)),stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error allocating H5dataset"
       return
    ENDIF
 !
    ! Select hyperslab in the dataset.
    !
    CALL h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, &
                                offset , count , ErrorFlag)
                                
     !
     ! Create memory dataspace.
     !
     
      CALL h5screate_simple_f(memrank, dims , memspace, errorFlag)

     !
     ! Select hyperslab in memory.
     !
      CALL h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
                                offset_out, count_out, errorFlag)               
                                            
    CALL DebugMessage(" --- Reading the data")
  
    CALL h5dread_f(d_id, H5T_NATIVE_INTEGER, H5dataset, count, ErrorFlag,memspace,dspace)
   
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error reading data"
       return
    ENDIF
 
    ALLOCATE(dataset(count(1),count(2)),stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error allocating dataset"
       return
    ENDIF

    dataset = H5dataset

    DEALLOCATE(H5dataset,stat = AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error deallocating H5dataset"
       return
    ENDIF

    DEALLOCATE(dims,stat = AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error deallocating dims"
       return
    ENDIF

    ! Close all that is open

    CALL H5Read_close

  END SUBROUTINE H5ReadInteger2DSub

! *******************************************************************
! *******************************************************************

  SUBROUTINE H5Read2Integer2DSub     &
       (                         &
       H5filename,               & !Coming in
       datasetname,              & !Coming in
       offset_in, &
       count_in, &
       dataset                   & !Going out
       )
    !
    !<<<<<<<<<<<<<<<<<<<<<<< Implicit statement >>>>>>>>>>>>>>>>>>>>>
    IMPLICIT NONE

    !<<<<<<<<<<<<<<<<<<<<<<< Arguments >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !Coming in
    CHARACTER (len = *), INTENT(in)                      :: H5filename
    CHARACTER (len = *), INTENT(in)                      :: datasetname
    
    INTEGER , DIMENSION(2) , INTENT(in) :: offset_in, count_in   

    !Going out:
    INTEGER(kind = 2),  DIMENSION(:,:), POINTER                    :: dataset 

    !<<<<<<<<<<<<<<<<<<<<<<< Local variables >>>>>>>>>>>>>>>>>>>>>>>>
    INTEGER, PARAMETER                                   :: maxdims = 4
    INTEGER                                              :: ndims
    INTEGER(hsize_t), DIMENSION(maxdims)                 :: datadims
    INTEGER(hsize_t), DIMENSION(maxdims)                 :: maxdatadims
    INTEGER(hsize_t), DIMENSION(:), ALLOCATABLE          :: dims
    INTEGER, DIMENSION(:,:), ALLOCATABLE, TARGET         :: H5dataset
   
    
    integer ( hsize_t ), dimension(2) :: offset , count
    integer ( hid_t ) :: memspace
    integer , parameter :: memrank = 2
    INTEGER(HSIZE_T), DIMENSION(2) :: offset_out       !offset writing
    INTEGER(HSIZE_T), DIMENSION(2) :: count_out       !offset writing 

    !<<<<<<<<<<<<<<<<<<<<<<< Start of routine code >>>>>>>>>>>>>>>>>>

    CALL DebugMessage(" === in H5Read2Integer2DSub ===")
    
    offset = offset_in
    count = count_in
    
    offset_out = [0 , 0]
    count_out = count
    

    ! Initialise the HDF5 file and open all levels up to
    ! the one that needs to be read.
   
    CALL H5Read_init ( H5filename, datasetname )
    IF (ErrorFlag.lt.0) return

    ! Get details of the dataset

    CALL DebugMessage(" --- Determining details of the data")
  
    CALL h5dget_space_f(d_id,dspace,ErrorFlag)
    
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace"
       return
    ENDIF

    CALL h5sget_simple_extent_ndims_f(dspace,ndims,ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace dimensionality"
       return
    ENDIF
    CALL DebugMessage("     > No. of dim.: ",ndims) 

    ALLOCATE(dims(ndims),stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error allocating dims"
       return
    ENDIF

    CALL h5sget_simple_extent_dims_f(dspace,datadims,maxdatadims,ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace size"
       return
    ENDIF
    
   
    
    dims(1) = count(1)
    dims(2) = count(2)
    CALL DebugMessage("     > Size dim. 1: ",int(dims(1)))
    CALL DebugMessage("     > Size dim. 2: ",int(dims(2)))

    ! Read the dataset and transfer the result

    ALLOCATE(H5dataset(count(1),count(2)),stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error allocating H5dataset"
       return
    ENDIF
 !
    ! Select hyperslab in the dataset.
    !
    CALL h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, &
                                offset , count , ErrorFlag)
                                
    !
    ! Create memory dataspace.
    !
     
      CALL h5screate_simple_f(memrank, dims , memspace, errorFlag)

    !
    ! Select hyperslab in memory.
    !
      CALL h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
                                offset_out, count_out, errorFlag)               
                                            
    CALL DebugMessage(" --- Reading the data")
  
    CALL h5dread_f(d_id, H5T_NATIVE_INTEGER, H5dataset, count, ErrorFlag,memspace,dspace)
   
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error reading data"
       return
    ENDIF
 
    ALLOCATE(dataset(count(1),count(2)),stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error allocating dataset"
       return
    ENDIF
    
   
    dataset = H5dataset
     
    DEALLOCATE(H5dataset,stat = AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error deallocating H5dataset"
       print *, errormessage
       return
    ENDIF

    DEALLOCATE(dims,stat = AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error deallocating dims"
       return
    ENDIF

    ! Close all that is open
    
    call h5sclose_f(memspace,errorFlag)
    call h5sclose_f(dspace,errorFlag)
    CALL H5Read_close

  END SUBROUTINE H5Read2Integer2DSub

! *******************************************************************


  SUBROUTINE H5ReadInteger3D     &
       (                         &
       H5filename,               & !Coming in
       datasetname,              & !Coming in
       dataset                   & !Going out
       )
    !
    !<<<<<<<<<<<<<<<<<<<<<<< Implicit statement >>>>>>>>>>>>>>>>>>>>>
    IMPLICIT NONE

    !<<<<<<<<<<<<<<<<<<<<<<< Arguments >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !Coming in
    CHARACTER (len = *), INTENT(in)                      :: H5filename
    CHARACTER (len = *), INTENT(in)                      :: datasetname

    !Going out:
    INTEGER,  DIMENSION(:,:,:), POINTER                  :: dataset

    !<<<<<<<<<<<<<<<<<<<<<<< Local variables >>>>>>>>>>>>>>>>>>>>>>>>
    INTEGER, PARAMETER                                   :: maxdims = 4
    INTEGER                                              :: ndims
    INTEGER(hsize_t), DIMENSION(maxdims)                 :: datadims
    INTEGER(hsize_t), DIMENSION(maxdims)                 :: maxdatadims
    INTEGER(hsize_t), DIMENSION(:), ALLOCATABLE          :: dims
    INTEGER, DIMENSION(:,:,:), ALLOCATABLE, TARGET       :: H5dataset
    INTEGER                                              :: ltype

    !<<<<<<<<<<<<<<<<<<<<<<< Start of routine code >>>>>>>>>>>>>>>>>>

    CALL DebugMessage(" === in H5ReadInteger3D ===")

    ! Initialise the HDF5 file and open all levels up to
    ! the one that needs to be read.

    CALL H5Read_init ( H5filename, datasetname )
    IF (ErrorFlag.lt.0) return

    ! Get details of the dataset

    CALL DebugMessage(" --- Determining details of the data")

    ltype=OpenLevels_type(NLevels)

    IF ( ltype == H5G_DATASET_F) THEN
       CALL h5dget_space_f(d_id,dspace,ErrorFlag)
    ELSE  !  ltype.eqv.3
       CALL h5aget_space_f(d_id,dspace,ErrorFlag)
    ENDIF
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace"
       return
    ENDIF

    CALL h5sget_simple_extent_ndims_f(dspace,ndims,ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace dimensionality"
       return
    ENDIF
    CALL DebugMessage("     > No. of dim.: ",ndims) 

    ALLOCATE(dims(ndims),stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error allocating dims"
       return
    ENDIF

    CALL h5sget_simple_extent_dims_f(dspace,datadims,maxdatadims,ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace size"
       return
    ENDIF

    dims(1) = datadims(1)
    dims(2) = datadims(2)
    dims(3) = datadims(3)
    CALL DebugMessage("     > Size dim. 1: ",int(dims(1)))
    CALL DebugMessage("     > Size dim. 2: ",int(dims(2)))
    CALL DebugMessage("     > Size dim. 3: ",int(dims(3)))

    ! Read the dataset and transfer the result

    ALLOCATE(H5dataset(dims(1),dims(2),dims(3)),stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error allocating H5dataset"
       return
    ENDIF

    CALL DebugMessage(" --- Reading the data")
    IF ( ltype == H5G_DATASET_F) THEN
       CALL h5dread_f(d_id, H5T_NATIVE_INTEGER, H5dataset, dims, ErrorFlag)
    ELSE  !  ltype.eqv.3
       CALL h5aread_f(d_id, H5T_NATIVE_INTEGER, H5dataset, dims, ErrorFlag)
    ENDIF
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error reading data"
       return
    ENDIF

    ALLOCATE(dataset(dims(1),dims(2),dims(3)),stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error allocating dataset"
       return
    ENDIF

    dataset = H5dataset

    DEALLOCATE(H5dataset,stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error deallocating H5dataset"
       return
    ENDIF

    DEALLOCATE(dims,stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error deallocating dims"
       return
    ENDIF

    ! Close all that is open

    CALL H5Read_close

  END SUBROUTINE H5ReadInteger3D


! *******************************************************************
! *** Reading 0-, 1-, 2- and 3-D Real Datasets or Attributes      ***
! *******************************************************************

  SUBROUTINE H5ReadReal0D        &
       (                         &
       H5filename,               & !Coming in
       datasetname,              & !Coming in
       dataset                   & !Going out
       )
    !
    !<<<<<<<<<<<<<<<<<<<<<<< Implicit statement >>>>>>>>>>>>>>>>>>>>>
    IMPLICIT NONE

    !<<<<<<<<<<<<<<<<<<<<<<< Arguments >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !Coming in
    CHARACTER (len = *), INTENT(in)                      :: H5filename
    CHARACTER (len = *), INTENT(in)                      :: datasetname

    !Going out:
    REAL(KIND=SP)                                        :: dataset

    !<<<<<<<<<<<<<<<<<<<<<<< Local variables >>>>>>>>>>>>>>>>>>>>>>>>
    INTEGER, PARAMETER                                   :: maxdims = 4
    
    INTEGER(hsize_t), DIMENSION(1)                       :: dims
    REAL(KIND=SP)                                        :: H5dataset
    INTEGER                                              :: ltype

    !<<<<<<<<<<<<<<<<<<<<<<< Start of routine code >>>>>>>>>>>>>>>>>>

    CALL DebugMessage(" === in H5ReadReal0D ===")

    ! Initialise the HDF5 file and open all levels up to
    ! the one that needs to be read.

    CALL H5Read_init ( H5filename, datasetname )
    IF (ErrorFlag.lt.0) return

    ! Get details of the dataset

    CALL DebugMessage(" --- Determining details of the data")

    ltype=OpenLevels_type(NLevels)

    dims=1

    ! Read the dataset and transfer the result

    CALL DebugMessage(" --- Reading the data")
    IF ( ltype == H5G_DATASET_F) THEN
       CALL h5dread_f(d_id, H5T_NATIVE_REAL, H5dataset, dims, ErrorFlag)
    ELSE  !  ltype.eqv.3
       CALL h5aread_f(d_id, H5T_NATIVE_REAL, H5dataset, dims, ErrorFlag)
    ENDIF
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error reading data"
       return
    ENDIF

    dataset = H5dataset

    ! Close all that is open

    CALL H5Read_close

  END SUBROUTINE H5ReadReal0D

! *******************************************************************

  SUBROUTINE H5ReadReal1D        &
       (                         &
       H5filename,               & !Coming in
       datasetname,              & !Coming in
       dataset                   & !Going out
       )
    !
    !<<<<<<<<<<<<<<<<<<<<<<< Implicit statement >>>>>>>>>>>>>>>>>>>>>
    IMPLICIT NONE

    !<<<<<<<<<<<<<<<<<<<<<<< Arguments >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !Coming in
    CHARACTER (len = *), INTENT(in)                      :: H5filename
    CHARACTER (len = *), INTENT(in)                      :: datasetname

    !Going out:
    REAL(KIND=SP), DIMENSION(:), POINTER                 :: dataset

    !<<<<<<<<<<<<<<<<<<<<<<< Local variables >>>>>>>>>>>>>>>>>>>>>>>>
    INTEGER, PARAMETER                                   :: maxdims = 4
    INTEGER                                              :: ndims
    INTEGER(hsize_t), DIMENSION(1)                       :: datadims
    INTEGER(hsize_t), DIMENSION(1)                       :: maxdatadims
    INTEGER(hsize_t), DIMENSION(1)                       :: dims
    REAL(KIND=SP), DIMENSION(:), ALLOCATABLE, TARGET     :: H5dataset
    INTEGER                                              :: ltype

    !<<<<<<<<<<<<<<<<<<<<<<< Start of routine code >>>>>>>>>>>>>>>>>>

    CALL DebugMessage(" === in H5ReadReal1D ===")

    ! Initialise the HDF5 file and open all levels up to
    ! the one that needs to be read.

    CALL H5Read_init ( H5filename, datasetname )
    IF (ErrorFlag.lt.0) return

    ! Get details of the dataset

    CALL DebugMessage(" --- Determining details of the data")

    ltype=OpenLevels_type(NLevels)

    CALL DebugMessage(" --- Determining details of the data")
    IF ( ltype == H5G_DATASET_F) THEN
       CALL h5dget_space_f(d_id,dspace,ErrorFlag)
    ELSE  !  ltype.eqv.3
       CALL h5aget_space_f(d_id,dspace,ErrorFlag)
    ENDIF
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace"
       return
    ENDIF

    CALL h5sget_simple_extent_ndims_f(dspace,ndims,ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace dimensionality"
       return
    ENDIF
    CALL DebugMessage("     > No. of dim.: ",ndims) 

    CALL h5sget_simple_extent_dims_f(dspace,datadims,maxdatadims,ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace size"
       return
    ENDIF

    dims(1) = datadims(1)
    CALL DebugMessage("     > Size dim. 1: ",int(dims(1)))

    ! Read the dataset and transfer the result

    ALLOCATE(H5dataset(dims(1)),stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error allocating H5dataset"
       return
    ENDIF

    CALL DebugMessage(" --- Reading the data")
    IF ( ltype == H5G_DATASET_F) THEN
       CALL h5dread_f(d_id, H5T_NATIVE_REAL, H5dataset, dims, ErrorFlag)
    ELSE  !  ltype.eqv.3
       CALL h5aread_f(d_id, H5T_NATIVE_REAL, H5dataset, dims, ErrorFlag)
    ENDIF
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error reading data"
       return
    ENDIF

    ALLOCATE(dataset(dims(1)),stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error allocating dataset"
       return
    ENDIF

    dataset = H5dataset

    DEALLOCATE(H5dataset,stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error deallocating H5dataset"
       return
    ENDIF

    ! Close all that is open

    CALL H5Read_close

  END SUBROUTINE H5ReadReal1D

! *******************************************************************

  SUBROUTINE H5ReadReal2D        &
       (                         &
       H5filename,               & !Coming in
       datasetname,              & !Coming in
       dataset                   & !Going out
       )
    !
    !<<<<<<<<<<<<<<<<<<<<<<< Implicit statement >>>>>>>>>>>>>>>>>>>>>
    IMPLICIT NONE

    !<<<<<<<<<<<<<<<<<<<<<<< Arguments >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !Coming in
    CHARACTER (len = *), INTENT(in)                      :: H5filename
    CHARACTER (len = *), INTENT(in)                      :: datasetname

    !Going out:
    REAL(KIND=SP), DIMENSION(:,:), POINTER               :: dataset

    !<<<<<<<<<<<<<<<<<<<<<<< Local variables >>>>>>>>>>>>>>>>>>>>>>>>
    INTEGER, PARAMETER                                   :: maxdims = 4
    INTEGER                                              :: ndims
    INTEGER(hsize_t), DIMENSION(maxdims)                 :: datadims
    INTEGER(hsize_t), DIMENSION(maxdims)                 :: maxdatadims
    INTEGER(hsize_t), DIMENSION(:), ALLOCATABLE          :: dims
    REAL(KIND=SP), DIMENSION(:,:), ALLOCATABLE, TARGET   :: H5dataset
    INTEGER                                              :: ltype

    !<<<<<<<<<<<<<<<<<<<<<<< Start of routine code >>>>>>>>>>>>>>>>>>

    CALL DebugMessage(" === in H5ReadReal2D ===")

    ! Initialise the HDF5 file and open all levels up to
    ! the one that needs to be read.

    CALL H5Read_init ( H5filename, datasetname )
    IF (ErrorFlag.lt.0) return

    ! Get details of the dataset

    CALL DebugMessage(" --- Determining details of the data")

    ltype=OpenLevels_type(NLevels)

    IF ( ltype == H5G_DATASET_F) THEN
       CALL h5dget_space_f(d_id,dspace,ErrorFlag)
    ELSE  !  ltype.eqv.3
       CALL h5aget_space_f(d_id,dspace,ErrorFlag)
    ENDIF
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace"
       return
    ENDIF

    CALL h5sget_simple_extent_ndims_f(dspace,ndims,ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace dimensionality"
       return
    ENDIF
    CALL DebugMessage("     > No. of dim.: ",ndims) 

    ALLOCATE(dims(ndims),stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error allocating dims"
       return
    ENDIF

    CALL h5sget_simple_extent_dims_f(dspace,datadims,maxdatadims,ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace size"
       return
    ENDIF

    dims(1) = datadims(1)
    dims(2) = datadims(2)
    CALL DebugMessage("     > Size dim. 1: ",int(dims(1)))
    CALL DebugMessage("     > Size dim. 2: ",int(dims(2)))

    ! Read the dataset and transfer the result

    ALLOCATE(H5dataset(dims(1),dims(2)),stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error allocating H5dataset"
       return
    ENDIF

    CALL DebugMessage(" --- Reading the data")
    IF ( ltype == H5G_DATASET_F) THEN
       CALL h5dread_f(d_id, H5T_NATIVE_REAL, H5dataset, dims, ErrorFlag)
    ELSE  !  ltype.eqv.3
       CALL h5aread_f(d_id, H5T_NATIVE_REAL, H5dataset, dims, ErrorFlag)
    ENDIF
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error reading data"
       return
    ENDIF

    ALLOCATE(dataset(dims(1),dims(2)),stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error allocating dataset"
       return
    ENDIF
   
    dataset = H5dataset
   
    DEALLOCATE(H5dataset,stat=AllocStat)
   
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error deallocating H5dataset"
       return
    ENDIF

    DEALLOCATE(dims,stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error deallocating dims"
       return
    ENDIF

    ! Close all that is open

    CALL H5Read_close

  END SUBROUTINE H5ReadReal2D

! *******************************************************************

SUBROUTINE H5ReadReal2DSub        &
       (                         &
       H5filename,               & !Coming in
       datasetname,              & !Coming in
       offset_in,                        & !Coming in
       count_in ,                      & !Coming in
       dataset                   & !Going out
       )
    !
    !<<<<<<<<<<<<<<<<<<<<<<< Implicit statement >>>>>>>>>>>>>>>>>>>>>
    IMPLICIT NONE

    !<<<<<<<<<<<<<<<<<<<<<<< Arguments >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !Coming in
    CHARACTER (len = *), INTENT(in)                      :: H5filename
    CHARACTER (len = *), INTENT(in)                      :: datasetname
    INTEGER , INTENT(in)                                        :: offset_in ( 2) 
    INTEGER , INTENT(in)                                        :: count_in ( 2) 

    !Going out:
    REAL(KIND=SP), DIMENSION(:,:), POINTER               :: dataset

    !<<<<<<<<<<<<<<<<<<<<<<< Local variables >>>>>>>>>>>>>>>>>>>>>>>>
    INTEGER, PARAMETER                                   :: maxdims = 4
 
    REAL, DIMENSION(:,:), ALLOCATABLE, TARGET   :: H5dataset
   
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T), DIMENSION(2) :: offset              !offset reading dataset
    INTEGER(HSIZE_T), DIMENSION(2) :: offset_out          !offset writing
    INTEGER(HSIZE_T), DIMENSION(2) :: count               !size of hyperslab
    INTEGER, PARAMETER :: memrank = 2
     INTEGER(HSIZE_T), DIMENSION(2) :: count_out       !offset writing 
    !<<<<<<<<<<<<<<<<<<<<<<< Start of routine code >>>>>>>>>>>>>>>>>>

    CALL DebugMessage(" === in H5ReadReal2DSub ===")
    
    count = count_in
    offset = offset_in
   
    offset_out = [0,0]
    count_out = count
    
    ! Initialise the HDF5 file and open all levels up to
    ! the one that needs to be read.

    CALL H5Read_init ( H5filename, datasetname )
    IF (ErrorFlag.lt.0) return

    ! Get details of the dataset

    CALL DebugMessage(" --- Determining details of the data")

  
    CALL h5dget_space_f(d_id,dspace,ErrorFlag)
   
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace"
       return
    ENDIF

   
    ! - hyperslab
    call h5sselect_hyperslab_f ( dspace , H5S_SELECT_SET_F , offset, count , ErrorFlag )
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error selecting hyperspace"
       return
    ENDIF
   
    call h5screate_simple_f(memrank, count , memspace, ErrorFlag)
    
    call h5sselect_hyperslab_f(memspace,H5S_SELECT_SET_F,offset_out,count_out,ErrorFlag)
    
    
   
    CALL DebugMessage("     > Size dim. 1: ",int(count(1)))
    CALL DebugMessage("     > Size dim. 2: ",int(count(2)))

    ! Read the dataset and transfer the result

    ALLOCATE(H5dataset(count(1),count(2)),stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error allocating H5dataset"
       return
    ENDIF

    CALL DebugMessage(" --- Reading the data ")
    
    CALL h5dread_f(d_id, H5T_NATIVE_REAL, H5dataset, count, ErrorFlag,memspace, dspace)
    
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error reading data"
       return
    ENDIF
   
    ALLOCATE(dataset(count(1),count(2)),stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error allocating dataset"
       return
    ENDIF
  
    dataset = H5dataset
   
    DEALLOCATE(H5dataset,stat=AllocStat)
  
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error deallocating H5dataset"
       return
    ENDIF
 
    CALL H5Read_close

  END SUBROUTINE H5ReadReal2DSub

! *******************************************************************


  SUBROUTINE H5ReadReal3D        &
       (                         &
       H5filename,               & !Coming in
       datasetname,              & !Coming in
       dataset                   & !Going out
       )
    !
    !<<<<<<<<<<<<<<<<<<<<<<< Implicit statement >>>>>>>>>>>>>>>>>>>>>
    IMPLICIT NONE

    !<<<<<<<<<<<<<<<<<<<<<<< Arguments >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !Coming in
    CHARACTER (len = *), INTENT(in)                      :: H5filename
    CHARACTER (len = *), INTENT(in)                      :: datasetname

    !Going out:
    REAL(KIND=SP), DIMENSION(:,:,:), POINTER             :: dataset

    !<<<<<<<<<<<<<<<<<<<<<<< Local variables >>>>>>>>>>>>>>>>>>>>>>>>
    INTEGER, PARAMETER                                   :: maxdims = 4
    INTEGER                                              :: ndims
    INTEGER(hsize_t), DIMENSION(maxdims)                 :: datadims
    INTEGER(hsize_t), DIMENSION(maxdims)                 :: maxdatadims
    INTEGER(hsize_t), DIMENSION(:), ALLOCATABLE          :: dims
    REAL(KIND=SP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: H5dataset
    INTEGER                                              :: ltype

    !<<<<<<<<<<<<<<<<<<<<<<< Start of routine code >>>>>>>>>>>>>>>>>>

    CALL DebugMessage(" === in H5ReadReal3D ===")

    ! Initialise the HDF5 file and open all levels up to
    ! the one that needs to be read.

    CALL H5Read_init ( H5filename, datasetname )
    IF (ErrorFlag.lt.0) return

    ! Get details of the dataset

    CALL DebugMessage(" --- Determining details of the data")

    ltype=OpenLevels_type(NLevels)

    IF ( ltype == H5G_DATASET_F) THEN
       CALL h5dget_space_f(d_id,dspace,ErrorFlag)
    ELSE  !  ltype.eqv.3
       CALL h5aget_space_f(d_id,dspace,ErrorFlag)
    ENDIF
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace"
       return
    ENDIF

    CALL h5sget_simple_extent_ndims_f(dspace,ndims,ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace dimensionality"
       return
    ENDIF
    CALL DebugMessage("     > No. of dim.: ",ndims) 

    ALLOCATE(dims(ndims),stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error allocating dims"
       return
    ENDIF

    CALL h5sget_simple_extent_dims_f(dspace,datadims,maxdatadims,ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace size"
       return
    ENDIF

    dims(1) = datadims(1)
    dims(2) = datadims(2)
    dims(3) = datadims(3)
    CALL DebugMessage("     > Size dim. 1: ",int(dims(1)))
    CALL DebugMessage("     > Size dim. 2: ",int(dims(2)))
    CALL DebugMessage("     > Size dim. 3: ",int(dims(3)))

    ! Read the dataset and transfer the result

    ALLOCATE(H5dataset(dims(1),dims(2),dims(3)),stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error allocating H5dataset"
       return
    ENDIF


    CALL DebugMessage(" --- Reading the data")
    IF ( ltype == H5G_DATASET_F) THEN
       CALL h5dread_f(d_id, H5T_NATIVE_REAL, H5dataset, dims, ErrorFlag)
    ELSE  !  ltype.eqv.3
       CALL h5aread_f(d_id, H5T_NATIVE_REAL, H5dataset, dims, ErrorFlag)
    ENDIF
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error reading data"
       return
    ENDIF

    ALLOCATE(dataset(dims(1),dims(2),dims(3)),stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error allocating dataset"
       return
    ENDIF

    dataset = H5dataset

    DEALLOCATE(H5dataset,stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error deallocating H5dataset"
       return
    ENDIF

    DEALLOCATE(dims,stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error deallocating dims"
       return
    ENDIF

    ! Close all that is open

    CALL H5Read_close

  END SUBROUTINE H5ReadReal3D


! *******************************************************************
! *** Reading 0-, 1-, 2- and 3-D Double Datasets or Attributes    ***
! *******************************************************************

  SUBROUTINE H5ReadDouble0D      &
       (                         &
       H5filename,               & !Coming in
       datasetname,              & !Coming in
       dataset                   & !Going out
       )
    !
    !<<<<<<<<<<<<<<<<<<<<<<< Implicit statement >>>>>>>>>>>>>>>>>>>>>
    IMPLICIT NONE

    !<<<<<<<<<<<<<<<<<<<<<<< Arguments >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !Coming in
    CHARACTER (len = *), INTENT(in)                      :: H5filename
    CHARACTER (len = *), INTENT(in)                      :: datasetname

    !Going out:
    REAL(KIND=DP)                                        :: dataset

    !<<<<<<<<<<<<<<<<<<<<<<< Local variables >>>>>>>>>>>>>>>>>>>>>>>>
    INTEGER, PARAMETER                                   :: maxdims = 4
  
    
   
    INTEGER(hsize_t), DIMENSION(1)                       :: dims
    REAL(KIND=DP)                                        :: H5dataset
    INTEGER                                              :: ltype

    !<<<<<<<<<<<<<<<<<<<<<<< Start of routine code >>>>>>>>>>>>>>>>>>

    CALL DebugMessage(" === in H5ReadDboule0D ===")

    ! Initialise the HDF5 file and open all levels up to
    ! the one that needs to be read.

    CALL H5Read_init ( H5filename, datasetname )
    IF (ErrorFlag.lt.0) return

    ! Get details of the dataset

    CALL DebugMessage(" --- Determining details of the data")

    ltype=OpenLevels_type(NLevels)

    dims=1

    ! Read the dataset and transfer the result

    CALL DebugMessage(" --- Reading the data")

    IF ( ltype == H5G_DATASET_F) THEN
       CALL h5dread_f(d_id, H5T_NATIVE_DOUBLE, H5dataset, dims, ErrorFlag)
    ELSE  !  ltype.eqv.3
       CALL h5aread_f(d_id, H5T_NATIVE_DOUBLE, H5dataset, dims, ErrorFlag)
    ENDIF
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error reading data"
       return
    ENDIF

    dataset = H5dataset

    ! Close all that is open

    CALL H5Read_close

  END SUBROUTINE H5ReadDouble0D

! *******************************************************************

  SUBROUTINE H5ReadDouble1D      &
       (                         &
       H5filename,               & !Coming in
       datasetname,              & !Coming in
       dataset                   & !Going out
       )
    !
    !<<<<<<<<<<<<<<<<<<<<<<< Implicit statement >>>>>>>>>>>>>>>>>>>>>
    IMPLICIT NONE

    !<<<<<<<<<<<<<<<<<<<<<<< Arguments >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !Coming in
    CHARACTER (len = *), INTENT(in)                      :: H5filename
    CHARACTER (len = *), INTENT(in)                      :: datasetname

    !Going out:
    REAL(KIND=DP), DIMENSION(:), POINTER                 :: dataset

    !<<<<<<<<<<<<<<<<<<<<<<< Local variables >>>>>>>>>>>>>>>>>>>>>>>>
    INTEGER, PARAMETER                                   :: maxdims = 4
    INTEGER                                              :: ndims
    INTEGER(hsize_t), DIMENSION(1)                       :: datadims
    INTEGER(hsize_t), DIMENSION(1)                       :: maxdatadims
    INTEGER(hsize_t), DIMENSION(1)                       :: dims
    REAL(KIND=DP), DIMENSION(:), ALLOCATABLE, TARGET     :: H5dataset
    INTEGER                                              :: ltype

    !<<<<<<<<<<<<<<<<<<<<<<< Start of routine code >>>>>>>>>>>>>>>>>>

    CALL DebugMessage(" === in H5ReadDouble1D ===")

    ! Initialise the HDF5 file and open all levels up to
    ! the one that needs to be read.

    CALL H5Read_init ( H5filename, datasetname )
    IF (ErrorFlag.lt.0) return

    ! Get details of the dataset

    CALL DebugMessage(" --- Determining details of the data")

    ltype=OpenLevels_type(NLevels)

    IF ( ltype == H5G_DATASET_F) THEN
       CALL h5dget_space_f(d_id,dspace,ErrorFlag)
    ELSE  !  ltype.eqv.3
       CALL h5aget_space_f(d_id,dspace,ErrorFlag)
    ENDIF
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace"
       return
    ENDIF

    CALL h5sget_simple_extent_ndims_f(dspace,ndims,ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace dimensionality"
       return
    ENDIF
    CALL DebugMessage("     > No. of dim.: ",ndims) 

    CALL h5sget_simple_extent_dims_f(dspace,datadims,maxdatadims,ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace size"
       return
    ENDIF

    dims(1) = datadims(1)
    CALL DebugMessage("     > Size dim. 1: ",int(dims(1)))

    ! Read the dataset and transfer the result

    ALLOCATE(H5dataset(dims(1)),stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error allocating H5dataset"
       return
    ENDIF

    CALL DebugMessage(" --- Reading the data")
    IF ( ltype == H5G_DATASET_F) THEN
       CALL h5dread_f(d_id, H5T_NATIVE_DOUBLE, H5dataset, dims, ErrorFlag)
    ELSE  !  ltype.eqv.3
       CALL h5aread_f(d_id, H5T_NATIVE_DOUBLE, H5dataset, dims, ErrorFlag)
    ENDIF
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error reading data"
       return
    ENDIF

    ALLOCATE(dataset(dims(1)),stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error allocating dataset"
       return
    ENDIF

    dataset = H5dataset

    DEALLOCATE(H5dataset,stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error deallocating H5dataset"
       return
    ENDIF

    ! Close all that is open

    CALL H5Read_close

  END SUBROUTINE H5ReadDouble1D

! *******************************************************************

  SUBROUTINE H5ReadDouble2D      &
       (                         &
       H5filename,               & !Coming in
       datasetname,              & !Coming in
       dataset                   & !Going out
       )
    !
    !<<<<<<<<<<<<<<<<<<<<<<< Implicit statement >>>>>>>>>>>>>>>>>>>>>
    IMPLICIT NONE

    !<<<<<<<<<<<<<<<<<<<<<<< Arguments >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !Coming in
    CHARACTER (len = *), INTENT(in)                      :: H5filename
    CHARACTER (len = *), INTENT(in)                      :: datasetname

    !Going out:
    REAL(KIND=DP),  DIMENSION(:,:), POINTER              :: dataset

    !<<<<<<<<<<<<<<<<<<<<<<< Local variables >>>>>>>>>>>>>>>>>>>>>>>>
    INTEGER, PARAMETER                                   :: maxdims = 4
    INTEGER                                              :: ndims
    INTEGER(hsize_t), DIMENSION(maxdims)                 :: datadims
    INTEGER(hsize_t), DIMENSION(maxdims)                 :: maxdatadims
    INTEGER(hsize_t), DIMENSION(:), ALLOCATABLE          :: dims
    REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE, TARGET   :: H5dataset
    INTEGER                                              :: ltype

    !<<<<<<<<<<<<<<<<<<<<<<< Start of routine code >>>>>>>>>>>>>>>>>>

    CALL DebugMessage(" === in H5ReadDouble2D ===")

    ! Initialise the HDF5 file and open all levels up to
    ! the one that needs to be read.

    CALL H5Read_init ( H5filename, datasetname )
    IF (ErrorFlag.lt.0) return

    ! Get details of the dataset

    CALL DebugMessage(" --- Determining details of the data")

    ltype=OpenLevels_type(NLevels)

    IF ( ltype == H5G_DATASET_F) THEN
       CALL h5dget_space_f(d_id,dspace,ErrorFlag)
    ELSE  !  ltype.eqv.3
       CALL h5aget_space_f(d_id,dspace,ErrorFlag)
    ENDIF
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace"
       return
    ENDIF

    CALL h5sget_simple_extent_ndims_f(dspace,ndims,ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace dimensionality"
       return
    ENDIF
    CALL DebugMessage("     > No. of dim.: ",ndims) 

    ALLOCATE(dims(ndims),stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error allocating dims"
       return
    ENDIF

    CALL h5sget_simple_extent_dims_f(dspace,datadims,maxdatadims,ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace size"
       return
    ENDIF

    dims(1) = datadims(1)
    dims(2) = datadims(2)
    CALL DebugMessage("     > Size dim. 1: ",int(dims(1)))
    CALL DebugMessage("     > Size dim. 2: ",int(dims(2)))

    ! Read the dataset and transfer the result

    ALLOCATE(H5dataset(dims(1),dims(2)),stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error allocating H5dataset"
       return
    ENDIF

    CALL DebugMessage(" --- Reading the data")
    IF ( ltype == H5G_DATASET_F) THEN
       CALL h5dread_f(d_id, H5T_NATIVE_DOUBLE, H5dataset, dims, ErrorFlag)
    ELSE  !  ltype.eqv.3
       CALL h5aread_f(d_id, H5T_NATIVE_DOUBLE, H5dataset, dims, ErrorFlag)
    ENDIF
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error reading data"
       return
    ENDIF

    ALLOCATE(dataset(dims(1),dims(2)),stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error allocating dataset"
       return
    ENDIF

    dataset = H5dataset

    DEALLOCATE(H5dataset,stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error deallocating H5dataset"
       return
    ENDIF

    DEALLOCATE(dims,stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error deallocating dims"
       return
    ENDIF

    ! Close all that is open

    CALL H5Read_close

  END SUBROUTINE H5ReadDouble2D

! *******************************************************************

  SUBROUTINE H5ReadDouble3D      &
       (                         &
       H5filename,               & !Coming in
       datasetname,              & !Coming in
       dataset                   & !Going out
       )
    !
    !<<<<<<<<<<<<<<<<<<<<<<< Implicit statement >>>>>>>>>>>>>>>>>>>>>
    IMPLICIT NONE

    !<<<<<<<<<<<<<<<<<<<<<<< Arguments >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !Coming in
    CHARACTER (len = *), INTENT(in)                      :: H5filename
    CHARACTER (len = *), INTENT(in)                      :: datasetname

    !Going out:
    REAL(KIND=DP), DIMENSION(:,:,:), POINTER             :: dataset

    !<<<<<<<<<<<<<<<<<<<<<<< Local variables >>>>>>>>>>>>>>>>>>>>>>>>
    INTEGER, PARAMETER                                   :: maxdims = 4
    INTEGER                                              :: ndims
    INTEGER(hsize_t), DIMENSION(maxdims)                 :: datadims
    INTEGER(hsize_t), DIMENSION(maxdims)                 :: maxdatadims
    INTEGER(hsize_t), DIMENSION(:), ALLOCATABLE          :: dims
    REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: H5dataset
    INTEGER                                              :: ltype

    !<<<<<<<<<<<<<<<<<<<<<<< Start of routine code >>>>>>>>>>>>>>>>>>

    CALL DebugMessage(" === in H5ReadDouble3D ===")

    ! Initialise the HDF5 file and open all levels up to
    ! the one that needs to be read.

    CALL H5Read_init ( H5filename, datasetname )
    IF (ErrorFlag.lt.0) return

    ! Get details of the dataset

    CALL DebugMessage(" --- Determining details of the data")

    ltype=OpenLevels_type(NLevels)

    IF ( ltype == H5G_DATASET_F) THEN
       CALL h5dget_space_f(d_id,dspace,ErrorFlag)
    ELSE  !  ltype.eqv.3
       CALL h5aget_space_f(d_id,dspace,ErrorFlag)
    ENDIF
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace"
       return
    ENDIF

    CALL h5sget_simple_extent_ndims_f(dspace,ndims,ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace dimensionality"
       return
    ENDIF
    CALL DebugMessage("     > No. of dim.: ",ndims) 

    ALLOCATE(dims(ndims),stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error allocating dims"
       return
    ENDIF

    CALL h5sget_simple_extent_dims_f(dspace,datadims,maxdatadims,ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace size"
       return
    ENDIF

    dims(1) = datadims(1)
    dims(2) = datadims(2)
    dims(3) = datadims(3)
    CALL DebugMessage("     > Size dim. 1: ",int(dims(1)))
    CALL DebugMessage("     > Size dim. 2: ",int(dims(2)))
    CALL DebugMessage("     > Size dim. 3: ",int(dims(3)))

    ! Read the dataset and transfer the result

    ALLOCATE(H5dataset(dims(1),dims(2),dims(3)),stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error allocating H5dataset"
       return
    ENDIF

    CALL DebugMessage(" --- Reading the data")
    IF ( ltype == H5G_DATASET_F) THEN
       CALL h5dread_f(d_id, H5T_NATIVE_DOUBLE, H5dataset, dims, ErrorFlag)
    ELSE  !  ltype.eqv.3
       CALL h5aread_f(d_id, H5T_NATIVE_DOUBLE, H5dataset, dims, ErrorFlag)
    ENDIF
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error reading data"
       return
    ENDIF

    ALLOCATE(dataset(dims(1),dims(2),dims(3)),stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error allocating H5dataset"
       return
    ENDIF

    dataset = H5dataset

    DEALLOCATE(H5dataset,stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error deallocating H5dataset"
       return
    ENDIF

    DEALLOCATE(dims,stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error deallocating dims"
       return
    ENDIF

    ! Close all that is open

    CALL H5Read_close

  END SUBROUTINE H5ReadDouble3D


! *******************************************************************
! *** Reading 0- and 1-D Character Datasets or Attributes         ***
! *******************************************************************
  ! Not made/tested for higher dimensional sets


  SUBROUTINE H5ReadString0D      &
       (                         &
       H5filename,               & !Coming in
       datasetname,              & !Coming in
       dataset                   & !Going out
       )
    !
    !<<<<<<<<<<<<<<<<<<<<<<< Implicit statement >>>>>>>>>>>>>>>>>>>>>
    implicit none

    !<<<<<<<<<<<<<<<<<<<<<<< Arguments >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !Coming in
    CHARACTER (len = *), INTENT(in)                      :: H5filename
    CHARACTER (len = *), INTENT(in)                      :: datasetname

    !Going out:
    CHARACTER(LEN=*)                                     :: dataset

    !<<<<<<<<<<<<<<<<<<<<<<< Local variables >>>>>>>>>>>>>>>>>>>>>>>>

    
    
    INTEGER(hsize_t), DIMENSION(1)                       :: dims
    CHARACTER(len=LEN(dataset))                          :: H5dataset
    INTEGER(SIZE_T)                                      :: strlen
    INTEGER                                              :: ltype

    !<<<<<<<<<<<<<<<<<<<<<<< Start of routine code >>>>>>>>>>>>>>>>>>

    CALL DebugMessage(" === in H5ReadString0D ===")

    ! Initialise the HDF5 file and open all levels up to
    ! the one that needs to be read.

    CALL H5Read_init( H5filename, datasetname)
    IF (ErrorFlag.lt.0) return

    ! Get details of the dataset

    CALL DebugMessage(" --- Determining details of the data")

    ltype=OpenLevels_type(NLevels)

    IF ( ltype == H5G_DATASET_F) THEN
       CALL h5dget_type_f(d_id,type_id,ErrorFlag)
    ELSE  !  ltype.eqv.3
       CALL h5aget_type_f(d_id,type_id,ErrorFlag)
    ENDIF
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining datatype"
       return
    ENDIF

    ! Strangely enough type_id (201326890) is not equal to either
    ! H5T_NATIVE_CHARACTER (201326866) or H5T_STRING (201326889)

    ! The variable 'strlen' determined with h5tget_size_f gives
    ! the defined length of the string, not the length of the
    ! text in that string: that text may be followed by spaces
    ! and/or a NULL character. To make sure that in particular
    ! the latter does not appear in the final string, use the
    ! real length of the text in the string.

    CALL h5tget_size_f(type_id,strlen,ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining string length"
       return
    ENDIF

    IF (strlen.gt.LEN(dataset)) THEN
       ErrorMessage=" *** Error: string of insufficient length defined"
       ErrorFlag=-1
       return
    ENDIF

    dims=1

    CALL DebugMessage("     > String length  :",int(strlen))
    CALL DebugMessage("     > Defined length :",LEN(H5dataset))

    H5dataset = ""

    ! Read the dataset and transfer the result

    CALL DebugMessage(" --- Reading the data")
    IF ( ltype == H5G_DATASET_F) THEN
       CALL h5dread_f(d_id, type_id, H5dataset, dims, ErrorFlag)
    ELSE  !  ltype.eqv.3
       CALL h5aread_f(d_id, type_id, H5dataset, dims, ErrorFlag)
    ENDIF
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error reading data"
       return
    ENDIF

    dataset=H5dataset(1:LengString(H5dataset))

    ! Close all that is open

    CALL h5tclose_f(type_id, ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error closing data type"
       return
    ENDIF

    CALL H5Read_close

  END SUBROUTINE H5ReadString0D

! *******************************************************************

  SUBROUTINE H5ReadString1D      &
       (                         &
       H5filename,               & !Coming in
       datasetname,              & !Coming in
       dataset                   & !Going out
       )
    !
    !<<<<<<<<<<<<<<<<<<<<<<< Implicit statement >>>>>>>>>>>>>>>>>>>>>
    implicit none

    !<<<<<<<<<<<<<<<<<<<<<<< Arguments >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !Coming in
    CHARACTER (len = *), INTENT(in)                      :: H5filename
    CHARACTER (len = *), INTENT(in)                      :: datasetname

    !Going out:
    CHARACTER(LEN=*), DIMENSION(:), POINTER              :: dataset

    !<<<<<<<<<<<<<<<<<<<<<<< Local variables >>>>>>>>>>>>>>>>>>>>>>>>
    INTEGER                                              :: ndims
    INTEGER(hsize_t), DIMENSION(1)                       :: datadims
    INTEGER(hsize_t), DIMENSION(1)                       :: maxdatadims
    INTEGER(hsize_t), DIMENSION(1)                       :: dims
    CHARACTER(len=LEN(dataset)), DIMENSION(:), ALLOCATABLE, TARGET :: H5dataset
    INTEGER(SIZE_T)                                      :: strlen
    INTEGER                                              :: ltype

    CHARACTER(len=LEN(dataset))                          :: substr
    INTEGER                                              :: indx, idts
    INTEGER                                              :: ibeg, iend
    INTEGER                                              :: deflen

    !<<<<<<<<<<<<<<<<<<<<<<< Start of routine code >>>>>>>>>>>>>>>>>>

    CALL DebugMessage(" === in H5ReadString1D ===")

    ! Initialise the HDF5 file and open all levels up to
    ! the one that needs to be read.

    CALL H5Read_init( H5filename, datasetname)
    IF (ErrorFlag.lt.0) return

    ! Get details of the dataset

    CALL DebugMessage(" --- Determining details of the data")

    ltype=OpenLevels_type(NLevels)

    IF ( ltype == H5G_DATASET_F) THEN
       CALL h5dget_type_f(d_id,type_id,ErrorFlag)
    ELSE  !  ltype.eqv.3
       CALL h5aget_type_f(d_id,type_id,ErrorFlag)
    ENDIF
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining datatype"
       return
    ENDIF

    ! Strangely enough type_id (201326891) is not equal to either
    ! H5T_NATIVE_CHARACTER (201326866) or H5T_STRING (201326889),
    ! nor is it equal to the type_id of a 0-dim string (201326890).

    ! The variable 'strlen' determined with h5tget_size_f gives
    ! the defined length of the string, not the length of the
    ! text in that string: that text may be followed by spaces
    ! and/or a NULL character.
    !
    ! Example: strlen=8, with as text "Nominal", which is
    ! 7 characters long; the 8-th character is then a NULL.
    !
    ! To make sure that in particular the NULL does not appear
    ! in the final string (since Fortran does not like that),
    ! use the real length of the text in the final string.

    CALL h5tget_size_f(type_id,strlen,ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining string length"
       return
    ENDIF

    IF (strlen.gt.LEN(dataset)) THEN
       ErrorMessage=" *** Error: string of insufficient length defined"
       ErrorFlag=-1
       return
    ENDIF

    IF ( ltype == H5G_DATASET_F) THEN
       CALL h5dget_space_f(d_id,dspace,ErrorFlag)
    ELSE  !  ltype.eqv.3
       CALL h5aget_space_f(d_id,dspace,ErrorFlag)
    ENDIF
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace"
       return
    ENDIF

    CALL h5sget_simple_extent_ndims_f(dspace,ndims,ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace dimensionality"
       return
    ENDIF

    CALL h5sget_simple_extent_dims_f(dspace,datadims,maxdatadims,ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace size"
       return
    ENDIF

    dims(1) = datadims(1)

    ! Read the dataset and transfer the result

    ALLOCATE(H5dataset(dims(1)),stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error allocating dims"
       return
    ENDIF

    CALL DebugMessage("     > String length  :",INT(strlen))
    CALL DebugMessage("     > Defined length :",LEN(H5dataset))
    CALL DebugMessage("     > Array dimension:",SIZE(H5dataset))

    H5dataset = ""

    CALL DebugMessage(" --- Reading the data")
    IF ( ltype == H5G_DATASET_F) THEN
       CALL h5dread_f(d_id, type_id, H5dataset, dims, ErrorFlag)
    ELSE  !  ltype.eqv.3
       CALL h5aread_f(d_id, type_id, H5dataset, dims, ErrorFlag)
    ENDIF
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error reading data"
       return
    ENDIF

    ! The texts are now stored in H5dataset. This is _not_
    ! in the form of one strlen long string per array (as 
    ! one would expect), but an X number of strlen strings
    ! pasted together to exactly fill the defined length of
    ! H5dataset, where X is a real, not an integer.
    ! Example (as above): strlen=8, with "Nominal" as text,
    ! followed by a NULL to fill the rest (see above). With
    ! a 30-character defined string:
    !    print*,H5dataset(1)
    ! gives (indicating NULL by @):
    !    Nominal@Nominal@Nominal@Nomina
    !   0....5....0....5....0....5....0
    ! The rest of the last text is then in H5dataset(2) ...
    ! This makes transferring the string to the final array
    ! in need of some bookkeeping.

    ALLOCATE(dataset(dims(1)),stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error allocating dataset"
       return
    ENDIF

    CALL DebugMessage(" --- Reconstructing the individual strings")
    deflen=LEN(H5dataset)        ! defined length of H5dataset
    idts=1                       ! index in the H5dataset
    ibeg=1                       ! start index along string at index idts

    do indx=1,dims(1)
       iend=ibeg+strlen-1        ! end index along string at index idts
       if (iend.le.deflen) then  ! still fully in current idts
          substr = H5dataset(idts)(ibeg:iend)
          ibeg=iend+1
       else                      ! partly in current idts, partly in next idts
          substr = H5dataset(idts)(ibeg:deflen) &
                // H5dataset(idts+1)(1:strlen-(deflen-ibeg+1))
          idts=idts+1
          ibeg=strlen-(deflen-ibeg+1)+1
       endif
       dataset(indx) = substr(1:LengString(substr))
    enddo

    DEALLOCATE(H5dataset,stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error deallocating H5dataset"
       return
    ENDIF

    ! Close all that is open

    CALL h5tclose_f(type_id, ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error closing data type"
       return
    ENDIF

    CALL H5Read_close

  END SUBROUTINE H5ReadString1D


! *******************************************************************
! *** Reading 0- and 1-D Integer Compound Datasets or Attributes  ***
! *******************************************************************
  ! Not made/tested for higher dimensional sets


  SUBROUTINE H5ReadCompInteger0D &
       (                         &
       H5filename,               & !Coming in
       datasetname,              & !Coming in
       compoundsetname,          & !Coming in
       dataset                   & !Going out
       )
    !
    !<<<<<<<<<<<<<<<<<<<<<<<<< Implicit statement >>>>>>>>>>>>>>>>>>>>>>>>>
    IMPLICIT NONE

    !<<<<<<<<<<<<<<<<<<<<<<<<< Arguments >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !Coming in
    CHARACTER (len = *), INTENT(in)                      :: H5filename
    CHARACTER (len = *), INTENT(in)                      :: datasetname
    CHARACTER (len = *), INTENT(in)                      :: compoundsetname

    !Going out:
    INTEGER                                              :: dataset

    !<<<<<<<<<<<<<<<<<<<<<<<<< Local variables >>>>>>>>>>>>>>>>>>>>>>>>>>>>
    INTEGER, PARAMETER                                   :: maxdims = 4
 
   
   
    INTEGER(hsize_t), DIMENSION(1)                       :: dims
    INTEGER                                              :: H5dataset
    INTEGER                                              :: ltype

    INTEGER(SIZE_T)                                      :: type_sizeI
    INTEGER (SIZE_T)                                              :: offset

    !<<<<<<<<<<<<<<<<<<<<<<<< Start of routine code >>>>>>>>>>>>>>>>>>>>>>>

    CALL DebugMessage(" === in H5ReadCompInteger0D ===")

    ! Initialise the HDF5 file and open all levels up to
    ! the one that needs to be read.

    CALL H5Read_init ( H5filename, datasetname )
    IF (ErrorFlag.lt.0) return

    ! Get details of the dataset

    CALL DebugMessage(" --- Determining details of the data")

    ltype=OpenLevels_type(NLevels)

    dims=1

    ! Create a field in the memory for the dataset and insert it

    CALL h5tget_size_f(H5T_NATIVE_INTEGER, type_sizeI, ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining data type"
       return
    ENDIF

    CALL h5tcreate_f(H5T_COMPOUND_F, type_sizeI, type_id, ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining creating field"
       return
    ENDIF

    offset = 0
    CALL h5tinsert_f(type_id, compoundsetname, offset, &
            H5T_NATIVE_INTEGER, ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining inserting field"
       return
    ENDIF

    ! Read the dataset and transfer the result

    CALL DebugMessage(" --- Reading the data")
    IF ( ltype == H5G_DATASET_F) THEN
       CALL h5dread_f(d_id, type_id, H5dataset, dims, ErrorFlag)
    ELSE  !  ltype.eqv.3
       CALL h5aread_f(d_id, type_id, H5dataset, dims, ErrorFlag)
    ENDIF
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error reading data"
       return
    ENDIF

    dataset = H5dataset

    ! Close all that is open

    CALL h5tclose_f(type_id,ErrorFlag)
    CALL H5Read_close

  END SUBROUTINE H5ReadCompInteger0D

! *******************************************************************

  SUBROUTINE H5ReadCompInteger1D &
       (                         &
       H5filename,               & !Coming in
       datasetname,              & !Coming in
       compoundsetname,          & !Coming in
       dataset                   & !Going out
       )
    !
    !<<<<<<<<<<<<<<<<<<<<<<<<< Implicit statement >>>>>>>>>>>>>>>>>>>>>>>>>
    IMPLICIT NONE

    !<<<<<<<<<<<<<<<<<<<<<<<<< Arguments >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !Coming in
    CHARACTER (len = *), INTENT(in)                      :: H5filename
    CHARACTER (len = *), INTENT(in)                      :: datasetname
    CHARACTER (len = *), INTENT(in)                      :: compoundsetname

    !Going out:
    INTEGER, DIMENSION(:), POINTER                       :: dataset

    !<<<<<<<<<<<<<<<<<<<<<<<<< Local variables >>>>>>>>>>>>>>>>>>>>>>>>>>>>
    INTEGER, PARAMETER                                   :: maxdims = 4
    INTEGER                                              :: ndims
    
    INTEGER(hsize_t), DIMENSION(1)                       :: datadims
    INTEGER(hsize_t), DIMENSION(1)                       :: maxdatadims
    INTEGER(hsize_t), DIMENSION(1)                       :: dims
    INTEGER, DIMENSION(:), ALLOCATABLE, TARGET           :: H5dataset
    INTEGER                                              :: ltype

    INTEGER(SIZE_T)                                      :: type_sizeI
    INTEGER (SIZE_T)                                             :: offset

    !<<<<<<<<<<<<<<<<<<<<<<<< Start of routine code >>>>>>>>>>>>>>>>>>>>>>>

    CALL DebugMessage(" === in H5ReadCompInteger1D ===")

    ! Initialise the HDF5 file and open all levels up to
    ! the one that needs to be read.

    CALL H5Read_init ( H5filename, datasetname )
    IF (ErrorFlag.lt.0) return

    ! Get details of the dataset

    CALL DebugMessage(" --- Determining details of the data")

    ltype=OpenLevels_type(NLevels)

    IF ( ltype == H5G_DATASET_F) THEN
       CALL h5dget_space_f(d_id,dspace,ErrorFlag)
    ELSE  !  ltype.eqv.3
       CALL h5aget_space_f(d_id,dspace,ErrorFlag)
    ENDIF
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace"
       return
    ENDIF

    CALL h5sget_simple_extent_ndims_f(dspace,ndims,ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace dimensionality"
       return
    ENDIF
    CALL DebugMessage("     > No. of dim.: ",ndims) 

    CALL h5sget_simple_extent_dims_f(dspace,datadims,maxdatadims,ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace size"
       return
    ENDIF

    dims(1) = datadims(1)
    CALL DebugMessage("     > Size dim. 1: ",int(dims(1)))

    ! Create a field in the memory for the dataset and insert it

    CALL h5tget_size_f(H5T_NATIVE_INTEGER, type_sizeI, ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining data type"
       return
    ENDIF

    CALL h5tcreate_f(H5T_COMPOUND_F, type_sizeI, type_id, ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining creating field"
       return
    ENDIF

    offset = 0
    CALL h5tinsert_f(type_id, compoundsetname, offset, &
            H5T_NATIVE_INTEGER, ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining inserting field"
       return
    ENDIF

    ! Read the dataset and transfer the result

    ALLOCATE(H5dataset(dims(1)),stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error allocating H5dataset"
       return
    ENDIF

    CALL DebugMessage(" --- Reading the data")
    IF ( ltype == H5G_DATASET_F) THEN
       CALL h5dread_f(d_id, type_id, H5dataset, dims, ErrorFlag)
    ELSE  !  ltype.eqv.3
       CALL h5aread_f(d_id, type_id, H5dataset, dims, ErrorFlag)
    ENDIF
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error reading data"
       return
    ENDIF

    ALLOCATE(dataset(dims(1)),stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error allocating dataset"
       return
    ENDIF

    dataset = H5dataset

    DEALLOCATE(H5dataset,stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error deallocating H5dataset"
       return
    ENDIF

    ! Close all that is open

    CALL h5tclose_f(type_id,ErrorFlag)
    CALL H5Read_close

  END SUBROUTINE H5ReadCompInteger1D


! *******************************************************************
! *** Reading 0- and 1-D Real Compound Datasets or Attributes     ***
! *******************************************************************
  ! Not made/tested for higher dimensional sets


  SUBROUTINE H5ReadCompReal0D    &
       (                         &
       H5filename,               & !Coming in
       datasetname,              & !Coming in
       compoundsetname,          & !Coming in
       dataset                   & !Going out
       )
    !
    !<<<<<<<<<<<<<<<<<<<<<<<<< Implicit statement >>>>>>>>>>>>>>>>>>>>>>>>>
    IMPLICIT NONE

    !<<<<<<<<<<<<<<<<<<<<<<<<< Arguments >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !Coming in
    CHARACTER (len = *), INTENT(in)                      :: H5filename
    CHARACTER (len = *), INTENT(in)                      :: datasetname
    CHARACTER (len = *), INTENT(in)                      :: compoundsetname

    !Going out:
    REAL(KIND=SP)                                        :: dataset

    !<<<<<<<<<<<<<<<<<<<<<<<<< Local variables >>>>>>>>>>>>>>>>>>>>>>>>>>>>
    INTEGER, PARAMETER                                   :: maxdims = 4
  
   
    INTEGER(hsize_t), DIMENSION(1)                       :: dims
    REAL(KIND=SP)                                        :: H5dataset
    INTEGER                                              :: ltype

    INTEGER(SIZE_T)                                      :: type_sizeR
    INTEGER (SIZE_T)                                              :: offset

    !<<<<<<<<<<<<<<<<<<<<<<<< Start of routine code >>>>>>>>>>>>>>>>>>>>>>>

    CALL DebugMessage(" === in H5ReadCompReal0D ===")

    ! Initialise the HDF5 file and open all levels up to
    ! the one that needs to be read.

    CALL H5Read_init ( H5filename, datasetname )
    IF (ErrorFlag.lt.0) return

    ! Get details of the dataset

    CALL DebugMessage(" --- Determining details of the data")

    ltype=OpenLevels_type(NLevels)

    dims=1

    ! Create a field in the memory for the dataset and insert it

    CALL h5tget_size_f(H5T_NATIVE_REAL, type_sizeR, ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining data type"
       return
    ENDIF

    CALL h5tcreate_f(H5T_COMPOUND_F, type_sizeR, type_id, ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining creating field"
       return
    ENDIF

    offset = 0
    CALL h5tinsert_f(type_id, compoundsetname, offset, &
            H5T_NATIVE_REAL, ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining inserting field"
       return
    ENDIF

    ! Read the dataset and transfer the result

    CALL DebugMessage(" --- Reading the data")
    IF ( ltype == H5G_DATASET_F) THEN
       CALL h5dread_f(d_id, type_id, H5dataset, dims, ErrorFlag)
    ELSE  !  ltype.eqv.3
       CALL h5aread_f(d_id, type_id, H5dataset, dims, ErrorFlag)
    ENDIF
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error reading data"
       return
    ENDIF

    dataset = H5dataset

    ! Close all that is open

    CALL h5tclose_f(type_id,ErrorFlag)
    CALL H5Read_close

  END SUBROUTINE H5ReadCompReal0D

! *******************************************************************

  SUBROUTINE H5ReadCompReal1D    &
       (                         &
       H5filename,               & !Coming in
       datasetname,              & !Coming in
       compoundsetname,          & !Coming in
       dataset                   & !Going out
       )
    !
    !<<<<<<<<<<<<<<<<<<<<<<<<< Implicit statement >>>>>>>>>>>>>>>>>>>>>>>>>
    IMPLICIT NONE

    !<<<<<<<<<<<<<<<<<<<<<<<<< Arguments >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !Coming in
    CHARACTER (len = *), INTENT(in)                      :: H5filename
    CHARACTER (len = *), INTENT(in)                      :: datasetname
    CHARACTER (len = *), INTENT(in)                      :: compoundsetname

    !Going out:
    REAL(KIND=SP), DIMENSION(:), POINTER                 :: dataset

    !<<<<<<<<<<<<<<<<<<<<<<<<< Local variables >>>>>>>>>>>>>>>>>>>>>>>>>>>>
    INTEGER, PARAMETER                                   :: maxdims = 4
    INTEGER                                              :: ndims
   
    INTEGER(hsize_t), DIMENSION(1)                       :: datadims
    INTEGER(hsize_t), DIMENSION(1)                       :: maxdatadims
    INTEGER(hsize_t), DIMENSION(1)                       :: dims
    REAL(KIND=SP), DIMENSION(:), ALLOCATABLE, TARGET     :: H5dataset
    INTEGER                                              :: ltype

    INTEGER(SIZE_T)                                      :: type_sizeR
    INTEGER(SIZE_T)                                               :: offset

    !<<<<<<<<<<<<<<<<<<<<<<<< Start of routine code >>>>>>>>>>>>>>>>>>>>>>>

    CALL DebugMessage(" === in H5ReadCompReal1D ===")

    ! Initialise the HDF5 file and open all levels up to
    ! the one that needs to be read.

    CALL H5Read_init ( H5filename, datasetname )
    IF (ErrorFlag.lt.0) return

    ! Get details of the dataset

    CALL DebugMessage(" --- Determining details of the data")

    ltype=OpenLevels_type(NLevels)

    IF ( ltype == H5G_DATASET_F) THEN
       CALL h5dget_space_f(d_id,dspace,ErrorFlag)
    ELSE  !  ltype.eqv.3
       CALL h5aget_space_f(d_id,dspace,ErrorFlag)
    ENDIF
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace"
       return
    ENDIF

    CALL h5sget_simple_extent_ndims_f(dspace,ndims,ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace dimensionality"
       return
    ENDIF
    CALL DebugMessage("     > No. of dim.: ",ndims) 

    CALL h5sget_simple_extent_dims_f(dspace,datadims,maxdatadims,ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace size"
       return
    ENDIF

    dims(1) = datadims(1)
    CALL DebugMessage("     > Size dim. 1: ",int(dims(1)))

    ! Create a field in the memory for the dataset and insert it

    CALL h5tget_size_f(H5T_NATIVE_REAL, type_sizeR, ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining data type"
       return
    ENDIF

    CALL h5tcreate_f(H5T_COMPOUND_F, type_sizeR, type_id, ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining creating field"
       return
    ENDIF

    offset = 0
    CALL h5tinsert_f(type_id, compoundsetname, offset, &
            H5T_NATIVE_REAL, ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining inserting field"
       return
    ENDIF

    ! Read the dataset and transfer the result

    ALLOCATE(H5dataset(dims(1)),stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error allocating H5dataset"
       return
    ENDIF

    CALL DebugMessage(" --- Reading the data")
    IF ( ltype == H5G_DATASET_F) THEN
       CALL h5dread_f(d_id, type_id, H5dataset, dims, ErrorFlag)
    ELSE  !  ltype.eqv.3
       CALL h5aread_f(d_id, type_id, H5dataset, dims, ErrorFlag)
    ENDIF
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error reading data"
       return
    ENDIF

    ALLOCATE(dataset(dims(1)),stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error allocating dataset"
       return
    ENDIF

    dataset = H5dataset

    DEALLOCATE(H5dataset,stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error deallocating H5dataset"
       return
    ENDIF

    ! Close all that is open

    CALL h5tclose_f(type_id,ErrorFlag)
    CALL H5Read_close

  END SUBROUTINE H5ReadCompReal1D


! *******************************************************************
! *** Reading 0- and 1-D Double Compound Datasets or Attributes   ***
! *******************************************************************
  ! Not made/tested for higher dimensional sets


  SUBROUTINE H5ReadCompDouble0D  &
       (                         &
       H5filename,               & !Coming in
       datasetname,              & !Coming in
       compoundsetname,          & !Coming in
       dataset                   & !Going out
       )
    !
    !<<<<<<<<<<<<<<<<<<<<<<<<< Implicit statement >>>>>>>>>>>>>>>>>>>>>>>>>
    IMPLICIT NONE

    !<<<<<<<<<<<<<<<<<<<<<<<<< Arguments >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !Coming in
    CHARACTER (len = *), INTENT(in)                      :: H5filename
    CHARACTER (len = *), INTENT(in)                      :: datasetname
    CHARACTER (len = *), INTENT(in)                      :: compoundsetname

    !Going out:
    REAL(KIND=DP)                                        :: dataset

    !<<<<<<<<<<<<<<<<<<<<<<<<< Local variables >>>>>>>>>>>>>>>>>>>>>>>>>>>>
    INTEGER, PARAMETER                                   :: maxdims = 4
  
    
    
    
    INTEGER(hsize_t), DIMENSION(1)                       :: dims
    REAL(KIND=DP)                                        :: H5dataset
    INTEGER                                              :: ltype

    INTEGER(SIZE_T)                                      :: type_sizeD
    INTEGER(SIZE_T)                                               :: offset

    !<<<<<<<<<<<<<<<<<<<<<<<< Start of routine code >>>>>>>>>>>>>>>>>>>>>>>

    CALL DebugMessage(" === in H5ReadCompDouble0D ===")

    ! Initialise the HDF5 file and open all levels up to
    ! the one that needs to be read.

    CALL H5Read_init ( H5filename, datasetname )
    IF (ErrorFlag.lt.0) return

    ! Get details of the dataset

    CALL DebugMessage(" --- Determining details of the data")

    ltype=OpenLevels_type(NLevels)

    dims=1

    ! Create a field in the memory for the dataset and insert it

    CALL h5tget_size_f(H5T_NATIVE_DOUBLE, type_sizeD, ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining data type"
       return
    ENDIF

    CALL h5tcreate_f(H5T_COMPOUND_F, type_sizeD, type_id, ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining creating field"
       return
    ENDIF

    offset = 0
    CALL h5tinsert_f(type_id, compoundsetname, offset, &
            H5T_NATIVE_DOUBLE, ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining inserting field"
       return
    ENDIF

    ! Read the dataset and transfer the result

    CALL DebugMessage(" --- Reading the data")
    IF ( ltype == H5G_DATASET_F) THEN
       CALL h5dread_f(d_id, type_id, H5dataset, dims, ErrorFlag)
    ELSE  !  ltype.eqv.3
       CALL h5aread_f(d_id, type_id, H5dataset, dims, ErrorFlag)
    ENDIF
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error reading data"
       return
    ENDIF

    dataset = H5dataset

    ! Close all that is open

    CALL h5tclose_f(type_id,ErrorFlag)
    CALL H5Read_close

  END SUBROUTINE H5ReadCompDouble0D

! *******************************************************************

  SUBROUTINE H5ReadCompDouble1D  &
       (                         &
       H5filename,               & !Coming in
       datasetname,              & !Coming in
       compoundsetname,          & !Coming in
       dataset                   & !Going out
       )
    !
    !<<<<<<<<<<<<<<<<<<<<<<<<< Implicit statement >>>>>>>>>>>>>>>>>>>>>>>>>
    IMPLICIT NONE

    !<<<<<<<<<<<<<<<<<<<<<<<<< Arguments >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !Coming in
    CHARACTER (len = *), INTENT(in)                      :: H5filename
    CHARACTER (len = *), INTENT(in)                      :: datasetname
    CHARACTER (len = *), INTENT(in)                      :: compoundsetname

    !Going out:
    REAL(KIND=DP), DIMENSION(:), POINTER                 :: dataset

    !<<<<<<<<<<<<<<<<<<<<<<<<< Local variables >>>>>>>>>>>>>>>>>>>>>>>>>>>>
    INTEGER, PARAMETER                                   :: maxdims = 4
    INTEGER                                              :: ndims
   
    INTEGER(hsize_t), DIMENSION(1)                       :: datadims
    INTEGER(hsize_t), DIMENSION(1)                       :: maxdatadims
    INTEGER(hsize_t), DIMENSION(1)                       :: dims
    REAL(KIND=DP), DIMENSION(:), ALLOCATABLE, TARGET     :: H5dataset

    INTEGER(SIZE_T)                                      :: type_sizeD
    INTEGER (SIZE_T)                                      :: offset
    INTEGER                                              :: ltype

    !<<<<<<<<<<<<<<<<<<<<<<<< Start of routine code >>>>>>>>>>>>>>>>>>>>>>>

    CALL DebugMessage(" === in H5ReadCompDouble1D ===")

    ! Initialise the HDF5 file and open all levels up to
    ! the one that needs to be read.

    CALL H5Read_init ( H5filename, datasetname )
    IF (ErrorFlag.lt.0) return

    ! Get details of the dataset

    CALL DebugMessage(" --- Determining details of the data")

    ltype=OpenLevels_type(NLevels)

    IF ( ltype == H5G_DATASET_F) THEN
       CALL h5dget_space_f(d_id,dspace,ErrorFlag)
    ELSE  !  ltype.eqv.3
       CALL h5aget_space_f(d_id,dspace,ErrorFlag)
    ENDIF
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace"
       return
    ENDIF

    CALL h5sget_simple_extent_ndims_f(dspace,ndims,ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace dimensionality"
       return
    ENDIF
    CALL DebugMessage("     > No. of dim.: ",ndims) 

    CALL h5sget_simple_extent_dims_f(dspace,datadims,maxdatadims,ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining dataspace size"
       return
    ENDIF

    dims(1) = datadims(1)
    CALL DebugMessage("     > Size dim. 1: ",int(dims(1)))

    ! Create a field in the memory for the dataset and insert it

    CALL h5tget_size_f(H5T_NATIVE_DOUBLE, type_sizeD, ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining data type"
       return
    ENDIF

    CALL h5tcreate_f(H5T_COMPOUND_F, type_sizeD, type_id, ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining creating field"
       return
    ENDIF

    offset = 0
    CALL h5tinsert_f(type_id, compoundsetname, offset, &
            H5T_NATIVE_DOUBLE, ErrorFlag)
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error determining inserting field"
       return
    ENDIF

    ! Read the dataset and transfer the result

    ALLOCATE(H5dataset(dims(1)),stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error allocating H5dataset"
       return
    ENDIF

    CALL DebugMessage(" --- Reading the data")
    IF ( ltype == H5G_DATASET_F) THEN
       CALL h5dread_f(d_id, type_id, H5dataset, dims, ErrorFlag)
    ELSE  !  ltype.eqv.3
       CALL h5aread_f(d_id, type_id, H5dataset, dims, ErrorFlag)
    ENDIF
    IF (ErrorFlag.lt.0) THEN
       ErrorMessage=" *** Error reading data"
       return
    ENDIF

    ALLOCATE(dataset(dims(1)),stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error allocating dataset"
       return
    ENDIF

    dataset = H5dataset

    DEALLOCATE(H5dataset,stat=AllocStat)
    IF ( AllocStat.ne.0 ) THEN
       ErrorFlag=-1
       ErrorMessage=" *** Error deallocating H5dataset"
       return
    ENDIF

    ! Close all that is open

    CALL h5tclose_f(type_id,ErrorFlag)
    CALL H5Read_close

  END SUBROUTINE H5ReadCompDouble1D


! *******************************************************************
! *******************************************************************

END MODULE ReadH5Dataset

! *******************************************************************
