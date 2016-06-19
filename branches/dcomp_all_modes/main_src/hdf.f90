! $Id$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: hdf.f90 (src)
!       HDF (program)
!
! PURPOSE:	Fortran header file for HDF routines
!
! DESCRIPTION:
!       This file is a modularized version of the 'hdf.f90' file from the 'include'
!       directory of the standard HDF library installation. This version of the file
!       is taken from the HDF4.2r0 library distribution.
!
! AUTHORS:
!  Andrew Heidinger, Andrew.Heidinger@noaa.gov
!  Andi Walther, CIMSS, andi.walther@ssec.wisc.edu
!  Denis Botambekov, CIMSS, denis.botambekov@ssec.wisc.edu
!  William Straka, CIMSS, wstraka@ssec.wisc.edu
!  Last modified by: Aleksandar Jelenak <Aleksandar.Jelenak@noaa.gov>
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
!   LAST MODIFIED: "Thu, 08 Apr 2004 09:50:18 -0400 (ajelenak)"
!
! * Contents: 
! *     Tag definitions
! *     Error return codes
! *    	Logical constants
! * Remarks: This file can be included with Fortran user programs.  As a
! *          general rule, don't use DFNT constants that don't include a
! *          number in their name.  E.g., don't use DFNT_FLOAT, use
! *          DFNT_FLOAT32 or DFNT_FLOAT64.  The DFNT constants that don't
! *          include numbers are for backward compatibility only.  Also,
! *          there are no current plans to support 128-bit number types.
! *          For more information about constants in this file, see the
! *          equivalent constant declarations in the C include file 'hdf.h'
!--------------------------------------------------------------------------------------
module HDF
   implicit none
   private
!	Error Return Codes 
   integer, parameter, public:: DFE_NOERROR = 0
   integer, parameter, public:: DFE_NONE = 0
   integer, parameter, public:: DFE_FNF = -1
   integer, parameter, public:: DFE_DENIED = -2
   integer, parameter, public:: DFE_ALROPEN = -3
   integer, parameter, public:: DFE_TOOMANY = -4
   integer, parameter, public:: DFE_BADNAME = -5
   integer, parameter, public:: DFE_BADACC = -6
   integer, parameter, public:: DFE_BADOPEN = -7
   integer, parameter, public:: DFE_NOTOPEN = -8
   integer, parameter, public:: DFE_CANTCLOSE = -9
   integer, parameter, public:: DFE_DFNULL = -10
   integer, parameter, public:: DFE_ILLTYPE = -11
   integer, parameter, public:: DFE_UNSUPPORTED = -12
   integer, parameter, public:: DFE_BADDDLIST = -13
   integer, parameter, public:: DFE_NOTDFFILE = -14
   integer, parameter, public:: DFE_SEEDTWICE = -15
   integer, parameter, public:: DFE_NOSPACE = -16
   integer, parameter, public:: DFE_NOSUCHTAG = -17
   integer, parameter, public:: DFE_READERROR = -18
 
   integer, parameter, public:: DFE_WRITEERROR = -19
   integer, parameter, public:: DFE_SEEKERROR = -20
   integer, parameter, public:: DFE_NOFREEDD = -21
   integer, parameter, public:: DFE_BADTAG = -22
   integer, parameter, public:: DFE_BADREF = -23
   integer, parameter, public:: DFE_RDONLY = -24
   integer, parameter, public:: DFE_BADCALL = -25
   integer, parameter, public:: DFE_BADPTR = -26
   integer, parameter, public:: DFE_BADLEN = -27
   integer, parameter, public:: DFE_BADSEEK = -28
   integer, parameter, public:: DFE_NOMATCH = -29
   integer, parameter, public:: DFE_NOTINSET = -30
   integer, parameter, public:: DFE_BADDIM = -31
   integer, parameter, public:: DFE_BADOFFSET = -32
   integer, parameter, public:: DFE_BADSCHEME = -33
   integer, parameter, public:: DFE_NODIM = -34
   integer, parameter, public:: DFE_NOTENOUGH = -35
   integer, parameter, public:: DFE_NOVALS = -36
   integer, parameter, public:: DFE_CORRUPT = -37
   integer, parameter, public:: DFE_BADFP = -38

   integer, parameter, public:: DFE_NOREF = -39
   integer, parameter, public:: DFE_BADDATATYPE = -40
   integer, parameter, public:: DFE_BADMCTYPE = -41
   integer, parameter, public:: DFE_BADNUMTYPE = -42
   integer, parameter, public:: DFE_BADORDER = -43
   integer, parameter, public:: DFE_ARGS = -44
   integer, parameter, public:: DFE_INTERNAL = -45
   integer, parameter, public:: DFE_DUPDD = -46
   integer, parameter, public:: DFE_CANTMOD = -47
   integer, parameter, public:: DFE_RANGE = -48
   integer, parameter, public:: DFE_BADTABLE = -49
   integer, parameter, public:: DFE_BADSDG = -50
   integer, parameter, public:: DFE_BADNDG = -51
   integer, parameter, public:: DFE_BADFIELDS = -52
   integer, parameter, public:: DFE_NORESET = -53
   integer, parameter, public:: DFE_NOVS = -54
   integer, parameter, public:: DFE_VGSIZE = -55
   integer, parameter, public:: DFE_DIFFFILES = -56
   integer, parameter, public:: DFE_VTAB = -57
   integer, parameter, public:: DFE_BADAID = -58

   integer, parameter, public:: DFE_OPENAID = -59
   integer, parameter, public:: DFE_BADCONV = -60
   integer, parameter, public:: DFE_GENAPP = -61
   integer, parameter, public:: DFE_CANTFLUSH = -62
   integer, parameter, public:: DFE_BADTYPE = -63
   integer, parameter, public:: DFE_SYMSIZE = -64
   integer, parameter, public:: DFE_BADATTACH = -65
   integer, parameter, public:: DFE_CANTDETACH = -66

! internal file access codes
   integer, parameter, public:: DFACC_READ = 1
   integer, parameter, public:: DFACC_WRITE = 2
   integer, parameter, public:: DFACC_CREATE = 4
   integer, parameter, public:: DFACC_ALL = 7
   integer, parameter, public:: DFACC_RDONLY = 1
   integer, parameter, public:: DFACC_RDWR = 3
   integer, parameter, public:: DFACC_CLOBBER = 4

!	Access types for SDsetaccecellype
   integer, parameter, public:: DFACC_DEFAULT = 0
   integer, parameter, public:: DFACC_SERIAL = 1
   integer, parameter, public:: DFACC_PARALLEL = 9

!	Constants for DFSDsetorder
   integer, parameter, public:: DFO_FORTRAN = 1
   integer, parameter, public:: DFO_C = 2

!	Definitions of storage convention
   integer, parameter, public:: DFNTF_IEEE = 1
   integer, parameter, public:: DFNTF_VAX = 2
   integer, parameter, public:: DFNTF_CRAY = 3
   integer, parameter, public:: DFNTF_PC = 4
   integer, parameter, public:: DFNTF_CONVEX = 5
   integer, parameter, public:: DFNTF_VP = 6

!       Masks for types
   integer, parameter, public:: DFNT_HDF = 0
   integer, parameter, public:: DFNT_NATIVE = 4096
   integer, parameter, public:: DFNT_CUSTOM = 8192
   integer, parameter, public:: DFNT_LITEND = 16384

!	Number type info codes 
   integer, parameter, public:: DFNT_NONE = 0
   integer, parameter, public:: DFNT_QUERY = 0
   integer, parameter, public:: DFNT_VERSION = 1
   
   integer, parameter, public:: DFNT_FLOAT32 = 5
   integer, parameter, public:: DFNT_FLOAT = 5
   integer, parameter, public:: DFNT_FLOAT64 = 6
   integer, parameter, public:: DFNT_DOUBLE = 6
   integer, parameter, public:: DFNT_FLOAT128 = 7
 
   integer, parameter, public:: DFNT_INT8 = 20
   integer, parameter, public:: DFNT_UINT8 = 21
   integer, parameter, public:: DFNT_INT16 = 22
   integer, parameter, public:: DFNT_UINT16 = 23
   integer, parameter, public:: DFNT_INT32 = 24
   integer, parameter, public:: DFNT_UINT32 = 25
   integer, parameter, public:: DFNT_INT64 = 26
   integer, parameter, public:: DFNT_UINT64 = 27
   integer, parameter, public:: DFNT_INT128 = 28
   integer, parameter, public:: DFNT_UINT128 = 29
 
   integer, parameter, public:: DFNT_UCHAR8 = 3
   integer, parameter, public:: DFNT_UCHAR = 3
   integer, parameter, public:: DFNT_CHAR8 = 4
   integer, parameter, public:: DFNT_CHAR = 4
   integer, parameter, public:: DFNT_CHAR16 = 42
   integer, parameter, public:: DFNT_UCHAR16 = 43

   integer, parameter, public:: DFNT_NFLOAT32 = 4101
   integer, parameter, public:: DFNT_NFLOAT = 4101
   integer, parameter, public:: DFNT_NFLOAT64 = 4102
   integer, parameter, public:: DFNT_NDOUBLE = 4102
   integer, parameter, public:: DFNT_NFLOAT128 = 4103
 
   integer, parameter, public:: DFNT_NINT8 = 4116
   integer, parameter, public:: DFNT_NUINT8 = 4117
   integer, parameter, public:: DFNT_NINT16 = 4118
   integer, parameter, public:: DFNT_NUINT16 = 4119
   integer, parameter, public:: DFNT_NINT32 = 4120
   integer, parameter, public:: DFNT_NUINT32 = 4121
   integer, parameter, public:: DFNT_NINT64 = 4122
   integer, parameter, public:: DFNT_NUINT64 = 4123
   integer, parameter, public:: DFNT_NINT128 = 4124
   integer, parameter, public:: DFNT_NUINT128 = 4125
 
   integer, parameter, public:: DFNT_NUCHAR8 = 4099
   integer, parameter, public:: DFNT_NUCHAR = 4099
   integer, parameter, public:: DFNT_NCHAR8 = 4100
   integer, parameter, public:: DFNT_NCHAR = 4100
   integer, parameter, public:: DFNT_NCHAR16 = 4138
   integer, parameter, public:: DFNT_NUCHAR16 = 4139

   integer, parameter, public:: DFNT_LFLOAT32 = 16389
   integer, parameter, public:: DFNT_LFLOAT = 16389
   integer, parameter, public:: DFNT_LFLOAT64 = 16390
   integer, parameter, public:: DFNT_LDOUBLE = 16390
   integer, parameter, public:: DFNT_LFLOAT128 = 16391
 
   integer, parameter, public:: DFNT_LINT8 = 16404
   integer, parameter, public:: DFNT_LUINT8 = 16405
   integer, parameter, public:: DFNT_LINT16 = 16406
   integer, parameter, public:: DFNT_LUINT16 = 16407
   integer, parameter, public:: DFNT_LINT32 = 16408
   integer, parameter, public:: DFNT_LUINT32 = 16409
   integer, parameter, public:: DFNT_LINT64 = 16410
   integer, parameter, public:: DFNT_LUINT64 = 16411
   integer, parameter, public:: DFNT_LINT128 = 16412
   integer, parameter, public:: DFNT_LUINT128 = 16413
 
   integer, parameter, public:: DFNT_LUCHAR8 = 16387
   integer, parameter, public:: DFNT_LUCHAR = 16387
   integer, parameter, public:: DFNT_LCHAR8 = 16388
   integer, parameter, public:: DFNT_LCHAR = 16388
   integer, parameter, public:: DFNT_LCHAR16 = 16426
   integer, parameter, public:: DFNT_LUCHAR16 = 16427

!	tags and refs
   integer, parameter, public:: DFREF_WILDCARD = 0, DFTAG_WILDCARD = 0
   integer, parameter, public:: DFTAG_NULL = 1, DFTAG_LINKED = 20
   integer, parameter, public:: DFTAG_VERSION = 30,DFTAG_COMPRESSED = 40


!	utility set
   integer, parameter, public:: DFTAG_FID = 100, DFTAG_FD = 101
   integer, parameter, public:: DFTAG_TID = 102, DFTAG_TD = 103
   integer, parameter, public:: DFTAG_DIL = 104, DFTAG_DIA = 105
   integer, parameter, public:: DFTAG_NT = 106, DFTAG_MT = 107

! 	raster-8 set 
   integer, parameter, public:: DFTAG_ID8 = 200, DFTAG_IP8 = 201
   integer, parameter, public:: DFTAG_RI8 = 202, DFTAG_CI8 = 203
   integer, parameter, public:: DFTAG_II8 = 204

!	Raster Image set
   integer, parameter, public:: DFTAG_ID = 300, DFTAG_LUT = 301
   integer, parameter, public:: DFTAG_RI = 302, DFTAG_CI = 303
   integer, parameter, public:: DFTAG_RIG = 306, DFTAG_LD = 307
   integer, parameter, public:: DFTAG_MD = 308, DFTAG_MA = 309
   integer, parameter, public:: DFTAG_CCN = 310, DFTAG_CFM = 311
   integer, parameter, public:: DFTAG_AR = 312
   integer, parameter, public:: DFTAG_DRAW = 400, DFTAG_RUN = 401
   integer, parameter, public:: DFTAG_XYP = 500, DFTAG_MTO = 501

!	Tektronix 
   integer, parameter, public:: DFTAG_T14 = 602, DFTAG_T105 = 603

!	Scientific Data set 
   integer, parameter, public:: DFTAG_SDG = 700, DFTAG_SDD = 701
   integer, parameter, public:: DFTAG_SD = 702, DFTAG_SDS = 703
   integer, parameter, public:: DFTAG_SDL = 704, DFTAG_SDU = 705
   integer, parameter, public:: DFTAG_SDF = 706, DFTAG_SDM = 707
   integer, parameter, public:: DFTAG_SDC = 708, DFTAG_SDT = 709
   integer, parameter, public:: DFTAG_SDLNK = 710, DFTAG_NDG = 720
   integer, parameter, public:: DFTAG_CAL = 731, DFTAG_FV = 732
   integer, parameter, public:: DFTAG_BREQ = 799, DFTAG_EREQ = 780

!	VSets 
   integer, parameter, public:: DFTAG_VG = 1965, DFTAG_VH = 1962
   integer, parameter, public:: DFTAG_VS = 1963

!	compression schemes 
   integer, parameter, public:: DFTAG_RLE = 11, DFTAG_IMC = 12
   integer, parameter, public:: DFTAG_IMCOMP = 12, DFTAG_JPEG = 13
   integer, parameter, public:: DFTAG_GREYJPEG = 14

!	SPECIAL CODES 
   integer, parameter, public:: SPECIAL_LINKED = 1, SPECIAL_EXT = 2

!	PARAMETERS 
   integer, parameter, public:: DF_MAXFNLEN = 256, SD_UNLIMITED = 0
   integer, parameter, public:: SD_DIMVAL_BW_COMP = 1, SD_DIMVAL_BW_INCOMP = 0
   integer, parameter, public:: SD_FILL = 0, SD_NOFILL = 256
   integer, parameter, public:: HDF_VDATA = -1

!       Standard return codes       
   integer, parameter, public:: SUCCEED = 0, FAIL = -1


!	Compression Types 
   integer, parameter, public:: COMP_NONE = 0, COMP_RLE = 11
   integer, parameter, public:: COMP_IMCOMP = 12, COMP_JPEG = 2

!       Fortran chunking (SD and GR interfaces) and compression routines use
!       the following compression types:
   integer, parameter, public:: COMP_CODE_NONE = 0
   integer, parameter, public:: COMP_CODE_RLE = 1
   integer, parameter, public:: COMP_CODE_NBIT = 2
   integer, parameter, public:: COMP_CODE_SKPHUFF = 3
   integer, parameter, public:: COMP_CODE_DEFLATE = 4
   integer, parameter, public:: COMP_CODE_JPEG = 6

!	Interlace Types 
   integer, parameter, public:: MFGR_INTERLACE_PIXEL = 0
   integer, parameter, public:: MFGR_INTERLACE_LINE = 1
   integer, parameter, public:: MFGR_INTERLACE_COMPONENT = 2

   integer, parameter, public:: FULL_INTERLACE = 0, NO_INTERLACE = 1

!       Vdata fields packing types
   integer, parameter, public:: HDF_VSPACK = 0, HDF_VSUNPACK = 1

!    Multi-file Annotation types
   integer, parameter, public:: AN_DATA_LABEL = 0, AN_DATA_DESC = 1
   integer, parameter, public:: AN_FILE_LABEL = 2, AN_FILE_DESC = 3
   
!    Maximum number of dimensions - mpav
   integer, parameter, public:: MAX_RANK_HDF = 10
end module HDF
