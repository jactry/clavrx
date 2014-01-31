! $Id$
! $URL:$

REAL(KIND=real4), PARAMETER, private :: SAT_ZEN_MAX = 82.0 !-
REAL(KIND=real4), PARAMETER, private :: SOL_ZEN_MAX = 82.0 !-

!----------------------------------------------------------------------
! optimal estimation fix values
!----------------------------------------------------------------------

REAL(KIND=real4), PARAMETER, PRIVATE :: CONV_CRIT   =  1.     !- convergence criteria
REAL(KIND=real4), PARAMETER, PRIVATE :: MAX_STEP_SIZE = 0.2   !- maximal size of delta_x
INTEGER, PARAMETER ,PRIVATE          :: MAX_ITER = 22         !- maximal number of iterations
REAL(kind=real4), PARAMETER, PRIVATE :: PP_ERR = 0.30	       !- assumed plan-parallel error
REAL(kind=real4), PARAMETER, PRIVATE :: CAL_ERR = 0.05	       !- assumed calibration error
REAL(kind=real4), PARAMETER, PRIVATE :: FM_ERR_WATER = 0.01	     !- assumed forward model error
REAL(kind=real4), PARAMETER, PRIVATE :: FM_ERR_ICE = 0.03	     !- assumed forward model error
REAL(kind=real4), PARAMETER, PRIVATE :: CPS_APRIORI_WATER = 1.0	 !- CPS aprori water 
REAL(kind=real4), PARAMETER, PRIVATE :: CPS_APRIORI_ICE = 1.3	     !- CPS aprori ice
REAL(kind=real4), PARAMETER, PRIVATE :: CPS_APRIORI_ERROR_WATER = 0.5
REAL(kind=real4), PARAMETER, PRIVATE :: CPS_APRIORI_ERROR_ICE = 0.75
REAL(kind=real4), PARAMETER, PRIVATE :: ERROR_OFFSET_VIS = 0.02	     !- CPS aprori ice
REAL(kind=real4), PARAMETER, PRIVATE :: ERROR_OFFSET_IR  = 0.02	     !- CPS aprori ice
REAL(kind=real4), PARAMETER, PRIVATE :: ERROR_SNOW_VIS  = 1000.	     !- Error for snow surface 
 
!-----------------------------------------------------------------------
!    atmospheric layer definitions ( below or above or all)
!-----------------------------------------------------------------------
 
INTEGER(KIND=int1),PARAMETER ,PRIVATE :: LAYER_BELOW_CLOUD = 1
INTEGER(KIND=int1),PARAMETER,PRIVATE :: LAYER_ABOVE_CLOUD = 2
INTEGER(KIND=int1),PARAMETER,PRIVATE :: LAYER_ALL_COLUMN = 3
  
!-----------------------------------------------------------------------
! phase definitions specific for this algorithm
!-----------------------------------------------------------------------

INTEGER(KIND=int1),PARAMETER,private :: WATER_PHASE = 1
INTEGER(KIND=int1),PARAMETER,private :: ICE_PHASE = 2

!-----------------------------------------------------------------------
! Aerosol defaults
!-----------------------------------------------------------------------

REAL(KIND=real4), PARAMETER, PRIVATE :: OD_AEROSOL_BACKGROUND_VIS = 0.01  ! optical depth
REAL(KIND=real4), PARAMETER, PRIVATE :: WO_AEROSOL_BACKGROUND_VIS = 0.9   ! single scattering albedo Zhao&Li 2006
REAL(KIND=real4), PARAMETER, PRIVATE :: G_AEROSOL_BACKGROUND_VIS  = 0.6   ! similarity PARAMETER
  
!-----------------------------------------------------------------------
! Rayleigh defaults
!-----------------------------------------------------------------------  

REAL(KIND=real4), PARAMETER, PRIVATE :: OD_RAYLEIGH_BACKGROUND = 0.044  ! optical depth

!-----------------------------------------------------------------------
! Quality Flags (to be done)
!----------------------------------------------------------------------- 

INTEGER(KIND=int4), PARAMETER, PRIVATE :: QC_CYCLE_VZA        = 2**0
INTEGER(KIND=int4), PARAMETER, PRIVATE :: QC_CYCLE_SZA        = 2**1
INTEGER(KIND=int4), PARAMETER, PRIVATE :: QC_CYCLE_NOCLOUD    = 2**2
INTEGER(KIND=int4), PARAMETER, PRIVATE :: QC_CYCLE_NOHEIGHT   = 2**3
INTEGER(KIND=int4), PARAMETER, PRIVATE :: QC_CYCLE_CLOUDTYPE  = 2**4
INTEGER(KIND=int4), PARAMETER, PRIVATE :: QC_CYCLE_SFC_PRESS  = 2**5
INTEGER(KIND=int4), PARAMETER, PRIVATE :: QC_CYCLE_NOCONVERGE = 2**6
INTEGER(KIND=int4), PARAMETER, PRIVATE :: QC_VALID_TWILIGHT   = 2**7
INTEGER(KIND=int4), PARAMETER, PRIVATE :: QC_CYCLE_WATER      = 2**8
INTEGER(KIND=int4), PARAMETER, PRIVATE :: QC_CYCLE_ICE        = 2**9
