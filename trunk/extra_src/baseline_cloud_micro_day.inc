!----------------------------------------------------------------------
!
! Name:
!   baseline_cloud_micro_day.inc
!
! Description:
!   Include file that defines structures for reflectance and 
!   emissivity tables
!
!----------------------------------------------------------------------
 Integer(kind=int4), parameter, private:: Chan_Idx_065um = 1
 Integer(kind=int4), parameter, private:: Chan_Idx_160um = 6
 Integer(kind=int4), parameter, private:: Chan_Idx_375um = 20

!---- fixed thresholds
 Real(kind=real4), parameter, private:: zen_max = 70.0
 Real(kind=real4), parameter, private:: solzen_max = 70.0
 Real(kind=real4), parameter, private:: pp_err = 0.30    
 Real(kind=real4), parameter, private:: cal_err = 0.05   

 Integer, parameter, private:: nphases = 3

!-----------------------------------------------------------------------
! structures for reflectance tables
!-----------------------------------------------------------------------
   Type, private :: Reflectance_Lookup_Table
     Integer (kind=int4) :: flag              !- flag whether populated or not
     character (len=72) :: header             !- descriptive string
     Integer (kind=int4) :: nlog10tau         !- size of tau dimension
     Integer (kind=int4) :: nlog10re          !- size of re dimension
     Integer (kind=int4) :: nzen           !- size of zenith dimension
     Integer (kind=int4) :: nsolzen        !- size of solar zenith dimension
     Integer (kind=int4) :: nrelaz            !- size of relative az. dimension
     Integer (kind=int4) :: nrelaz_in_pp      !- size of relative az. dimension in pp zone
     Integer (kind=int4) :: nrelaz_out_pp     !- size of relative az. dimension out pp zone
     Real (kind=real4) :: delta_log10tau      !- spacing between tau points
     Real (kind=real4) :: delta_log10re       !- spacing between re points
     Real (kind=real4) :: delta_zen           !- spacing between zen points
     Real (kind=real4) :: delta_solzen        !- spacing between solzen points
     Real (kind=real4) :: delta_relaz_in_pp   !- spacing between relaz points in pp zone
     Real (kind=real4) :: delta_relaz_out_pp  !- spacing between relaz points out pp zone
     Real (kind=real4) :: relaz_switch_pp     !- spacing between relaz points out pp zone
     Real (kind=real4) :: lambda              !- wavelength of this channel
     Real (kind=real4) :: lambda_ref          !- reference wavelength for optical depth
     Real (kind=real4) :: rho                 !- density
     Real (kind=real4), dimension(:), allocatable :: zen       !zenith vector
     Real (kind=real4), dimension(:), allocatable :: solzen    !solar zenith vector
     Real (kind=real4), dimension(:), allocatable :: relaz     !relative azimuth vector
     Real (kind=real4), dimension(:), allocatable :: log10tau  !log10 optical depth vector
     Real (kind=real4), dimension(:), allocatable:: log10re    !log10 effective radius vector
     Real (kind=real4), dimension(:), allocatable :: ext_ratio !extinction ratio (lambda/lambda_ref)
     Real (kind=real4), dimension(:), allocatable :: Qext_ref  !extinction efficiency at lambda_ref
     Real (kind=real4), dimension(:,:), allocatable :: sph_alb !spherical albedo
     Real (kind=real4), dimension(:,:,:), allocatable :: trn   !flux transmission
     Real (kind=real4), dimension(:,:,:), allocatable :: alb   !flux albedo
     Real (kind=real4), dimension(:,:,:,:,:), allocatable :: ref  !reflectance
  End Type Reflectance_Lookup_Table

  Type, private:: Phase_Reflectance_Tables
      Integer (kind=int1):: flag 
      Type(Reflectance_Lookup_Table),dimension(nphases):: phase
  End Type Phase_Reflectance_Tables

  Type, private:: Channel_Reflectance_Tables
      Integer (kind=int1):: flag
      Type(Phase_Reflectance_Tables), dimension(nchan_max):: channel
  End Type Channel_Reflectance_Tables

  Type (Channel_Reflectance_Tables), private, save, dimension(:), allocatable:: cld_ref


!-----------------------------------------------------------------------
! structures for emissivity tables
!-----------------------------------------------------------------------
   Type, private :: Emissivity_Lookup_Table
     Integer (kind=int4) :: flag
     character (len=72) :: header
     Integer (kind=int4) :: nlog10tau
     Integer (kind=int4) :: nlog10re
     Integer (kind=int4) :: nzen
     Real (kind=real4) :: delta_log10tau
     Real (kind=real4) :: delta_log10re
     Real (kind=real4) :: delta_zen
     Real (kind=real4) :: lambda
     Real (kind=real4) :: lambda_ref
     Real (kind=real4) :: rho
     Real (kind=real4), dimension(:), allocatable :: zen 
     Real (kind=real4), dimension(:), allocatable :: log10tau 
     Real (kind=real4), dimension(:), allocatable :: log10re 
     Real (kind=real4), dimension(:), allocatable :: Qext_ref 
     Real (kind=real4), dimension(:), allocatable :: ext_ratio 
     Real (kind=real4), dimension(:,:,:), allocatable :: trans 
     Real (kind=real4), dimension(:,:,:), allocatable :: emiss 
  End Type Emissivity_Lookup_Table

  Type, private:: Phase_Emissivity_Tables
      Integer (kind=int1):: flag
      Type(Emissivity_Lookup_Table),dimension(nphases):: phase
  End Type Phase_Emissivity_Tables

  Type, private:: Channel_Emissivity_Tables
      Integer (kind=int1):: flag
      Type(Phase_Emissivity_Tables), dimension(nchan_max):: channel
  End Type Channel_Emissivity_Tables

  Type (Channel_Emissivity_Tables), private, save, dimension(:), allocatable:: cld_ems




