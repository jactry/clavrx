  implicit none
  
  ! Includes nécessaires à l'utilisation des librairies HDF

  include 'icadffunc.f90'
  include 'hdf.f90'
 ! include 'netcdf.f90'
   
  integer MAXNCDIM 
  integer MAXNCNAM
   
  integer, parameter :: &
       UNCLBRTD = -1, &
       DECLBRTD =  0, &
       CALIBRTD =  1

  character(len=15), dimension(-1:1) :: &
       clbrstat = (/'SDS non calibré', 'SDS décalibré  ', 'SDS calibré    '/)

  ! Déclaration du nouveau type de donnée hdf_att
   
   parameter(MAXNCDIM = 32)
   parameter(MAXNCNAM = 128)
   
  type hdf_data

     integer                                    :: &
          type,                        & ! Type des données
          datasize,                    & ! Taille en octets du type des données
          size,                        & ! Taille en octets des données
          nval,                        & ! Nombre de valeurs des données
          rank,                        & ! Nombre de dimensions des données
          dimsize(MAXNCDIM),           & ! Valeurs des dimensions des données
          utype,                       & ! Type des données décalibrées
          calbrtd                        ! Statut de calibration des données

     real(kind=8) &
          calibr(4)                      ! Coefficients éventuels de calibration des données

     character,       dimension(:), allocatable :: &
          c1values                       ! Valeurs des données

     integer(kind=1), dimension(:), allocatable :: &
          i1values                       ! Valeurs des données

     integer(kind=2), dimension(:), allocatable :: &
          i2values                       ! Valeurs des données

     integer(kind=4), dimension(:), allocatable :: &
          i4values                       ! Valeurs des données

! En prévision...
!     integer(kind=8), dimension(:), allocatable :: &
!          i8values                       ! Valeurs des données

     real   (kind=4), dimension(:), allocatable :: &
          r4values                       ! Valeurs des données

     real   (kind=8), dimension(:), allocatable :: &
          r8values                       ! Valeurs des données

  end type hdf_data

  ! Déclaration du nouveau type de donnée hdf_att

  type hdf_att

     character (len=MAXNCNAM)                   :: &
          name                           ! Nom de l'attribut

     type(hdf_data)                             :: &
          data

  end type hdf_att

  ! Déclaration du nouveau type de donnée sds_struct

  type hdf_sds

     character (len=MAXNCNAM)                   :: &
          name                           ! Nom du SDS

     integer                                    :: &
          nattr                          ! Nombre d'attributs du SDS                     

     type(hdf_data)                             :: &
          data                           ! Données du SDS

     type(hdf_att), dimension(:), allocatable   :: &
          attr                           ! Attributs du SDS

  end type hdf_sds
