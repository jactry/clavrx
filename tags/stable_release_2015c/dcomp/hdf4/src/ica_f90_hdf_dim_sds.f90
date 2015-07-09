module ica_f90_hdf_dim_sds
  
  use ica_f90_hdf_sds
  
contains
  
!-------------------------------------------------------------------------------
  
  function hdf_get_sds(hdf_file, sds_name)

    ! Déclaration des arguments de la function

    character (len=*), intent( in)                   :: &
         hdf_file               ! Chemin du fichier HDF

    character (len=*), intent( in)                   :: &
         sds_name               ! Nom du SDS à extraire

    integer*1, dimension(:), pointer                 :: &
         hdf_get_sds            ! Valeur retournée par la fonction

    ! Déclaration des variables locales

    type(hdf_sds), dimension(:), allocatable, target :: &
         sds                    ! Tableau des SDS extraits

    integer                                          :: &
         nsds                   ! Nombre de SDS effectivement extraits

    nullify(hdf_get_sds)

    if (hdf_get_file_sds(hdf_file, nsds, sds, 1, (/sds_name/), .false., .false., outtype=DFNT_UCHAR8) >= 0) then

       allocate(hdf_get_sds(sds(1)%data%size))
       hdf_get_sds = sds(1)%data%i1values

    endif

    call dealloc_sds(nsds, sds)

  end function hdf_get_sds

!-------------------------------------------------------------------------------

  function hdf_get_int_sds(hdf_file, sds_name, raw, cal_sub)

    ! Déclaration des arguments de la function

    character (len=*), intent( in)                   :: &
         hdf_file               ! Chemin du fichier HDF

    character (len=*), intent( in)                   :: &
         sds_name               ! Chemin du fichier HDF

    logical, optional                                :: &
         raw                    ! true si données brutes, .false. ou omis sinon

    integer(kind=4), dimension(:), pointer           :: &
         hdf_get_int_sds        ! Valeur retournée par la fonction

    optional :: &
         cal_sub

    external cal_sub

    ! Déclaration des variables locales

    type(hdf_sds), dimension(:), allocatable, target :: &
         sds                    ! Tableau des SDS extraits

    integer                                          :: &
         nsds                   ! Nombre de SDS effectivement extraits

    logical                                          :: &
         dclb                   ! Flag réel de décalibration

    dclb = .true.
    if (present(raw)) dclb = (.not.raw)
    nullify(hdf_get_int_sds)

    if (hdf_get_file_sds(hdf_file, nsds, sds, 1, (/sds_name/), dclb, .false., cal_sub=cal_sub) >= 0) then

       allocate(hdf_get_int_sds(sds(1)%data%dimsize(1)))

       if     (allocated(sds(1)%data%c1values)) then
          hdf_get_int_sds = transfer(sds(1)%data%c1values, (/1,1/))
       elseif (allocated(sds(1)%data%i1values)) then
          hdf_get_int_sds = sds(1)%data%i1values
       elseif (allocated(sds(1)%data%i2values)) then
          hdf_get_int_sds = sds(1)%data%i2values
       elseif (allocated(sds(1)%data%i4values)) then
          hdf_get_int_sds = sds(1)%data%i4values
! En prévision...
!       elseif (allocated(sds(1)%data%i8values)) then
!          hdf_get_int_sds = sds(1)%data%i8values
       elseif (allocated(sds(1)%data%r4values)) then
          hdf_get_int_sds = sds(1)%data%r4values
       elseif (allocated(sds(1)%data%r8values)) then
          hdf_get_int_sds = sds(1)%data%r8values
       endif
       
    endif

    call dealloc_sds(nsds, sds)

  end function hdf_get_int_sds

!-------------------------------------------------------------------------------
  
  function hdf_get_int_sds_2D(hdf_file, sds_name, raw, cal_sub)

    ! Déclaration des arguments de la function

    character (len=*), intent( in)                   :: &
         hdf_file               ! Chemin du fichier HDF

    character (len=*), intent( in)                   :: &
         sds_name               ! Chemin du fichier HDF

    logical, optional                                :: &
         raw                    ! true si données brutes, .false. ou omis sinon

    integer(kind=4), dimension(:,:), pointer         :: &
         hdf_get_int_sds_2D     ! Valeur retournée par la fonction

    optional :: &
         cal_sub

    external cal_sub

    ! Déclaration des variables locales

    type(hdf_sds), dimension(:), allocatable, target :: &
         sds                    ! Tableau des SDS extraits

    integer                                          :: &
         nsds,                & ! Nombre de SDS effectivement extraits
         dim1,                & ! 1ère dimension du SDS
         dim2                   ! 2ème dimension du SDS

    logical                                          :: &
         dclb                   ! Flag réel de décalibration

    dclb = .true.
    if (present(raw)) dclb = (.not.raw)
    nullify(hdf_get_int_sds_2D)

    if (hdf_get_file_sds(hdf_file, nsds, sds, 1, (/sds_name/), dclb, .false., cal_sub=cal_sub) >= 0) then

       dim1 = sds(1)%data%dimsize(1)
       dim2 = sds(1)%data%dimsize(2)
       allocate(hdf_get_int_sds_2D(dim1, dim2))

       if     (allocated(sds(1)%data%c1values)) then
          hdf_get_int_sds_2D = reshape(transfer(sds(1)%data%c1values, (/1,1/)), (/dim1, dim2/))
       elseif (allocated(sds(1)%data%i1values)) then
          hdf_get_int_sds_2D = reshape(sds(1)%data%i1values, (/dim1, dim2/))
       elseif (allocated(sds(1)%data%i2values)) then
          hdf_get_int_sds_2D = reshape(sds(1)%data%i2values, (/dim1, dim2/))
       elseif (allocated(sds(1)%data%i4values)) then
          hdf_get_int_sds_2D = reshape(sds(1)%data%i4values, (/dim1, dim2/))
! En prévision...
!       elseif (allocated(sds(1)%data%i8values)) then
!          hdf_get_int_sds_2D = reshape(sds(1)%data%i8values, (/dim1, dim2/))
       elseif (allocated(sds(1)%data%r4values)) then
          hdf_get_int_sds_2D = reshape(sds(1)%data%r4values, (/dim1, dim2/))
       elseif (allocated(sds(1)%data%r8values)) then
          hdf_get_int_sds_2D = reshape(sds(1)%data%r8values, (/dim1, dim2/))
       endif
      
    endif
       
    call dealloc_sds(nsds, sds)

  end function hdf_get_int_sds_2D

!-------------------------------------------------------------------------------
    
  function hdf_get_int_sds_3D(hdf_file, sds_name, raw, cal_sub)

    ! Déclaration des arguments de la function

    character (len=*), intent( in)                   :: &
         hdf_file               ! Chemin du fichier HDF

    character (len=*), intent( in)                   :: &
         sds_name               ! Chemin du fichier HDF

    logical, optional                                :: &
         raw                    ! true si données brutes, .false. ou omis sinon

    integer(kind=4), dimension(:,:,:), pointer       :: &
         hdf_get_int_sds_3D     ! Valeur retournée par la fonction

    optional :: &
         cal_sub

    external cal_sub

    ! Déclaration des variables locales

    type(hdf_sds), dimension(:), allocatable, target :: &
         sds                    ! Tableau des SDS extraits

    integer                                          :: &
         nsds,                & ! Nombre de SDS effectivement extraits
         dim1,                & ! 1ère dimension du SDS
         dim2,                & ! 2ème dimension du SDS
         dim3                   ! 3ème dimension du SDS

    logical                                          :: &
         dclb                   ! Flag réel de décalibration

    dclb = .true.
    if (present(raw)) dclb = (.not.raw)
    nullify(hdf_get_int_sds_3D)

    if (hdf_get_file_sds(hdf_file, nsds, sds, 1, (/sds_name/), dclb, .false., cal_sub=cal_sub) >= 0) then

       dim1 = sds(1)%data%dimsize(1)
       dim2 = sds(1)%data%dimsize(2)
       dim3 = sds(1)%data%dimsize(3)
       allocate(hdf_get_int_sds_3D(dim1, dim2, dim3))

       if     (allocated(sds(1)%data%c1values)) then
          hdf_get_int_sds_3D = reshape(transfer(sds(1)%data%c1values, (/1,1/)), (/dim1, dim2, dim3/))
       elseif (allocated(sds(1)%data%i1values)) then
          hdf_get_int_sds_3D = reshape(sds(1)%data%i1values, (/dim1, dim2, dim3/))
       elseif (allocated(sds(1)%data%i2values)) then
          hdf_get_int_sds_3D = reshape(sds(1)%data%i2values, (/dim1, dim2, dim3/))
       elseif (allocated(sds(1)%data%i4values)) then
          hdf_get_int_sds_3D = reshape(sds(1)%data%i4values, (/dim1, dim2, dim3/))
! En prévision...
!       elseif (allocated(sds(1)%data%i8values)) then
!          hdf_get_int_sds_3D = reshape(sds(1)%data%i8values, (/dim1, dim2, dim3/))
       elseif (allocated(sds(1)%data%r4values)) then
          hdf_get_int_sds_3D = reshape(sds(1)%data%r4values, (/dim1, dim2, dim3/))
       elseif (allocated(sds(1)%data%r8values)) then
          hdf_get_int_sds_3D = reshape(sds(1)%data%r8values, (/dim1, dim2, dim3/))
       endif
      
    endif
       
    call dealloc_sds(nsds, sds)

  end function hdf_get_int_sds_3D

!-------------------------------------------------------------------------------
    
  function hdf_get_dbl_sds(hdf_file, sds_name, raw, cal_sub)

    ! Déclaration des arguments de la function

    character (len=*), intent( in)                   :: &
         hdf_file               ! Chemin du fichier HDF

    character (len=*), intent( in)                   :: &
         sds_name               ! Chemin du fichier HDF

    logical, optional                                :: &
         raw                    ! true si données brutes, .false. ou omis sinon

    real(kind=8), dimension(:), pointer              :: &
         hdf_get_dbl_sds        ! Valeur retournée par la fonction

    ! Déclaration des variables locales

    type(hdf_sds), dimension(:), allocatable, target :: &
         sds                    ! Tableau des SDS extraits

    integer                                          :: &
         nsds                   ! Nombre de SDS effectivement extraits

    logical                                          :: &
         dclb                   ! Flag réel de décalibration

    optional :: &
         cal_sub

    external cal_sub

    dclb = .true.
    if (present(raw)) dclb = (.not.raw)
    nullify(hdf_get_dbl_sds)

    if (hdf_get_file_sds(hdf_file, nsds, sds, 1, (/sds_name/), dclb, .false., cal_sub=cal_sub) >= 0) then

       allocate(hdf_get_dbl_sds(sds(1)%data%dimsize(1)))

       if     (allocated(sds(1)%data%c1values)) then
          hdf_get_dbl_sds = transfer(sds(1)%data%c1values, (/1,1/))
       elseif (allocated(sds(1)%data%i1values)) then
          hdf_get_dbl_sds = sds(1)%data%i1values
       elseif (allocated(sds(1)%data%i2values)) then
          hdf_get_dbl_sds = sds(1)%data%i2values
       elseif (allocated(sds(1)%data%i4values)) then
          hdf_get_dbl_sds = sds(1)%data%i4values
! En prévision...
!       elseif (allocated(sds(1)%data%i8values)) then
!          hdf_get_dbl_sds = sds(1)%data%i8values, (/dim1, dim2/))
       elseif (allocated(sds(1)%data%r4values)) then
          hdf_get_dbl_sds = sds(1)%data%r4values
       elseif (allocated(sds(1)%data%r8values)) then
          hdf_get_dbl_sds = sds(1)%data%r8values
       endif
      
    endif

    call dealloc_sds(nsds, sds)

  end function hdf_get_dbl_sds

!-------------------------------------------------------------------------------
    
  function hdf_get_dbl_sds_2D(hdf_file, sds_name, raw, cal_sub)

    ! Déclaration des arguments de la function

    character (len=*), intent( in)                   :: &
         hdf_file               ! Chemin du fichier HDF

    character (len=*), intent( in)                   :: &
         sds_name               ! Chemin du fichier HDF

    logical, optional                                :: &
         raw                    ! true si données brutes, .false. ou omis sinon

    real(kind=8), dimension(:,:), pointer            :: &
         hdf_get_dbl_sds_2D

    ! Déclaration des variables locales

    type(hdf_sds), dimension(:), allocatable, target :: &
         sds                    ! Tableau des SDS extraits

    integer                                          :: &
         nsds,                & ! Nombre de SDS effectivement extraits
         dim1,                & ! 1ère dimension du SDS
         dim2                   ! 2ème dimension du SDS

    logical                                          :: &
         dclb                   ! Flag réel de décalibration

    optional :: &
         cal_sub

    external cal_sub

    dclb = .true.
    if (present(raw)) dclb = (.not.raw)
    nullify(hdf_get_dbl_sds_2D)

    if (hdf_get_file_sds(hdf_file, nsds, sds, 1, (/sds_name/), dclb, .false., cal_sub=cal_sub) >= 0) then

       dim1 = sds(1)%data%dimsize(1)
       dim2 = sds(1)%data%dimsize(2)
       allocate(hdf_get_dbl_sds_2D(dim1, dim2))

       if     (allocated(sds(1)%data%c1values)) then
          hdf_get_dbl_sds_2D = reshape(transfer(sds(1)%data%c1values, (/1,1/)), (/dim1, dim2/))
       elseif (allocated(sds(1)%data%i1values)) then
          hdf_get_dbl_sds_2D = reshape(sds(1)%data%i1values, (/dim1, dim2/))
       elseif (allocated(sds(1)%data%i2values)) then
          hdf_get_dbl_sds_2D = reshape(sds(1)%data%i2values, (/dim1, dim2/))
       elseif (allocated(sds(1)%data%i4values)) then
          hdf_get_dbl_sds_2D = reshape(sds(1)%data%i4values, (/dim1, dim2/))
! En prévision...
!       elseif (allocated(sds(1)%data%i8values)) then
!          hdf_get_dbl_sds_2D = reshape(sds(1)%data%i8values, (/dim1, dim2/))
       elseif (allocated(sds(1)%data%r4values)) then
          hdf_get_dbl_sds_2D = reshape(sds(1)%data%r4values, (/dim1, dim2/))
       elseif (allocated(sds(1)%data%r8values)) then
          hdf_get_dbl_sds_2D = reshape(sds(1)%data%r8values, (/dim1, dim2/))
       endif
      
    endif

    call dealloc_sds(nsds, sds)

  end function hdf_get_dbl_sds_2D

!-------------------------------------------------------------------------------
    
  function hdf_get_dbl_sds_3D(hdf_file, sds_name, raw, cal_sub)

    ! Déclaration des arguments de la function

    character (len=*), intent( in)                   :: &
         hdf_file               ! Chemin du fichier HDF

    character (len=*), intent( in)                   :: &
         sds_name               ! Chemin du fichier HDF

    logical, optional                                :: &
         raw                    ! true si données brutes, .false. ou omis sinon

    real(kind=8), dimension(:,:,:), pointer          :: &
         hdf_get_dbl_sds_3D

    ! Déclaration des variables locales

    type(hdf_sds), dimension(:), allocatable, target :: &
         sds                    ! Tableau des SDS extraits

    integer                                          :: &
         nsds,                & ! Nombre de SDS effectivement extraits
         dim1,                & ! 1ère dimension du SDS
         dim2,                & ! 2ème dimension du SDS
         dim3                   ! 3ème dimension du SDS

    logical                                          :: &
         dclb                   ! Flag réel de décalibration

    optional :: &
         cal_sub

    external cal_sub

    dclb = .true.
    if (present(raw)) dclb = (.not.raw)
    nullify(hdf_get_dbl_sds_3D)

    if (hdf_get_file_sds(hdf_file, nsds, sds, 1, (/sds_name/), dclb, .false., cal_sub=cal_sub) >= 0) then

       dim1 = sds(1)%data%dimsize(1)
       dim2 = sds(1)%data%dimsize(2)
       dim3 = sds(1)%data%dimsize(3)
       allocate(hdf_get_dbl_sds_3D(dim1, dim2, dim3))

       if     (allocated(sds(1)%data%c1values)) then
          hdf_get_dbl_sds_3D = reshape(transfer(sds(1)%data%c1values, (/1,1/)), (/dim1, dim2, dim3/))
       elseif (allocated(sds(1)%data%i1values)) then
          hdf_get_dbl_sds_3D = reshape(sds(1)%data%i1values, (/dim1, dim2, dim3/))
       elseif (allocated(sds(1)%data%i2values)) then
          hdf_get_dbl_sds_3D = reshape(sds(1)%data%i2values, (/dim1, dim2, dim3/))
       elseif (allocated(sds(1)%data%i4values)) then
          hdf_get_dbl_sds_3D = reshape(sds(1)%data%i4values, (/dim1, dim2, dim3/))
! En prévision...
!       elseif (allocated(sds(1)%data%i8values)) then
!          hdf_get_dbl_sds_3D = reshape(sds(1)%data%i8values, (/dim1, dim2, dim3/))
       elseif (allocated(sds(1)%data%r4values)) then
          hdf_get_dbl_sds_3D = reshape(sds(1)%data%r4values, (/dim1, dim2, dim3/))
       elseif (allocated(sds(1)%data%r8values)) then
          hdf_get_dbl_sds_3D = reshape(sds(1)%data%r8values, (/dim1, dim2, dim3/))
       endif
      
    endif

    call dealloc_sds(nsds, sds)

  end function hdf_get_dbl_sds_3D

 end module ica_f90_hdf_dim_sds
