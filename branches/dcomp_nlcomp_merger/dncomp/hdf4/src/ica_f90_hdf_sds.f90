module ica_f90_hdf_sds
  
  include 'icahdfsdsinc.f90'
  
contains
  
!-------------------------------------------------------------------------------
  
  function hdf_get_finfo(hdf_file, nsds, sds_name, natt, att_name)

    integer :: hdf_get_finfo      ! 0 si OK, -1 sinon

    ! Déclaration des arguments de la function

    character(len=*), intent(in)                                     :: &
         hdf_file                 ! Chemin du fichier HDF

    integer, intent(out)                                             :: &
         nsds,                  & ! Nombre de SDS du fichier
         natt                     ! Nombre d'attributs du fichier

    character (len=MAXNCNAM), intent(out), dimension(:), allocatable :: &
         sds_name,              & ! Noms des SDS du fichier
         att_name                 ! Noms des attributs du fichier
    
    ! Déclaration des variables locales

    integer                                                          :: &
         sd_id,                 & ! Identificateur HDF du fichier
         sds_id,                & ! Identificateur HDF d'un SDS
         ierr,                  & ! Code d'erreur
         iobj,                  & ! Indice de boucle
         dummy(MAXNCDIM)          ! Pour arguments inutilisés

    hdf_get_finfo = -1

    ! Ouverture du fichier HDF

    sd_id = sfstart(hdf_file, DFACC_READ)
    if (hdf_error(sd_id)) return
    
    ! Récupération du nombre de SDS et d'attributs du fichier

    if (hdf_error(sffinfo(sd_id, nsds, natt))) goto 99999
    
    if (nsds > 0) then

       ! Allocation dynamique du tableau des noms de SDS

       if (allocated(sds_name)) then
          deallocate(sds_name, stat=ierr)
          if (std_error(ierr /= 0, "Dynamic memory desallocation error 1")) goto 99999
       endif
       allocate (sds_name(1:nsds), stat=ierr)
       if (std_error(ierr /= 0, "Dynamic memory allocation error")) goto 99999

       ! Récupération des noms des SDS

       do iobj = 1, nsds
          sds_id = sfselect(sd_id, iobj - 1)
          if (hdf_error(sds_id)) goto 99999
          if (hdf_error(sfginfo(sds_id, sds_name(iobj), dummy, dummy, dummy, dummy))) goto 99999
       enddo

    endif

    if (natt > 0) then
 
       ! Allocation dynamique du tableau des noms d'attributs

       if (allocated(att_name)) then
          deallocate(att_name, stat=ierr)
          if (std_error(ierr /= 0, "Dynamic memory desallocation error")) goto 99999
       endif
       allocate (att_name(1:natt), stat=ierr)
       if (std_error(ierr /= 0, "Dynamic memory allocation error")) goto 99999

       ! Récupération des noms des SDS

       do iobj = 1, natt
          if (hdf_error(sfgainfo(sd_id, iobj - 1, att_name(iobj), dummy, dummy))) goto 99999
       enddo

    endif

    hdf_get_finfo = 0

99999 continue

    ! Fermeture du fichier HDF

    if (hdf_error(sfend(sd_id))) hdf_get_finfo = -1

  end function hdf_get_finfo
  
!-------------------------------------------------------------------------------

  function hdf_get_file_att(hdf_file, natt, attrs, nattn, att_name)

    integer :: hdf_get_file_att

    ! Déclaration des arguments de la function

    character(len=*), intent(in)                          :: &
         hdf_file               ! Chemin du fichier HDF

    integer, intent(out)                                  :: &
         natt                   ! Nombre d'attributs réellement retournés du fichier

    type(hdf_att), intent(out), dimension(:), allocatable :: &
         attrs                  ! Tableau des structures d'attributs du fichier

    integer, intent(in), optional                         :: &
         nattn                  ! Nombre (optionnel) d'attributs requis
    
    character(len=*), intent(in), optional                :: &
         att_name(*)            ! Tableau (optionnel) des noms des attributs requis

    ! Déclaration des variables locales

    integer                                               :: &
         sd_id,               & ! Identificateur HDF du fichier
         nsds                   ! Nombre de SDS du fichier
    
    hdf_get_file_att = -1

    ! Ouverture du fichier HDF

    sd_id = sfstart(hdf_file, DFACC_READ)
    if (hdf_error(sd_id)) return
    
    ! Récupération du nombre d'attributs du fichier

    if (hdf_error(sffinfo(sd_id, nsds, natt))) goto 99999
   
    ! Récupération des attributs du fichier

    if (natt > 0) then
       if (hdf_get_obj_att(sd_id, natt, attrs, nattn, att_name) < 0) goto 99999
    endif

    hdf_get_file_att = 0
   
99999 continue

    ! Fermeture du fichier HDF

    if (hdf_error(sfend(sd_id))) hdf_get_file_att = -1

  end function hdf_get_file_att
  
!-------------------------------------------------------------------------------

  function hdf_get_obj_att(obj_id, natt, attrs, nattn, att_name)
    
    integer :: hdf_get_obj_att

    ! Déclaration des arguments de la function

    integer, intent(in)                                   :: &
         obj_id                 ! Identificateur HDF de l'objet
    
    integer, intent(out)                                  :: &
         natt                   ! Nombre d'attributs réellement retournés

    type(hdf_att), intent(out), dimension(:), allocatable :: &
         attrs                  ! Tableau des structures d'attributs de l'objet
    
    integer, intent( in), optional                        :: &
         nattn                  ! Nombre (optionnel) d'attributs requis
    
    character(len=*), intent( in), optional               :: &
         att_name(*)            ! Tableau (optionnel) des noms des attributs requis
    
    ! Déclaration des variables locales

    integer                                               :: &
         iatt,                & ! Indice de l'attribut en cours
         ierr                   ! Code d'erreur
 
    integer, dimension(:), allocatable                    :: &
         attind                 ! Table des index des attributs à extraire

    hdf_get_obj_att = -1

    ! Contrôle de validité des paramètres optionnels

!    if (std_error((present(nattn).or.present(att_name)).ne.(present(nattn).and.present(att_name)),  &
    if (std_error((present(nattn).or.present(att_name)).and.(.not.(present(nattn).and.present(att_name))),  &
         "Optionnal arguments must be both defined or undefined")) return

    if (present(nattn)) natt = nattn
    
    ! On regarde s'il y a quelque-chose à faire

    if (natt <= 0) return

    ! Allocation de la table des index des attributs à extraire

    allocate (attind(1:natt), stat=ierr)
    if (std_error(ierr /= 0, "Dynamic memory allocation error")) return

    ! Remplissage de la table des index des attributs à extraire

    if (.not.present(nattn)) then
       do iatt = 1, natt
          attind(iatt) = iatt - 1
       enddo
    else
       natt = nattn
       do iatt = 1, natt
          attind(iatt) = sffattr(obj_id, att_name(iatt))
          if (hdf_error(attind(iatt))) goto 99999
       enddo
    endif
    
    ! Allocation du tableau des structures d'attributs

    if (allocated(attrs)) then
       deallocate (attrs, stat=ierr)
       if (std_error(ierr /= 0, "Dynamic memory desallocation error")) goto 99999
    endif
    allocate (attrs(1:natt), stat=ierr)
    if (std_error(ierr /= 0, "Dynamic memory allocation error")) goto 99999

    do iatt = 1, natt                ! DÉBUT DE BOUCLE SUR LES ATTRIBUTS À EXTRAIRE

       ! Récupération des infos sur l'attribut

       if (hdf_error(sfgainfo(obj_id, attind(iatt), attrs(iatt)%name, attrs(iatt)%data%type, &
            attrs(iatt)%data%dimsize))) goto 99999

       ! On renseigne les infos sur l'attribut dans la structure

       attrs(iatt)%data%rank = 1
       attrs(iatt)%data%datasize = hdf_typesize(attrs(iatt)%data%type)
       attrs(iatt)%data%size = attrs(iatt)%data%dimsize(1)*attrs(iatt)%data%datasize
       attrs(iatt)%data%nval = attrs(iatt)%data%dimsize(1)

       ! Allocation et remplissage des valeurs de l'attribut selon son type

       select case (attrs(iatt)%data%type)

       case (DFNT_CHAR8)
          allocate (attrs(iatt)%data%c1values(1:attrs(iatt)%data%dimsize(1)), stat=ierr)
          if (std_error(ierr /= 0, "Dynamic memory allocation error")) goto 99999
          if (hdf_error(sfrcatt(obj_id, iatt - 1, attrs(iatt)%data%c1values))) goto 99999

       case (DFNT_UCHAR8, DFNT_UINT8, DFNT_INT8)
          allocate (attrs(iatt)%data%i1values(1:attrs(iatt)%data%dimsize(1)), stat=ierr)
          if (std_error(ierr /= 0, "Dynamic memory allocation error")) goto 99999
          if (hdf_error(sfrcatt(obj_id, iatt - 1, attrs(iatt)%data%i1values))) goto 99999

       case (DFNT_UINT16, DFNT_INT16)
          allocate (attrs(iatt)%data%i2values(1:attrs(iatt)%data%dimsize(1)), stat=ierr)
          if (std_error(ierr /= 0, "Dynamic memory allocation error")) goto 99999
          if (hdf_error(sfrcatt(obj_id, iatt - 1, attrs(iatt)%data%i2values))) goto 99999

       case (DFNT_UINT32, DFNT_INT32)
          allocate (attrs(iatt)%data%i4values(1:attrs(iatt)%data%dimsize(1)), stat=ierr)
          if (std_error(ierr /= 0, "Dynamic memory allocation error")) goto 99999
          if (hdf_error(sfrcatt(obj_id, iatt - 1, attrs(iatt)%data%i4values))) goto 99999

! En prévision...
!       case (DFNT_UINT64, DFNT_INT64)
!          allocate (attrs(iatt)%data%i8values(1:attrs(iatt)%data%dimsize(1)), stat=ierr)
!          if (std_error(ierr /= 0, "Dynamic memory allocation error")) goto 99999
!          if (hdf_error(sfrcatt(obj_id, iatt - 1, attrs(iatt)%data%i8values))) goto 99999

       case (DFNT_FLOAT32)
          allocate (attrs(iatt)%data%r4values(1:attrs(iatt)%data%dimsize(1)), stat=ierr)
          if (std_error(ierr /= 0, "Dynamic memory allocation error")) goto 99999
          if (hdf_error(sfrcatt(obj_id, iatt - 1, attrs(iatt)%data%r4values))) goto 99999

       case (DFNT_FLOAT64)
          allocate (attrs(iatt)%data%r8values(1:attrs(iatt)%data%dimsize(1)), stat=ierr)
          if (std_error(ierr /= 0, "Dynamic memory allocation error")) goto 99999
          if (hdf_error(sfrcatt(obj_id, iatt - 1, attrs(iatt)%data%r8values))) goto 99999

       case default
          print *, "Unimplemented data type: ", attrs(iatt)%data%type; goto 99999

       end select

    end do                           ! FIN DE BOUCLE SUR LES ATTRIBUTS À EXTRAIRE

    hdf_get_obj_att = 0
    
99999 continue

    ! Libération de la table d'index

    if (allocated(attind)) then
       deallocate (attind, stat=ierr)
       if (std_error(ierr /= 0, "Dynamic memory desallocation error"))  hdf_get_obj_att = -1
    endif

  end function hdf_get_obj_att
  
!----------------------------------------------------------------------------------------

  function hdf_get_file_sds(hdf_file, nsds, sdata, nsdsn, sds_name, dclb, attr, cal_sub, outtype)

    integer :: hdf_get_file_sds

    ! Déclaration des arguments de la function

    character (len=*), intent( in)                                :: &
         hdf_file               ! Chemin du fichier HDF

    integer, intent(out)                                          :: &
         nsds                   ! Nombre de SDS effectivement extraits

    type(hdf_sds), intent(out), dimension(:), allocatable, target :: &
         sdata                  ! Tableau des SDS extraits

    integer, optional, intent( in)                                :: &
         nsdsn                  ! Nombre (optionnel) de SDS requis

    character (len=*), intent( in), optional                      :: &
         sds_name(*)            ! Tableau (optionnel) des noms des SDS requis

    logical, optional :: &
         dclb, &
         attr

    integer, optional, intent( in)                                :: &
         outtype                ! Type éventuellement imposé des données du SDS en sortie

    optional :: &
         cal_sub

    external cal_sub

    ! Déclaration des variables locales

    integer                                                      :: &
         sd_id,               & ! Identificateur HDF du fichier
         sds_id,              & ! Identificateur HDF du SDS courant
         start(MAXNCDIM),     & ! Table des indices de début pour l'extraction d'un SDS
         stride(MAXNCDIM),    & ! Table des pas par dimension pour l'extraction d'un SDS
         ierr,                & ! Code d'erreur
         isds,                & ! Compteur sur les SDS à extraire
         idim,                & ! Compteur sur les dimensions du SDS courant
         att_idx,             &
         rtype,               & ! Type réel des données extraites (suivant décalibration)
         natt

    logical                                                       :: &
         att_sw

    integer, dimension(:), allocatable                            :: &
         sdsind                 ! Table des index des attributs à extraire

    integer*1, dimension(:), allocatable                          :: &
         bdata

    real(kind=8), dimension(:), allocatable                       :: &
         r8data

    type(hdf_sds), pointer                                        :: &
         ps

    type(hdf_data), pointer                                       :: &
         pd

    hdf_get_file_sds = -1

    ! Contrôle de validité des paramètres optionnels

    if (std_error((present(nsdsn).or.present(sds_name)).and.(.not.(present(nsdsn).and.present(sds_name))), &
         "Optionnal arguments must be both defined or undefined")) return

    ! Ouverture du fichier HDF

    sd_id = sfstart(hdf_file, DFACC_READ)
    if (hdf_error(sd_id)) return

    if (hdf_error(sffinfo(sd_id, nsds, natt))) goto 99999

    if (present(nsdsn)) nsds = nsdsn
    att_sw = .true.
    if (present(attr)) att_sw = attr

    ! On regarde s'il y a quelque-chose à faire

    if (nsds <= 0) goto 99999

    ! Allocation de la table des index des SDS à extraire

    allocate (sdsind(1:nsds), stat=ierr)
    if (std_error(ierr /= 0, "Dynamic memory allocation error")) goto 99999

    ! Remplissage de la table des index des SDS à extraire

    if (.not.present(nsdsn)) then
       do isds = 1, nsds
          sdsind(isds) = isds - 1
       enddo
    else
       do isds = 1, nsds
          sdsind(isds) = sfn2index(sd_id, sds_name(isds))
          if (hdf_error(sdsind(isds))) goto 99999
       enddo
    endif

    ! Allocation du tableau des SDS

    allocate (sdata(1:nsds), stat=ierr)
    if (std_error(ierr /= 0, "Dynamic memory allocation error")) goto 99999

    do isds = 1, nsds                ! DEBUT DE BOUCLE SUR LES SDS A EXTRAIRE

       ! Initialisation de l'accès au SDS

       sds_id = sfselect(sd_id, sdsind(isds))
       if (hdf_error(sds_id)) goto 99999

       ! Création de pointeurs pour le confort d'écriture

       ps => sdata(isds); pd => sdata(isds)%data

       ! Récupération d'informations sur le SDS

       if (hdf_error(sfginfo(sds_id, ps%name, pd%rank, pd%dimsize, pd%type, ps%nattr))) &
            goto 99999

       ! Récupération des éventuels paramètres de calibration des données du SDS.
       ! Il est impossible de faire la distinction entre erreur HDF et absence de
       ! coefficients de calibration (dans les 2 cas, le code de retour de sfgcal
       ! vaut -1). Si cela se produit, on essaye d'atteindre les coeff. de
       ! calibration en tant qu'attributs par leur nom. Si on les trouve, alors
       ! c'est que sfgcal a retourné une erreur et on fait la même chose; sinon,
       ! c'est qu'il n'y a pas de calibration et on initialise les coefficients
       ! respectifs à 1 et 0

       pd%calbrtd = DECLBRTD
       if (sfgcal(sds_id, pd%calibr(1), pd%calibr(2), pd%calibr(3), pd%calibr(4), pd%utype) < 0) then
          att_idx = sffattr(sds_id, 'scale_factor')
          if (att_idx >= 0) then
             att_idx = sffattr(sds_id, 'add_offset')
             if (std_error(att_idx >= 0, 'Error while getting calibration coefficients')) goto 99999
          endif
          pd%calibr = 0.; pd%calibr(1) = 1.
          pd%calbrtd = UNCLBRTD
          pd%utype = pd%type
       elseif (present(dclb)) then
          if (.not.dclb) pd%calbrtd = CALIBRTD
       endif

       if (pd%calbrtd == DECLBRTD) then
          rtype = pd%utype;
       else
          rtype = pd%type;
       endif
       
       ! Allocation dynamique et/ou initialisation des variables nécessaires à l'extraction:
       ! - remplissage à 0 du tableau start et à 1 du tableau stride
       ! - allocation de la zone pd après calcul de sa taille

       stride(:) = 1; start (:) = 0

       pd%datasize = hdf_typesize(rtype)
       pd%nval = 1
       do idim = 1, pd%rank
          pd%nval = pd%nval*pd%dimsize(idim)
       end do

       if (present(outtype)) then
          pd%nval = pd%nval*hdf_typesize(rtype)/hdf_typesize(outtype)
          pd%datasize = hdf_typesize(outtype)
          pd%type = outtype
          rtype = outtype
       endif

       pd%size = pd%nval*pd%datasize
       
       if (pd%calbrtd == DECLBRTD) then
          if (allocated(bdata)) then
             deallocate(bdata, stat=ierr)
             if (std_error(ierr /= 0, "Dynamic memory desallocation error")) goto 99999
          endif
          allocate(bdata(pd%size), stat=ierr)
          if (std_error(ierr /= 0, "Dynamic memory allocation error")) goto 99999

          if (hdf_error(sfrdata(sds_id, start, stride, pd%dimsize, bdata))) goto 99999

          if (allocated(r8data)) then
             deallocate(r8data, stat=ierr)
             if (std_error(ierr /= 0, "Dynamic memory desallocation error")) goto 99999
          endif
          allocate(r8data(pd%nval), stat=ierr)
          if (std_error(ierr /= 0, "Dynamic memory allocation error")) goto 99999

          select case (pd%type)
          case (DFNT_CHAR8)
             print *, "Character data cannot be decalibrated"; goto 99999
          case (DFNT_UCHAR8, DFNT_UINT8, DFNT_INT8)
             r8data = transfer(bdata, pd%i1values)
          case (DFNT_UINT16, DFNT_INT16)
             r8data = transfer(bdata, pd%i2values)
          case (DFNT_UINT32, DFNT_INT32)
             r8data = transfer(bdata, pd%i4values)
! En prévision...
!          case (DFNT_UINT64, DFNT_INT64)
!             r8data = transfer(bdata, pd%i8values)
          case (DFNT_FLOAT32)
             r8data = transfer(bdata, pd%r4values)
          case (DFNT_FLOAT64)
             r8data = transfer(bdata, pd%r8values)
          case default
             print *, "Unimplemented data type: ", pd%type; goto 99999
          end select

          if (present(cal_sub)) then
             call cal_sub(pd%nval, r8data, pd%calibr(1), pd%calibr(3))
          else
             call default_dclb(pd%nval, r8data, pd%calibr(1), pd%calibr(3))
          endif

          if (allocated(bdata)) then
             deallocate(bdata, stat=ierr)
             if (std_error(ierr /= 0, "Dynamic memory desallocation error")) goto 99999
          endif

       endif

       ! Allocation et remplissage des valeurs du SDS selon son type

       select case (rtype)

       case (DFNT_CHAR8)
          allocate (pd%c1values(1:pd%nval), stat=ierr)
          if (std_error(ierr /= 0, "Dynamic memory allocation error")) goto 99999
          if (pd%calbrtd == DECLBRTD) then
             pd%i1values = r8data
          else
             if (hdf_error(sfrdata(sds_id, start, stride, pd%dimsize, pd%c1values))) goto 99999
          endif

       case (DFNT_UCHAR8, DFNT_UINT8, DFNT_INT8)
          allocate (pd%i1values(1:pd%nval), stat=ierr)
          if (std_error(ierr /= 0, "Dynamic memory allocation error")) goto 99999
          if (pd%calbrtd == DECLBRTD) then
             pd%i1values = r8data
          else
             if (hdf_error(sfrdata(sds_id, start, stride, pd%dimsize, pd%i1values))) goto 99999
          endif
          
       case (DFNT_UINT16, DFNT_INT16)
          allocate (pd%i2values(1:pd%nval), stat=ierr)
          if (std_error(ierr /= 0, "Dynamic memory allocation error")) goto 99999
          if (pd%calbrtd == DECLBRTD) then
             pd%i2values = r8data
          else
             if (hdf_error(sfrdata(sds_id, start, stride, pd%dimsize, pd%i2values))) goto 99999
          endif
          
       case (DFNT_UINT32, DFNT_INT32)
          allocate (pd%i4values(1:pd%nval), stat=ierr)
          if (std_error(ierr /= 0, "Dynamic memory allocation error")) goto 99999
          if (pd%calbrtd == DECLBRTD) then
             pd%i4values = r8data
          else
             if (hdf_error(sfrdata(sds_id, start, stride, pd%dimsize, pd%i4values))) goto 99999
          endif
             
! En prévision...
!       case (DFNT_UINT64, DFNT_INT64)
!          allocate (pd%i8values(1:pd%nval), stat=ierr)
!          if (std_error(ierr /= 0, "Dynamic memory allocation error")) goto 99999
!          if (pd%calbrtd == DECLBRTD) then
!             pd%i8values = r8data
!          else
!             if (hdf_error(sfrdata(sds_id, start, stride, pd%dimsize, pd%i8values))) goto 99999
!          endif
             
       case (DFNT_FLOAT32)
          allocate (pd%r4values(1:pd%nval), stat=ierr)
          if (std_error(ierr /= 0, "Dynamic memory allocation error")) goto 99999
          if (pd%calbrtd == DECLBRTD) then
             pd%r4values = r8data
          else
             if (hdf_error(sfrdata(sds_id, start, stride, pd%dimsize, pd%r4values))) goto 99999
          endif
          
       case (DFNT_FLOAT64)
          allocate (pd%r8values(1:pd%nval), stat=ierr)
          if (std_error(ierr /= 0, "Dynamic memory allocation error")) goto 99999
          if (pd%calbrtd == DECLBRTD) then
             pd%r8values = r8data
          else
             if (hdf_error(sfrdata(sds_id, start, stride, pd%dimsize, pd%r8values))) goto 99999
          endif

       case default
          print *, "Unimplemented data type: ", pd%type; goto 99999

       end select

!          select case (pd%type)
!          case (DFNT_CHAR8)
!             print *, "Character data cannot be decalibrated"; goto 99999
!          case (DFNT_UCHAR8, DFNT_UINT8, DFNT_INT8)
!             r8data = (transfer(bdata, pd%i1values) - pd%calibr(3))*pd%calibr(1)
!          case (DFNT_UINT16, DFNT_INT16)
!             r8data = (transfer(bdata, pd%i2values) - pd%calibr(3))*pd%calibr(1)
!          case (DFNT_UINT32, DFNT_INT32)
!             r8data = (transfer(bdata, pd%i4values) - pd%calibr(3))*pd%calibr(1)
!! En prévision...
!!          case (DFNT_UINT64, DFNT_INT64)
!!             r8data = (transfer(bdata, pd%i8values) - pd%calibr(3))*pd%calibr(1)
!          case (DFNT_FLOAT32)
!             r8data = (transfer(bdata, pd%r4values) - pd%calibr(3))*pd%calibr(1)
!          case (DFNT_FLOAT64)
!             r8data = (transfer(bdata, pd%r8values) - pd%calibr(3))*pd%calibr(1)
!          case default
!             print *, "Unimplemented data type: ", pd%type; goto 99999
!          end select
!
!       ! Allocation et remplissage des valeurs du SDS selon son type
!
!       select case (rtype)
!
!       case (DFNT_CHAR8)
!          allocate (pd%c1values(1:pd%nval), stat=ierr)
!          if (std_error(ierr /= 0, "Dynamic memory allocation error")) goto 99999
!          if (pd%calbrtd == DECLBRTD) then
!             pd%i1values = r8data
!          else
!             if (hdf_error(sfrdata(sds_id, start, stride, pd%dimsize, pd%c1values))) goto 99999
!          endif
!
!       case (DFNT_UCHAR8, DFNT_UINT8, DFNT_INT8)
!          allocate (pd%i1values(1:pd%nval), stat=ierr)
!          if (std_error(ierr /= 0, "Dynamic memory allocation error")) goto 99999
!          if (pd%calbrtd == DECLBRTD) then
!             pd%i1values = r8data
!          else
!             if (hdf_error(sfrdata(sds_id, start, stride, pd%dimsize, pd%i1values))) goto 99999
!          endif
!          
!       case (DFNT_UINT16, DFNT_INT16)
!          allocate (pd%i2values(1:pd%nval), stat=ierr)
!          if (std_error(ierr /= 0, "Dynamic memory allocation error")) goto 99999
!          if (pd%calbrtd == DECLBRTD) then
!             pd%i2values = r8data
!          else
!             if (hdf_error(sfrdata(sds_id, start, stride, pd%dimsize, pd%i2values))) goto 99999
!          endif
!          
!       case (DFNT_UINT32, DFNT_INT32)
!          allocate (pd%i4values(1:pd%nval), stat=ierr)
!          if (std_error(ierr /= 0, "Dynamic memory allocation error")) goto 99999
!          if (pd%calbrtd == DECLBRTD) then
!             pd%i4values = r8data
!          else
!             if (hdf_error(sfrdata(sds_id, start, stride, pd%dimsize, pd%i4values))) goto 99999
!          endif
!             
!! En prévision...
!!       case (DFNT_UINT64, DFNT_INT64)
!!          allocate (pd%i8values(1:pd%nval), stat=ierr)
!!          if (std_error(ierr /= 0, "Dynamic memory allocation error")) goto 99999
!!          if (pd%calbrtd == DECLBRTD) then
!!             pd%i8values = r8data
!!          else
!!             if (hdf_error(sfrdata(sds_id, start, stride, pd%dimsize, pd%i8values))) goto 99999
!!          endif
!             
!       case (DFNT_FLOAT32)
!          allocate (pd%r4values(1:pd%nval), stat=ierr)
!          if (std_error(ierr /= 0, "Dynamic memory allocation error")) goto 99999
!          if (pd%calbrtd == DECLBRTD) then
!             pd%r4values = r8data
!          else
!             if (hdf_error(sfrdata(sds_id, start, stride, pd%dimsize, pd%r4values))) goto 99999
!          endif
!          
!       case (DFNT_FLOAT64)
!          allocate (pd%r8values(1:pd%nval), stat=ierr)
!          if (std_error(ierr /= 0, "Dynamic memory allocation error")) goto 99999
!          if (pd%calbrtd == DECLBRTD) then
!             pd%r8values = r8data
!          else
!             if (hdf_error(sfrdata(sds_id, start, stride, pd%dimsize, pd%r8values))) goto 99999
!          endif
!
!       case default
!          print *, "Unimplemented data type: ", pd%type; goto 99999
!
!       end select

       if (allocated(r8data)) then
          deallocate(r8data, stat=ierr)
          if (std_error(ierr /= 0, "Dynamic memory desallocation error")) goto 99999
       endif

       ! Récupération des éventuels attributs du SDS si demandée
       if (att_sw .and. sdata(isds)%nattr > 0) then
          if (hdf_get_obj_att(sds_id, sdata(isds)%nattr, sdata(isds)%attr) < 0) goto 99999
       endif

       ! Fin d'accès au SDS
       if (hdf_error(sfendacc(sds_id))) goto 99999

    end do                           ! FIN DE BOUCLE SUR LES SDS A EXTRAIRE

    hdf_get_file_sds = 0

99999 continue

    ! Libération mémoire devenue inutile
    if (allocated(sdsind)) then
       deallocate (sdsind, stat=ierr)
       if (std_error(ierr /= 0, "Dynamic memory desallocation error")) hdf_get_file_sds = -1
    endif

    ! Fermeture du fichier HDF
    if (hdf_error(sfend(sd_id))) hdf_get_file_sds = -1

  end function hdf_get_file_sds

!-------------------------------------------------------------------------------
  
  subroutine default_dclb(n, data, pente, offset)

    integer :: n
    real(kind=8), dimension(n) :: data
    real(kind=8) :: pente, offset

    data = (data - offset)*pente

  end subroutine default_dclb

!-------------------------------------------------------------------------------
  
  function hdf_error(hdfcode)
    logical :: hdf_error
    integer :: hdfcode, dummy
    hdf_error = (hdfcode < 0)
    if (hdf_error) dummy = heprnt(0)
    return
  end function hdf_error

!-------------------------------------------------------------------------------
  
  function std_error(ierr, msg)
    logical :: std_error, ierr
    character (len=*) :: msg
    std_error = ierr
    if (ierr) print *, '***ERROR: '//msg
  end function std_error

!-------------------------------------------------------------------------------
  
  function hdf_typesize(t)
    
    ! La fonction hdf_typesize retourne la taille en octets d'un type HDF passé en argument

    integer            :: &
         hdf_typesize,    & ! Taille retournée
         t                  ! Type de donnée HDF
  
    select case (t)
    case (DFNT_UCHAR8, DFNT_CHAR8, DFNT_UINT8, DFNT_INT8); hdf_typesize = 1
    case (DFNT_UINT16, DFNT_INT16)                       ; hdf_typesize = 2
    case (DFNT_UINT32, DFNT_INT32, DFNT_FLOAT32)         ; hdf_typesize = 4
    case (DFNT_FLOAT64)                                  ; hdf_typesize = 8
    case default                                         ; hdf_typesize = 0
    end select
  
  end function hdf_typesize

!-------------------------------------------------------------------------------

  function hdf_typedesc(htyp)

    ! La fonction hdf_typedesc donne la description d'un type de donnée HDF passé en argument.
    ! Si le type de donnée est invalide, 'unknow data type' est retourné.

    integer           :: &
         htyp            ! (IN)  Type de donnée HDF

    character(len=60) :: &
         hdf_typedesc    ! (OUT) Description du type de donnée HDF

     select case (htyp)
     case (DFNT_UCHAR8 ); hdf_typedesc = '8-bit unsigned character / integer*1'
     case (DFNT_CHAR8  ); hdf_typedesc = '8-bit signed character / character*1'
     case (DFNT_UINT8  ); hdf_typedesc = '8-bit unsigned integer / not supported'
     case (DFNT_INT8   ); hdf_typedesc = '8-bit signed integer / integer*1'
     case (DFNT_UINT16 ); hdf_typedesc = '16-bit unsigned integer / not supported'
     case (DFNT_INT16  ); hdf_typedesc = '16-bit signed integer / integer*2'
     case (DFNT_UINT32 ); hdf_typedesc = '32-bit unsigned integer / not supported'
     case (DFNT_INT32  ); hdf_typedesc = '32-bit signed integer / integer*4'
     case (DFNT_UINT64 ); hdf_typedesc = '64-bit unsigned integer / not supported'
     case (DFNT_INT64  ); hdf_typedesc = '64-bit signed integer / integer*4'
     case (DFNT_FLOAT32); hdf_typedesc = '32-bit floating point number / real*4'
     case (DFNT_FLOAT64); hdf_typedesc = '64-bit floating point number / real*8'
     case default       ; hdf_typedesc = 'unsupported data type'
     end select

     return
   end function hdf_typedesc

!-------------------------------------------------------------------------------

  subroutine dealloc_sds(nsds, sds)

    integer :: isds, nsds, iatt

    type(hdf_sds), dimension(:), allocatable :: &
         sds                  ! Tableau des SDS extraits

    do isds = 1, nsds
       if (allocated(sds(isds)%attr)) then
          do iatt = 1, sds(isds)%nattr
             call dealloc_hdata(sds(isds)%attr(iatt)%data)
          enddo
          deallocate(sds(isds)%attr)
       endif
       call dealloc_hdata(sds(isds)%data)
    enddo
    if (allocated(sds)) deallocate(sds)

  end subroutine dealloc_sds

!-------------------------------------------------------------------------------

  subroutine dealloc_hdata(data)

    type(hdf_data) :: &
         data           

    if (allocated(data%c1values)) deallocate(data%c1values)
    if (allocated(data%i1values)) deallocate(data%i1values)
    if (allocated(data%i2values)) deallocate(data%i2values)
    if (allocated(data%i4values)) deallocate(data%i4values)
    !if (allocated(data%i8values)) deallocate(data%i8values)
    if (allocated(data%r4values)) deallocate(data%r4values)
    if (allocated(data%r8values)) deallocate(data%r8values)

  end subroutine dealloc_hdata

!-------------------------------------------------------------------------------

end module ica_f90_hdf_sds
