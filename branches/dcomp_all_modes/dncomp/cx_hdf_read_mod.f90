!  started to make this icar tool ( thanks) in one file and a bit better readable
!  AW
!   


module cx_hdf_read_mod
   implicit none
   private
   
   include 'hdf.f90'
   include 'cx_hdf_standard_declarations.inc'
   
   integer, parameter :: MAXNCDIM = 32
   integer, parameter :: MAXNCNAM = 128
   
   integer, parameter :: UNCLBRTD = -1
   integer, parameter ::DECLBRTD =  0
   integer, parameter :: CALIBRTD =  1
      
   type hdf_data
      integer :: type
      integer :: datasize
      integer :: size
      integer :: nval
      integer :: rank
      integer :: dimsize(MAXNCDIM)
      integer :: utype
      integer :: calbrtd
      
      real ( kind = 8 ) :: calibr(4)
      character, dimension(:), allocatable :: c1values
      integer(kind=1), dimension(:), allocatable :: i1values
      integer(kind=2), dimension(:), allocatable :: i2values
      integer(kind=4), dimension(:), allocatable :: i4values
      real(kind=4), dimension(:), allocatable :: r4values
      real(kind=8), dimension(:), allocatable :: r8values
      
   end type hdf_data
   
   type hdf_att
      character ( len = MAXNCNAM) :: name
      type(hdf_data) :: data
   end type hdf_att
   
   type hdf_sds
      character( len = MAXNCNAM) :: name
      integer :: nattr
      type ( hdf_data) :: data
      type(hdf_att), dimension(:),allocatable :: attr
   end type hdf_sds
   
   
   public :: hdf_data
   public :: hdf_sds
   public :: MAXNCDIM
   public :: MAXNCNAM
   public :: hdf_get_file_sds
   
contains
   !
   !
   !  
   function hdf_get_finfo(hdf_file, nsds, sds_name, natt, att_name)

      integer :: hdf_get_finfo     
      
      
      character(len=*), intent(in) :: hdf_file
      integer, intent(out) :: nsds
      integer, intent(out) :: natt                 
      character ( len = MAXNCNAM), intent(out), allocatable :: sds_name(:)
      character ( len = MAXNCNAM), intent(out), allocatable :: att_name(:)
   
      integer :: sd_id
      integer :: sds_id
      integer :: ierr
      integer :: iobj
      integer :: dummy (MAXNCDIM)  

      hdf_get_finfo = -1

      sd_id = sfstart(hdf_file, DFACC_READ)
      if (hdf_error(sd_id)) return
    
      if (hdf_error(sffinfo(sd_id, nsds, natt))) goto 99999
    
      if (nsds > 0) then

         if (allocated(sds_name)) then
            deallocate(sds_name, stat=ierr)
            if (std_error(ierr /= 0, "Dynamic memory desallocation error 1")) goto 99999
         end if
         
         allocate (sds_name(1:nsds), stat=ierr)
         if (std_error(ierr /= 0, "Dynamic memory allocation error")) goto 99999

         do iobj = 1, nsds
            sds_id = sfselect(sd_id, iobj - 1)
            if (hdf_error(sds_id)) goto 99999
            if (hdf_error(sfginfo(sds_id, sds_name(iobj), dummy, dummy, dummy, dummy))) goto 99999
         end do

      end if

      if (natt > 0) then

         if (allocated(att_name)) then
            deallocate(att_name, stat=ierr)
            if (std_error(ierr /= 0, "Dynamic memory desallocation error")) goto 99999
         end if
         
         allocate (att_name(1:natt), stat=ierr)
         if (std_error(ierr /= 0, "Dynamic memory allocation error")) goto 99999

         do iobj = 1, natt
            if (hdf_error(sfgainfo(sd_id, iobj - 1, att_name(iobj), dummy, dummy))) goto 99999
         end do

      end if

      hdf_get_finfo = 0

99999 continue

      if (hdf_error(sfend(sd_id))) hdf_get_finfo = -1

   end function hdf_get_finfo
  
   !-------------------------------------------------------------------------------
   !
   !
   !
   function hdf_get_file_att(hdf_file, natt, attrs, nattn, att_name)

      integer :: hdf_get_file_att
      character(len=*), intent(in) :: hdf_file
      integer, intent(out) :: natt
      type(hdf_att), intent(out), allocatable :: attrs(:)
      integer, intent(in), optional:: nattn  
      character(len=*), intent(in), optional :: att_name(*)
      
      integer :: sd_id
      integer :: nsds
    
      hdf_get_file_att = -1

      sd_id = sfstart(hdf_file, DFACC_READ)
      if (hdf_error(sd_id)) return

      if (hdf_error(sffinfo(sd_id, nsds, natt))) goto 99999
   
      if (natt > 0) then
         if (hdf_get_obj_att(sd_id, natt, attrs, nattn, att_name) < 0) goto 99999
      end if

      hdf_get_file_att = 0
   
99999 continue

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

    end do                          

    hdf_get_obj_att = 0
    
99999 continue

    

    if (allocated(attind)) then
       deallocate (attind, stat=ierr)
       if (std_error(ierr /= 0, "Dynamic memory desallocation error"))  hdf_get_obj_att = -1
    endif

  end function hdf_get_obj_att
  
   !----------------------------------------------------------------------------------------

   function hdf_get_file_sds(hdf_file, nsds, sdata, nsdsn, sds_name, dclb, attr, cal_sub, outtype)

      integer :: hdf_get_file_sds

      character (len=*), intent( in) :: hdf_file
      integer, intent(out)  :: nsds
      type(hdf_sds), intent(out),  allocatable, target :: sdata(:)
      integer, optional, intent( in) :: nsdsn
      character (len=*), intent( in), optional :: sds_name(*)

      logical, optional :: dclb
      logical, optional :: attr
       
      integer, optional, intent( in) :: outtype
      optional :: cal_sub
      external cal_sub

      integer :: sd_id           
      integer :: sds_id          
      integer :: start(MAXNCDIM)
      integer :: stride(MAXNCDIM)
      integer :: ierr            
      integer :: isds            
      integer :: idim            
      integer :: att_idx         
      integer :: rtype          
      integer :: natt
   
      logical:: att_sw  
      integer, allocatable :: sdsind(:)
      integer*1, allocatable :: bdata(:)
      real(kind=8), allocatable :: r8data(:)
      type(hdf_sds), pointer :: ps
      type(hdf_data), pointer :: pd  

      hdf_get_file_sds = -1

      if (std_error((present(nsdsn).or.present(sds_name)).and.(.not.(present(nsdsn).and.present(sds_name))), &
         "Optionnal arguments must be both defined or undefined")) return

    
      sd_id = sfstart(hdf_file, DFACC_READ)
      if (hdf_error(sd_id)) return

      if (hdf_error(sffinfo(sd_id, nsds, natt))) goto 99999

      if (present(nsdsn)) nsds = nsdsn
      att_sw = .true.
      if (present(attr)) att_sw = attr

      if (nsds <= 0) goto 99999

      allocate (sdsind(1:nsds), stat=ierr)
      if (std_error(ierr /= 0, "Dynamic memory allocation error")) goto 99999

      if (.not.present(nsdsn)) then
         do isds = 1, nsds
            sdsind(isds) = isds - 1
         end do
      else
         do isds = 1, nsds
            sdsind(isds) = sfn2index(sd_id, sds_name(isds))
            if (hdf_error(sdsind(isds))) goto 99999
         end do
      end if

      allocate (sdata(1:nsds), stat=ierr)
      if (std_error(ierr /= 0, "Dynamic memory allocation error")) goto 99999

      do isds = 1, nsds        

         sds_id = sfselect(sd_id, sdsind(isds))
         if (hdf_error(sds_id)) goto 99999

          ps => sdata(isds)
          pd => sdata(isds)%data

 
         if (hdf_error(sfginfo(sds_id, ps%name, pd%rank, pd%dimsize, pd%type, ps%nattr))) &
                   & goto 99999


         pd % calbrtd = DECLBRTD
         
         if (sfgcal(sds_id, pd%calibr(1), pd%calibr(2), pd%calibr(3), pd%calibr(4), pd%utype) < 0) then
            att_idx = sffattr(sds_id, 'scale_factor')
            if (att_idx >= 0) then
               att_idx = sffattr(sds_id, 'add_offset')
               if (std_error(att_idx >= 0, 'Error while getting calibration coefficients')) goto 99999
            end if
            
            pd%calibr = 0.; pd%calibr(1) = 1.
            pd%calbrtd = UNCLBRTD
            pd%utype = pd%type
         
         else if (present(dclb)) then
            if (.not.dclb) pd%calbrtd = CALIBRTD
         end if

         if (pd%calbrtd == DECLBRTD) then
            rtype = pd%utype;
         else
            rtype = pd%type;
         endif
       
         stride(:) = 1
         start (:) = 0

         pd % datasize = hdf_typesize(rtype)
         pd % nval = 1
         do idim = 1, pd%rank
            pd%nval = pd%nval*pd%dimsize(idim)
         end do

         if (present(outtype)) then
            pd%nval = pd%nval*hdf_typesize(rtype)/hdf_typesize(outtype)
            pd%datasize = hdf_typesize(outtype)
            pd%type = outtype
            rtype = outtype
         end if

         pd % size = pd%nval*pd%datasize
       
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
            end if

            if (allocated(bdata)) then
               deallocate(bdata, stat=ierr)
               if (std_error(ierr /= 0, "Dynamic memory desallocation error")) goto 99999
            endif

         end if


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


       if (allocated(r8data)) then
          deallocate(r8data, stat=ierr)
          if (std_error(ierr /= 0, "Dynamic memory desallocation error")) goto 99999
       endif

       ! Récupération des éventuels attributs du SDS si demandée
       if (att_sw .and. sdata(isds)%nattr > 0) then
          if (hdf_get_obj_att(sds_id, sdata(isds)%nattr, sdata(isds)%attr) < 0) goto 99999
       endif

      
       if (hdf_error(sfendacc(sds_id))) goto 99999

    end do                           !

    hdf_get_file_sds = 0

99999 continue

    
    if (allocated(sdsind)) then
       deallocate (sdsind, stat=ierr)
       if (std_error(ierr /= 0, "Dynamic memory desallocation error")) hdf_get_file_sds = -1
    endif

   
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
  
   function hdf_get_sds(hdf_file, sds_name)

      character (len=*), intent( in)                   :: &
         hdf_file               

      character (len=*), intent( in)                   :: &
         sds_name              
      integer*1, dimension(:), pointer                 :: &
         hdf_get_sds           

      type(hdf_sds), dimension(:), allocatable, target :: &
         sds                   
      integer                                          :: &
         nsds                  

      nullify(hdf_get_sds)

      if (hdf_get_file_sds(hdf_file, nsds, sds, 1, (/sds_name/), .false., .false., outtype=DFNT_UCHAR8) >= 0) then

         allocate(hdf_get_sds(sds(1)%data%size))
         hdf_get_sds = sds(1)%data%i1values

      endif

      call dealloc_sds(nsds, sds)

  end function hdf_get_sds

   !-------------------------------------------------------------------------------

   function hdf_get_int_sds(hdf_file, sds_name, raw, cal_sub)

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

       elseif (allocated(sds(1)%data%r4values)) then
          hdf_get_dbl_sds_3D = reshape(sds(1)%data%r4values, (/dim1, dim2, dim3/))
       elseif (allocated(sds(1)%data%r8values)) then
          hdf_get_dbl_sds_3D = reshape(sds(1)%data%r8values, (/dim1, dim2, dim3/))
       endif
      
    endif

    call dealloc_sds(nsds, sds)

  end function hdf_get_dbl_sds_3D

end module cx_hdf_read_mod
