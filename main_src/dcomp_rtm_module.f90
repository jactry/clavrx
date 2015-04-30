!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: dcomp_rtm_module.f90 (src)
!       dcomp_rtm_module (program)
!
! PURPOSE: performs all dcomp related rtm calculations
!
! DESCRIPTION: 
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
! NOTES:
! -   performs all dcomp related rtm calculations 
! - 
! - 
! -   
!  -  output derived type:
!
! -  
!  - trans_ir_ac ( : , : )
!  - trans_ir_ac_nadir ( : , : )
! - tpw_ac ( : , : ) 
!  sfc_nwp ( : , : )  
!  rad_clear_sky_toc_ch20  ( : , : )
!  rad_clear_sky_toa_ch20 ( : , : )
!   ocone_path (:,:)
!--------------------------------------------------------------------------------------

module dcomp_rtm_module
   
   use constants, only: &
      INT1 , REAL4 , PI , INT4
   
      
   private
   public :: perform_rtm_dcomp
   public :: dcomp_rtm_type
   
   type dcomp_rtm_type
      real (kind= REAL4), allocatable :: trans_ir_ac ( : , : )
      real (kind= REAL4), allocatable :: trans_ir_ac_nadir ( : , : )
      real (kind= REAL4), allocatable :: tpw_ac ( : , : ) 
      real (kind= REAL4), allocatable :: sfc_nwp ( : , : )  
      real (kind= REAL4), allocatable :: rad_clear_sky_toc_ch20  ( : , : )
      real (kind= REAL4), allocatable :: rad_clear_sky_toa_ch20 ( : , : )
      real (kind= REAL4), allocatable :: ozone_path ( : , : )
      contains
      procedure :: deallocate_it
      
   end type dcomp_rtm_type
   
     
   
   contains
   
   subroutine  perform_rtm_dcomp (dcomp_rtm)
      
      use rtm_common, only: &
          nlevels_rtm &
          , rtm &
          , p_std_rtm 
          
      use pixel_common, only: &
           image &
         , geo &
         , acha &
         , i_nwp &
         , j_nwp &
         , zen_idx_rtm &
         , ch &
         , bad_pixel_mask &
         , Rad_Clear_Ch20_Solar_Rtm
        
      use nwp_common, only: &
           t_prof_nwp &
          , z_prof_nwp &
          , p_std_nwp &
          , sfc_level_nwp &
          , tropo_level_nwp &
          , inversion_level_nwp &
          , tpw_prof_nwp &
          , psfc_nwp &
          , ozone_nwp
      
      
      implicit none  
      type ( dcomp_rtm_type ) , intent (out) :: dcomp_rtm
      integer :: line_idx
      integer :: elem_idx
      integer :: dim_1 , dim_2
      
      integer :: x_nwp
      integer :: y_nwp
      
      real( kind = real4 )   :: placeholder_cld
      
      
      real(kind=real4)   :: cld_temp_loc
      real(kind=real4)   :: cld_press_loc
      real(kind=real4)   :: cld_height_loc
      
      integer( kind = int4 ) :: idx_lev_nwp
      integer( kind = int4 ) :: idx_lev_rtm
      real( kind = real4 ) :: prof_wgt_nwp
      real( kind = real4 ) :: prof_wgt_rtm
      
      integer( kind = int4 ) :: ivza
      
      integer(kind=int4), parameter :: num_levels_rtm_prof = nlevels_rtm

      real, dimension(num_levels_rtm_prof):: clear_trans_prof_rtm
      real, dimension(num_levels_rtm_prof):: clear_rad_prof_rtm
      
      
       ! -- pointers to nwp data structures 
      real( kind = real4 ), dimension(:), pointer :: temp_prof_nwp
      real( kind = real4 ), dimension(:), pointer :: hgt_prof_nwp
      real( kind = real4 ), dimension(:), pointer :: wv_prof_nwp 
      
      integer :: dim_outdata(2)
      
      ! -----------------------------------------
      
      
      dim_1 = Image%Number_Of_Elements
      dim_2 = Image%Number_Of_Lines_Read_This_Segment
      
           
      call allocate_dcomp_rtm ( dcomp_rtm , dim_1 , dim_2 )
     
       line_loop: do line_idx = 1 , dim_2   
         element_loop: do elem_idx= 1 ,dim_1
           
            if ( bad_pixel_mask(  elem_idx , line_idx ) == 1 ) cycle
            
            ! - alias local variables	
            cld_height_loc = acha % zc (elem_idx,line_idx)
            cld_temp_loc   = acha % tc (elem_idx,line_idx)
            cld_press_loc  = acha % pc (elem_idx,line_idx)
           
             ! - for convenience, save nwp indices to local variables
            x_nwp  =  i_nwp ( elem_idx , line_idx )       ! - nwp longitude cell     
            y_nwp  =  j_nwp ( elem_idx , line_idx )        ! - nwp latitude cell 
            
            dcomp_rtm % ozone_path (elem_idx,line_idx) = ozone_nwp(x_nwp,y_nwp)  
            
             ! - compute cloud level (idx_lev and prof_wgt) in nwp profiles
            temp_prof_nwp  => t_prof_nwp(:,x_nwp,y_nwp)    ! - temperature profile
            hgt_prof_nwp   => z_prof_nwp(:,x_nwp,y_nwp)     ! - height profile
            
              ! - level indicies and weights in nwp and rtm profiles
            placeholder_cld = cld_height_loc
            call t_to_pz_from_profile ( cld_temp_loc , &
                                temp_prof_nwp , &
                                p_std_nwp , &
                                hgt_prof_nwp, &
                                sfc_level_nwp(x_nwp,y_nwp), &
                                tropo_level_nwp(x_nwp,y_nwp), &
                                inversion_level_nwp(x_nwp,y_nwp), &
                                placeholder_cld, cld_height_loc, idx_lev_nwp , prof_wgt_nwp)
 
             call t_to_pz_from_profile ( cld_temp_loc , &
                               rtm(x_nwp,y_nwp)%t_prof , &
                               p_std_rtm , &
                               rtm(x_nwp,y_nwp)%z_prof, &
                               rtm(x_nwp,y_nwp)%sfc_level, &
                               rtm(x_nwp,y_nwp)%tropo_level, &
                               rtm(x_nwp,y_nwp)%inversion_level, &
                               placeholder_cld, cld_height_loc, idx_lev_rtm , prof_wgt_rtm)
                             
            temp_prof_nwp  => null()
            hgt_prof_nwp   => null()
            
            wv_prof_nwp    => tpw_prof_nwp(:,x_nwp,y_nwp)   ! - total water path profile     
            
            if (idx_lev_nwp == size(wv_prof_nwp)) then
               dcomp_rtm % tpw_ac(elem_idx,line_idx) = wv_prof_nwp ( idx_lev_nwp )  ! - catch invalid data
            else
               dcomp_rtm % tpw_ac(elem_idx,line_idx) = wv_prof_nwp(idx_lev_nwp)  + prof_wgt_nwp*(wv_prof_nwp(idx_lev_nwp+1) &
                                                  - wv_prof_nwp(idx_lev_nwp))
            end if
            wv_prof_nwp    => null()
            
            ! - viewing angle index
            ivza = zen_idx_rtm(elem_idx,line_idx) 
            if ( ivza .lt. 1) cycle
           
            
            dcomp_rtm % trans_ir_ac(elem_idx,line_idx) =  &
                                  & rtm(x_nwp,y_nwp) % d(ivza) % ch(20) % trans_atm_profile(idx_lev_rtm) +    &
                                  & prof_wgt_rtm * ( rtm(x_nwp,y_nwp) % d(ivza) % ch(20) % trans_atm_profile(idx_lev_rtm+1) -  &
                                  & rtm(x_nwp,y_nwp) % d(ivza) % ch(20) % trans_atm_profile(idx_lev_rtm))
                  
            dcomp_rtm % trans_ir_ac_nadir(elem_idx,line_idx) = &
                                 & dcomp_rtm % trans_ir_ac(elem_idx,line_idx) ** (cos ( geo % satzen(elem_idx, line_idx) * PI / 180.  ))
            
            dcomp_rtm % sfc_nwp  (elem_idx,line_idx)  = psfc_nwp(x_nwp,y_nwp)
            
            clear_trans_prof_rtm = rtm(x_nwp,y_nwp) % d(ivza) % ch(20) % trans_atm_profile
            clear_rad_prof_rtm   = rtm(x_nwp,y_nwp) % d(ivza) % ch(20) % rad_atm_profile

            dcomp_rtm % rad_clear_sky_toc_ch20 (elem_idx,line_idx) = clear_rad_prof_rtm (idx_lev_rtm)
        !-->dcomp_rtm % rad_clear_sky_toa_ch20 (elem_idx,line_idx) = ch(20)%rad_toa_clear(elem_idx,line_idx)             
            dcomp_rtm % rad_clear_sky_toa_ch20 (elem_idx,line_idx) = Rad_Clear_Ch20_Solar_Rtm(Elem_Idx,Line_Idx)
         
         end do element_loop
      end do line_loop
         
   end subroutine perform_rtm_dcomp
   
   
   subroutine allocate_dcomp_rtm ( dcomp_rtm, x_dim, y_dim)
      type ( dcomp_rtm_type) , intent(inout) :: dcomp_rtm
      integer , intent (in) :: x_dim, y_dim
      
      
      allocate ( dcomp_rtm % trans_ir_ac (x_dim , y_dim))
      allocate ( dcomp_rtm % trans_ir_ac_nadir(x_dim , y_dim))
      allocate ( dcomp_rtm % tpw_ac (x_dim , y_dim))
      allocate ( dcomp_rtm % sfc_nwp (x_dim , y_dim))
      allocate ( dcomp_rtm % rad_clear_sky_toc_ch20 (x_dim , y_dim))
      allocate ( dcomp_rtm % rad_clear_sky_toa_ch20 (x_dim , y_dim))
      allocate ( dcomp_rtm % ozone_path (x_dim , y_dim))
      
   end subroutine allocate_dcomp_rtm
   
   
   subroutine deallocate_it ( self )
      class ( dcomp_rtm_type)  :: self
      
      deallocate ( self % trans_ir_ac )
      deallocate ( self % trans_ir_ac_nadir)
      deallocate ( self % tpw_ac )
      deallocate ( self % sfc_nwp )
      deallocate ( self % rad_clear_sky_toc_ch20 )
      deallocate ( self % rad_clear_sky_toa_ch20 )
      deallocate ( self % ozone_path )
   
   end subroutine deallocate_it
   
   
   
   ! -----------------------------------------------------------------
   ! computes pressure and height from temperature from given profiles
   ! considers possible inversion tropopause 
   !
   !  
   subroutine t_to_pz_from_profile ( temp , &
                        t_prof , p_prof , z_prof , &
                        sfc_idx, trp_idx, inv_idx, &
                        press, height, lev_idx , lev_wgt) 
  
      implicit none
      
      real, intent(in) :: temp
      real, dimension(:), intent(in) :: t_prof
      real, dimension(:), intent(in) :: p_prof
      real, dimension(:), intent(in) :: z_prof
      integer (kind = int1) , intent(in) :: trp_idx
      integer (kind = int1), intent(in) :: sfc_idx
      integer (kind = int1), intent(in) :: inv_idx
  
      real, intent(out) :: press
      real, intent(out) :: height
      integer, intent(out) :: lev_idx
      real, intent(out) :: lev_wgt
  
      logical :: inv_flag 
      integer :: n_prof
      integer :: n_levs_temp
      real :: d_press
      real :: d_height
      real :: d_temp
 
      n_prof = size(p_prof)
 
      !-- temperature warmer than surface?
      if (temp > maxval (t_prof(trp_idx : sfc_idx) ) ) then
         press = p_prof(sfc_idx)
         height = z_prof(sfc_idx)
         lev_idx = sfc_idx
         lev_wgt = 0.
         return
      end if 

      if (temp .lt. minval (t_prof(trp_idx : sfc_idx) ) ) then
         press = p_prof(trp_idx)
         height = z_prof(trp_idx)
         lev_idx = trp_idx
         lev_wgt = 0.
         return
      end if 
 
      inv_flag = .false.
      if (inv_idx .gt. 0) then
         n_levs_temp = sfc_idx - inv_idx + 1
         call dcomp_locate ( t_prof (inv_idx : sfc_idx) , n_levs_temp , temp, lev_idx  )
         if ( ( lev_idx .gt. 0 ) .and. ( lev_idx .lt. n_levs_temp -1 ) ) then
            lev_idx = lev_idx + inv_idx - 1
            inv_flag = .true.
         end if 
      end if
 
      !--- if no solution within an inversion, look above
      if ( inv_flag .eqv. .false.) then
         n_levs_temp = sfc_idx - trp_idx + 1
         call dcomp_locate ( t_prof (trp_idx : sfc_idx) , n_levs_temp , temp, lev_idx  )
         lev_idx = lev_idx + trp_idx - 1
         lev_idx = max ( 1, min ( n_prof -1 , lev_idx  ))
      end if

      !-- if solution is above trop, set to trop values
      if (lev_idx .lt. trp_idx) then
         press = p_prof(trp_idx)
         height = z_prof(trp_idx)
         lev_idx = trp_idx
         lev_wgt = 0
         return
      end if

      !--- determine derivatives
      d_press = p_prof(lev_idx+1) - p_prof(lev_idx)
      d_height = z_prof(lev_idx+1) - z_prof(lev_idx)
      d_temp = t_prof(lev_idx+1) - t_prof(lev_idx)
   
      if ( d_temp .ne. 0.0 ) then
         lev_wgt = 1./d_temp * (temp - t_prof(lev_idx))
         press = p_prof (lev_idx) + d_press * lev_wgt
         height  = z_prof (lev_idx) + d_height * lev_wgt 
      else
         press = p_prof(lev_idx)
         height = z_prof(lev_idx)
         lev_wgt = 0
      end if
   end subroutine t_to_pz_from_profile

   !--------------------------------------------------------------------------
   subroutine dcomp_locate(xx, n, x, j)
      integer,                        intent(in)  :: n
      integer,                        intent(out) :: j
      real , intent(in)  :: x
      real , dimension(:), intent(in)  :: xx
      integer :: i, jl, jm, ju

      jl = 0
      ju = n + 1
      do i = 1, 2 * n
         if (ju-jl <= 1) exit
      
         jm = (ju + jl) / 2
         if ((xx(n) >= xx(1)) .eqv. (x >= xx(jm))) then
            jl = jm
         else
            ju = jm
         end if
      end do
      if (x == xx(1)) then
         j =1
      else if (x == xx(n)) then
         j = n - 1
      else
         j = jl
      end if

   end subroutine dcomp_locate



end module dcomp_rtm_module
