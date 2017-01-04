! $Header: https://svn.ssec.wisc.edu/repos/cloud_team_clavrx/trunk/main_src/cloud_type_bridge_module.f90 412 2014-06-11 17:35:49Z heidinger $
!
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 
!
! NAME: cloud_type_bridge_module.f90 (src)
!       cloud_type_bridge (program)
!
! PURPOSE: 
!       Builds the interface to Cloud Type algorithm
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
!  HISTORY:
!      2014/04/30:    new interface (AW)
!          /05/02:    add cld_phase (AW)
!
!  GLOBAL VARIABLES:
!
!    FROM PIXEL_COMMON:
!      1. work as input
!        1.1 configuration
!             Chan_On_Flag_Default                       integer (42)
!             Image%Number_Of_Elements                                    integer
!             Image%Number_Of_Lines_Read_This_Segment                             integer
!        1.2 geo data:                   
!              solzen                                    real (:,:)
!              satzen                                    real (:,:)
!        1.3 surface 
!        1.4 rtm / statistics 
!
!              bt_ch31_max_3x3                           real (:,:)
!              bt_Ch31_Std_3x3                           real (:,:)
!              Covar_Ch27_Ch31_5x5                       real (:,:)
!              Beta_11um_12um_Tropo_Rtm                  real (:,:)
!              Beta_11um_133um_Tropo_Rtm                  real (:,:)
!              zen_idx_rtm                               integer (:,:)
!              i_nwp
!              j_nwp
!              i_lrc
!              j_lrc
!        1.5 observations
!              ch                                        type ( observations )  
!        1.6 products
!             cld_mask    
!       
!      2. work as output  
!              cld_type                                  integer (:,:) 
!           
!     FROM RTM_COMMON:
!            rtm                                         type ( ) 
!            p_std_rtm                                   real(:)
!
!     CLOUD_TYPE_ALGO_MODULE
!             cloud_type_pixel                     subroutine
!                  subroutine which retrieves cloud type
!             cloud_type_input_type                      type definition
!                   structure definition of cloud_type retrieval input  
!
!     NAIVE_BAYESIAN_CLOUD_MASK_MODULE
!              et_cloudiness_class                type with enumerators for cloud mask
!
!--------------------------------------------------------------------------------------

module cloud_type_bridge_module

   
   use PIXEL_COMMON, only : &
        image &
      , sensor &
      , geo &
      , zen_idx_rtm &
      , ch &
      , i_nwp &
      , j_nwp &
      , Bt_Ch27_Max_3x3 &
      , Bt_Ch31_Max_3x3 &
      , Bt_Ch31_Std_3x3 &
      , Covar_Ch27_Ch31_5x5 &
      , cld_type &
      , cld_phase &
      , CLDMASK &
      , Dust_Mask &
      , Smoke_Mask &
      , Fire_Mask &
      , i_lrc, j_lrc &
      , Beta_11um_12um_Tropo_Rtm &
      , Beta_11um_133um_Tropo_Rtm &
      , Diag_Pix_Array_1 &
      , Diag_Pix_Array_2 &
      , Diag_Pix_Array_3 &
      , bad_pixel_mask
                 
   use CONSTANTS, only : &
        Cloud_Type_Version

   use CLOUD_TYPE_ALGO_MODULE, only : &
       cloud_type_pixel &
       , cloud_type_input_type &
       , cloud_type_diag_type &
       , et_cloud_type &
       , et_cloudiness_class &
       , set_cloud_phase
       
   use RTM_COMMON , only: &
      p_std_rtm &
      , rtm

   use CLAVRX_MESSAGE_MODULE, only: MESG
   
   implicit none
   
   public :: CLOUD_TYPE_BRIDGE  
   public :: SET_CLOUD_TYPE_VERSION

contains

   !====================================================================
   !  record svn version as a global variable for output to hdf
   !====================================================================
   subroutine SET_CLOUD_TYPE_VERSION()
      Cloud_Type_Version = "$Id: cloud_type_bridge_module.f90 412 2014-06-11 17:35:49Z heidinger $"
   end subroutine Set_CLOUD_TYPE_VERSION


   !====================================================================
   ! universal cloud type bridge
   !====================================================================
   subroutine cloud_type_bridge
      implicit none
      
      type ( cloud_type_input_type) :: type_inp
      type ( cloud_type_diag_type) :: diag_out
      
      integer :: ctype
      integer :: i , j   
      integer :: ii , jj
      real :: ice_prob 
      integer :: cld_type_lrc
      logical, save:: First_Call = .true.
      
      
      ! ------  Executable  ------------------------------------
      if (First_Call .eqv. .true.) then
        call MESG('Cloud Type starts ', color = 46)
      endif

      ice_prob = -999.0
      
      type_inp % sat % chan_on = Sensor%Chan_On_Flag_Default == 1
      
      !-----------    loop over LRC core pixels to get ice probabbilty -----         
      elem_loop: do  j = 1,Image%Number_Of_Lines_Read_This_Segment
         line_loop: do i = 1, Image%Number_Of_Elements  
            
            if ( bad_pixel_mask (i,j) == 1 ) then
               cld_type (i,j ) = et_cloud_type % MISSING
               cycle
            end if 
            
            if (CLDMASK%Cld_Mask (i,j) == et_cloudiness_class % CLEAR ) then
               cld_type (i,j ) = et_cloud_type % CLEAR
                cycle
            end if
                
            if (CLDMASK%Cld_Mask (i,j) == et_cloudiness_class % PROB_CLEAR ) then
               cld_type (i,j ) = et_cloud_type % PROB_CLEAR
               cycle
            end if

            if (Dust_Mask (i,j) == 1 ) then
               cld_type (i,j ) = et_cloud_type % DUST
               cycle
            end if

            if (Smoke_Mask (i,j) == 1 ) then
               cld_type (i,j ) = et_cloud_type % SMOKE
               cycle
            end if

            if (Fire_Mask (i,j) == 1 ) then
               cld_type (i,j ) = et_cloud_type % FIRE
               cycle
            end if
            
            ! - take only LRC cores
            if ( i /= i_lrc (i,j) .or. j /= j_lrc (i,j) ) cycle
                             
            call POPULATE_INPUT ( i, j , type_inp )
            call CLOUD_TYPE_PIXEL  ( type_inp, ctype , diag_out, ice_prob_out = ice_prob )
            cld_type (i,j)  = ctype

            call DEALLOCATE_INP ( type_inp )
         
         end do   line_loop
      end do elem_loop

 
      ! - now loop over all non lrc-cores
      elem_loop1: do  j = 1,Image%Number_Of_Lines_Read_This_Segment
         line_loop1: do i = 1, Image%Number_Of_Elements  
            
            if ( bad_pixel_mask (i,j) == 1 ) then
               cld_type (i,j ) = et_cloud_type % MISSING
               cycle
            end if 
            
            if (CLDMASK%Cld_Mask ( i,j) == et_cloudiness_class % CLEAR ) then
               cld_type (i , j ) = et_cloud_type % CLEAR
                cycle
            end if
                
            if (CLDMASK%Cld_Mask ( i,j) == et_cloudiness_class % PROB_CLEAR ) then
               cld_type (i , j ) = et_cloud_type % PROB_CLEAR
               cycle
            end if   

            if (Dust_Mask (i,j) == 1 ) then
               cld_type (i,j ) = et_cloud_type % DUST
               cycle
            end if

            if (Smoke_Mask (i,j) == 1 ) then
               cld_type (i,j ) = et_cloud_type % SMOKE
               cycle
            end if

            if (Fire_Mask (i,j) == 1 ) then
               cld_type (i,j ) = et_cloud_type % FIRE
               cycle
            end if
            
            ii = i_lrc (i,j)
            jj = j_lrc (i,j)
                        
            !- we dont need the lrc cores again
            if ( i == ii .and. j == jj ) cycle

            call POPULATE_INPUT ( i, j , type_inp )
            call CLOUD_TYPE_PIXEL  ( type_inp, ctype , diag_out, ice_prob_out = ice_prob )
            
            !--- set lrc value
            cld_type_lrc = et_cloud_type % UNKNOWN
            if ( ii > 0 .and. jj > 0) then
               cld_type_lrc = cld_type (ii , jj )
            endif

            ! - compare this ctype with LRC
           
            !  - identical or lrc is not valid => take the current
            if ( ctype == cld_type_lrc .or. ii < 1 .or. jj < 1 .or. cld_type_lrc == et_cloud_type%UNKNOWN) then
               cld_type (i,j)  = ctype
            else

               ! - if LRC core is water phase ==> use lrc
               if ( cld_type_lrc >= et_cloud_type % FIRST_WATER &
                  .and. cld_type_lrc <= et_cloud_type % LAST_WATER ) then
               
                   ! AKH says not to overwrite WATER,  FOG or SUPERCOOLED
                   if ( ctype >= et_cloud_type % FIRST_WATER &
                        .and. ctype <= et_cloud_type % LAST_WATER ) then
                        
                      cld_type(i,j) = ctype
                   else
                     ! - the original ice pixels should also be check on supercool, fog or water.
                     call CLOUD_TYPE_PIXEL  ( type_inp, ctype , diag_out, force_water = .true. )
                     cld_type (i,j)  = ctype
                   end if
                  
               ! - LRC core is ice phase
               else if ( (ctype  == et_cloud_type % FOG &
                  .or. ctype == et_cloud_type % WATER) &
                  .and. ( cld_type_lrc == et_cloud_type % CIRRUS &
                  .or. cld_type_lrc == et_cloud_type % OVERLAP &
                  .or. cld_type_lrc == et_cloud_type % OPAQUE_ICE)) then
                  
                     cld_type (i , j ) = et_cloud_type % CIRRUS
               
               ! - LRC core is ice phase and current is supercooled => switch to ice
               else if ( ( cld_type_lrc == et_cloud_type % CIRRUS & 
                         .or. cld_type_lrc == et_cloud_type % OPAQUE_ICE &
                         .or. cld_type_lrc == et_cloud_type % OVERLAP ) &
                        .and. ctype ==  et_cloud_type % SUPERCOOLED ) then
                  
                  call CLOUD_TYPE_PIXEL  ( type_inp, ctype , diag_out, force_ice = .true. )
                     cld_type (i,j)  = ctype
                    
               ! -- this is mainly cirrus / opaque ice => keep current
               else 
                  cld_type (i,j)  = ctype                                                 
               end if                     
            end if

            call DEALLOCATE_INP ( type_inp )
         
         end do   line_loop1
      end do elem_loop1      
      
      call set_cloud_phase ( cld_type, cld_phase) 


      First_Call = .false.
      
   end subroutine cloud_type_bridge
   
   ! --------- --------------- ---
   !
   !
   ! --------- ---------------
   subroutine POPULATE_INPUT ( i, j , type_inp)
      integer, intent(in) :: i
      integer, intent(in) :: j
      type ( cloud_type_input_type) :: type_inp
      
      integer :: nwp_lon_idx , nwp_lat_idx
      integer :: vza_idx
      
      integer :: n_rtm_prof
      
      
      Nwp_Lon_Idx = I_Nwp( i , j )
      Nwp_Lat_Idx = J_Nwp( i , j )
      Vza_Idx = zen_Idx_Rtm( i , j )            
      
      !-----------------------------------------------------------------------------------
      ! - sat
      !-----------------------------------------------------------------------------------
      if (sensor % chan_on_flag_default(31) == 1 ) type_inp % sat % rad_ch31 = ch(31) % rad_toa ( i,j )
      if (sensor % chan_on_flag_default(31) == 1 ) type_inp % sat % bt_ch31 =  ch(31) % bt_toa  ( i,j )
      if (sensor % chan_on_flag_default(32) == 1 ) type_inp % sat % bt_ch32 =  ch(32) % bt_toa  ( i,j )
      if (sensor % chan_on_flag_default(6) == 1 ) type_inp % sat % ref_ch6 =  ch(6)  % ref_toa  ( i,j )
      if (sensor % chan_on_flag_default(20) == 1 ) type_inp % sat % ref_ch20 = ch(20) % ref_toa ( i,j )

      if (sensor % chan_on_flag_default(27) == 1 ) then
         type_inp % sat % rad_ch27 = ch(27) % rad_toa (i,j)
         type_inp % sat % bt_ch27 =  ch(27) % bt_toa  (i,j)
      end if  
      
      if (sensor % chan_on_flag_default(29) == 1 ) then
         type_inp % sat % rad_ch29 = ch(29) % rad_toa (i,j)
         type_inp % sat % bt_ch29 =  ch(29) % bt_toa  (i,j)
      end if         
      
      !-----------------------------------------------------------------------------------
      ! - rtm
      !-----------------------------------------------------------------------------------
      n_rtm_prof = size ( rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%T_prof )
      allocate ( type_inp % rtm % t_prof(n_rtm_prof ) , source = rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%T_prof )
      allocate ( type_inp % rtm % z_prof(n_rtm_prof ) , source = rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%z_prof )
      type_inp % rtm % tropo_lev = rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%Tropo_Level
      type_inp % rtm % sfc_lev = rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%sfc_Level
      
      
      
      if (sensor % chan_on_flag_default(6) == 1 )  then
         type_inp % rtm % ref_ch6_clear       = ch(6)%Ref_Toa_Clear( i,j )
      endif
      if (sensor % chan_on_flag_default(31) == 1 )  then
         
         
         allocate ( type_inp % rtm % rad_ch31_bb_prof (n_rtm_prof ) &
              , source = Rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%d(Vza_Idx)%ch(31)%Rad_BB_Cloud_Profile)
          
         type_inp % rtm % bt_ch31_3x3_max    = Bt_Ch31_Max_3x3( i,j )
         type_inp % rtm % bt_ch31_3x3_std    = Bt_Ch31_Std_3x3( i,j )
         type_inp % rtm % rad_ch31_atm_sfc   = ch(31)%Rad_Toa_Clear(i,j)
         type_inp % rtm % bt_ch31_atm_sfc    = ch(31)%Bt_Toa_Clear( i,j )
         type_inp % rtm % emiss_tropo_ch31   = ch(31)%Emiss_Tropo( i,j )
         if (sensor % chan_on_flag_default(27) == 1 )  then
            type_inp % rtm % Covar_Ch27_Ch31_5x5 = Covar_Ch27_Ch31_5x5( i,j )
         endif
         if (sensor % chan_on_flag_default(32) ==1 )  then
            type_inp % rtm % Beta_11um_12um_Tropo  = Beta_11um_12um_Tropo_Rtm( i,j )
            type_inp % rtm % bt_ch32_atm_sfc       = ch(32)%Bt_Toa_Clear( i,j )
         endif
         if (sensor % chan_on_flag_default(33) ==1 )  then
            type_inp % rtm % Beta_11um_133um_Tropo = Beta_11um_133um_Tropo_Rtm( i,j )
         endif
      endif
      
      if (sensor % chan_on_flag_default(27) == 1) then
         type_inp % rtm % bt_ch27_3x3_max    = Bt_Ch27_Max_3x3( i,j )
         type_inp % rtm % rad_ch27_atm_sfc = ch(27)%Rad_Toa_Clear(i,j)
         allocate ( type_inp % rtm % rad_ch27_bb_prof(n_rtm_prof ) &
         , source = Rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%d(Vza_Idx)%ch(27)%Rad_BB_Cloud_Profile)
      end if   
      
      !-----------------------------------------------------------------------------------
      ! - geo
      !-----------------------------------------------------------------------------------
      type_inp % geo % sol_zen     = geo % solzen ( i , j )
      type_inp % geo % sat_zen     = geo % satzen ( i , j )
      
      !-----------------------------------------------------------------------------------
      !- sfc
      !-----------------------------------------------------------------------------------
      type_inp % sfc % emiss_ch20     = ch(20) % sfc_emiss ( i , j )
            
   end subroutine POPULATE_INPUT
   
   !----------------------------------------------------------------------------------------- 
   !
   !----------------------------------------------------------------------------------------- 
   subroutine DEALLOCATE_INP ( type_inp)
      type ( cloud_type_input_type) :: type_inp
       
       deallocate ( type_inp % rtm % rad_ch31_bb_prof )
         if ( allocated( type_inp % rtm % rad_ch27_bb_prof)) deallocate ( type_inp % rtm % rad_ch27_bb_prof )
       deallocate ( type_inp % rtm % t_prof )
       deallocate ( type_inp % rtm % z_prof )
   
   end subroutine DEALLOCATE_INP

end module cloud_type_bridge_module
