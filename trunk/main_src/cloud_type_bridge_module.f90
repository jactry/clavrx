! $Header$
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
!             num_pix                                    integer
!             num_scans_read                             integer
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
!              ET_cloudiness_class                type with enumerators for cloud mask
!
!--------------------------------------------------------------------------------------

module cloud_type_bridge_module

   
   use PIXEL_COMMON, only : &
        solzen  &
      , satzen &   
      , zen_idx_rtm &
      , Chan_On_Flag_Default &
      , ch &
      , num_pix &
      , num_scans_read &
      , i_nwp &
      , j_nwp &
      , Bt_Ch31_Max_3x3 &
      , Bt_Ch31_Std_3x3 &
      , Covar_Ch27_Ch31_5x5 &
      , cld_type &
      , cld_phase &
      , cld_mask &
      , i_lrc, j_lrc &
      , Beta_11um_12um_Tropo_Rtm &
      , Beta_11um_133um_Tropo_Rtm
                 
  
   use CLOUD_TYPE_ALGO_MODULE, only : &
       cloud_type_pixel &
       , cloud_type_input_type &
       , ET_cloud_type &
       , set_cloud_phase
       
   use NAIVE_BAYESIAN_CLOUD_MASK_MODULE, only: &
      ET_cloudiness_class
   
   use RTM_COMMON , only: &
      p_std_rtm &
      , rtm
   
   implicit none
   
   public :: CLOUD_TYPE_BRIDGE  

contains
   subroutine cloud_type_bridge
      implicit none
      
      type ( cloud_type_input_type) :: type_inp
      
      integer :: ctype
      integer :: i , j   
      integer :: ii , jj
      real :: ice_prob 
      integer :: cld_type_lrc
      
      ! ------  Executable  ------------------------------------
      
      allocate ( type_inp % rtm % p_prof ,source = p_std_rtm ) 
      
      
      ice_prob = -999.
      
      type_inp % sat % chan_on = Chan_On_Flag_Default == 1
      
      
      ! -----------    loop over LRC core pixels to get ice probabbilty -----         
      elem_loop: do  j = 1,num_scans_read
         line_loop: do i = 1, num_pix  
         
            if (    cld_mask ( i,j) == ET_cloudiness_class % CLEAR ) then
               cld_type (i , j ) = ET_cloud_type % CLEAR
                cycle
            end if
                
            if (    cld_mask ( i,j) == ET_cloudiness_class % PROB_CLEAR ) then
               cld_type (i , j ) = ET_cloud_type % PROB_CLEAR
               cycle
            end if
            
            ! - take only LRC cores
            if ( i /= i_lrc (i,j) .or. j /= j_lrc (i,j) ) cycle
                             
            call POPULATE_INPUT ( i, j , type_inp )
            call CLOUD_TYPE_PIXEL  ( type_inp, ctype , ice_prob_out = ice_prob )
            cld_type (i,j)  = ctype
            if ( ctype < 0 ) print*,'lrc',i,j
            
            deallocate ( type_inp % rtm % rad_ch31_bb_prof )
            deallocate ( type_inp % rtm % t_prof )
            deallocate ( type_inp % rtm % z_prof )
                    
         
         end do   line_loop
      end do elem_loop
      
      
      ! - now loop over all non lrc-cores
      elem_loop1: do  j = 1,num_scans_read
         line_loop1: do i = 1, num_pix  
            
            if (    cld_mask ( i,j) == ET_cloudiness_class % CLEAR ) then
               cld_type (i , j ) = ET_cloud_type % CLEAR
                cycle
            end if
                
            if (    cld_mask ( i,j) == ET_cloudiness_class % PROB_CLEAR ) then
               cld_type (i , j ) = ET_cloud_type % PROB_CLEAR
               cycle
            end if   
               
            
            ii = i_lrc (i,j)
            jj = j_lrc (i,j)
                        
            !- we dont need the lrc cores again
            if ( i == ii .and. j == jj ) cycle
            
            cld_type_lrc = cld_type (ii , jj )
                             
            call POPULATE_INPUT ( i, j , type_inp )
            call CLOUD_TYPE_PIXEL  ( type_inp, ctype , ice_prob_out = ice_prob )
            if ( ctype < 0 ) print*,i,j
            
            ! - compare this ctype with LRC
           
            !  - identical or lrc is not valid => take the current
            if ( ctype == cld_type_lrc .or. ii < 1 .or. jj < 1 .or. cld_type_lrc <= 0 ) then
               cld_type (i,j)  = ctype
            else
               ! - if LRC core is water phase ==> use lrc
               if ( cld_type (ii,jj) >= ET_cloud_type % FIRST_WATER &
                  .and. cld_type (ii,jj) <= ET_cloud_type % LAST_WATER ) then
               
                  cld_type (i , j ) = cld_type_lrc
                  
               ! - LRC core is ice phase
               else if ( (ctype  == ET_cloud_type % FOG &
                  .or. ctype == ET_cloud_type % WATER) &
                  .and. ( cld_type_lrc == ET_cloud_type % CIRRUS &
                  .or. cld_type_lrc == ET_cloud_type % OVERLAP &
                  .or. cld_type_lrc == ET_cloud_type % OPAQUE_ICE)) then
                  
                     cld_type (i , j ) = ET_cloud_type % CIRRUS
               
               ! - LRC core is ice phase and current is supercooled => switch to ice
               else if ( ( cld_type_lrc == ET_cloud_type % CIRRUS & 
                         .or. cld_type_lrc == ET_cloud_type % OPAQUE_ICE ) &
                        .and. ctype ==  ET_cloud_type % SUPERCOOLED ) then
                  
                  call CLOUD_TYPE_PIXEL  ( type_inp, ctype , force_ice = .true. )
                     cld_type (i,j)  = ctype
                    
               ! -- this is mainly cirrus / opaque ice => keep current
               else 
                  cld_type (i,j)  = ctype               
                                  
               end if      
               
            end if
            
            deallocate ( type_inp % rtm % rad_ch31_bb_prof )
            deallocate ( type_inp % rtm % t_prof )
            deallocate ( type_inp % rtm % z_prof )
                    
         
         end do   line_loop1
      end do elem_loop1      
      
      
      
      deallocate ( type_inp % rtm % p_prof )
      
      call set_cloud_phase ( cld_type, cld_phase) 
           
      
   end subroutine cloud_type_bridge
   
   ! --------- --------------- ---
   !
   ! --------- ---------------
   subroutine populate_input ( i, j , type_inp)
      integer, intent(in) :: i
      integer, intent(in) :: j
      type ( cloud_type_input_type) :: type_inp
      
      integer :: nwp_lon_idx , nwp_lat_idx
      integer :: vza_idx
      
      Nwp_Lon_Idx = I_Nwp( i , j )
      Nwp_Lat_Idx = J_Nwp( i , j )
      Vza_Idx = zen_Idx_Rtm( i , j )            
      
      ! - sat
      type_inp % sat % rad_ch31 = ch(31) % rad_toa ( i,j )
      type_inp % sat % bt_ch31 =  ch(31) % bt_toa  ( i,j )
      type_inp % sat % bt_ch32 =  ch(32) % bt_toa  ( i,j )
      type_inp % sat % ref_ch6 =  ch(6)  % ref_toa  ( i,j )
      type_inp % sat % ref_ch20 = ch(20) % ref_toa ( i,j )
      
      ! - rtm
      allocate ( type_inp % rtm % rad_ch31_bb_prof &
         , source = Rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%d(Vza_Idx)%ch(31)%Rad_BB_Cloud_Profile)
      
      allocate ( type_inp % rtm % t_prof , source = rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%T_prof )
      allocate ( type_inp % rtm % z_prof , source = rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%z_prof )
      
      type_inp % rtm % tropo_lev = rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%Tropo_Level
      type_inp % rtm % sfc_lev = rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%sfc_Level
      
      type_inp % rtm % bt_ch31_3x3_max     = Bt_Ch31_Max_3x3( i,j )
      type_inp % rtm % bt_ch31_3x3_std     = Bt_Ch31_Std_3x3( i,j )
      type_inp % rtm % Beta_11um_12um_Tropo = Beta_11um_12um_Tropo_Rtm( i,j )
      type_inp % rtm % Beta_11um_133um_Tropo = Beta_11um_133um_Tropo_Rtm( i,j )
      
      type_inp % rtm % Covar_Ch27_Ch31_5x5 = -999.
      
      if ( chan_on_flag_default(27) == 1 ) then
         type_inp % rtm % Covar_Ch27_Ch31_5x5 = Covar_Ch27_Ch31_5x5( i,j )
         type_inp % sat % rad_ch27 = ch(27) % rad_toa ( i,j )
         type_inp % sat % bt_ch27 =  ch(27) % bt_toa  ( i,j )
      end if   
      
      type_inp % rtm % ref_ch6_clear       = ch(6)%Ref_Toa_Clear( i,j )
      type_inp % rtm % bt_ch31_atm_sfc     = ch(31)%Bt_Toa_Clear( i,j )
      type_inp % rtm % bt_ch32_atm_sfc     = ch(32)%Bt_Toa_Clear( i,j )
      type_inp % rtm % emiss_tropo_ch31    = ch(31)%Emiss_Tropo( i,j )
      
      ! - geo
      type_inp % geo % sol_zen     = solzen ( i , j )
       type_inp % geo % sat_zen     = satzen ( i , j )
      
      !- sfc
      type_inp % sfc % emiss_ch20     = ch(20) % sfc_emiss ( i , j )
            
   
   
   end subroutine populate_input

end module cloud_type_bridge_module
