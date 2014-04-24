! $Header:$


module cloud_type_bridge_module

   
   use PIXEL_COMMON, only : &
      solzen  &
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
      , cld_phase &
      , cld_type &
      , cld_mask &
      , i_lrc, j_lrc
                 
   use CLOUD_PHASE_ALGO_MODULE
   use CLOUD_TYPE_ALGO_MODULE 
   use NAIVE_BAYESIAN_CLOUD_MASK_MODULE
   
   use RTM_COMMON , only: &
      p_std_rtm &
      , rtm
   
   implicit none
   
   public :: CLOUD_TYPE_BRIDGE  

contains
   subroutine cloud_type_bridge
      implicit none
      
      type ( cloud_type_input_type) :: type_inp
      type ( cloud_phase_input_type) :: phase_inp
      integer :: phase
      integer :: ctype
      integer :: i , j
      integer :: nwp_lon_idx , nwp_lat_idx
      integer :: vza_idx
      integer :: ii , jj
      
      allocate ( phase_inp % rtm % p_prof ,source = p_std_rtm ) 
      phase_inp % sat % chan_on = Chan_On_Flag_Default == 1
      Cld_Phase  = ET_cloud_phase % CLEAR
      
      ! -----------    loop over pixels -----   
      line_loop: do i = 1, num_pix
         elem_loop: do  j = 1,num_scans_read
           
            
            if (    cld_mask ( i,j) == ET_cloudiness_class % CLEAR &
               .or. cld_mask ( i,j) == ET_cloudiness_class % PROB_CLEAR) cycle   
               
            Nwp_Lon_Idx = I_Nwp( i , j )
            Nwp_Lat_Idx = J_Nwp( i , j )
            Vza_Idx = zen_Idx_Rtm( i , j )            
            
            ! - sat
            phase_inp % sat % rad_ch31 = ch(31) % rad_toa ( i,j )
            phase_inp % sat % bt_ch31 =  ch(31) % bt_toa  ( i,j )
            phase_inp % sat % bt_ch32 =  ch(32) % bt_toa  ( i,j )
            phase_inp % sat % ref_ch6 =  ch(6)  % ref_toa  ( i,j )
            phase_inp % sat % ref_ch20 = ch(20) % ref_toa ( i,j )
            
            ! - rtm
            allocate ( phase_inp % rtm % rad_ch31_bb_prof &
               , source = Rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%d(Vza_Idx)%ch(31)%Rad_BB_Cloud_Profile)
            
            allocate ( phase_inp % rtm % t_prof , source = rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%T_prof )
            allocate ( phase_inp % rtm % z_prof , source = rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%z_prof )
            
            phase_inp % rtm % tropo_lev = rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%Tropo_Level
            phase_inp % rtm % sfc_lev = rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%sfc_Level
            
            phase_inp % rtm % bt_ch31_3x3_max     = Bt_Ch31_Max_3x3( i,j )
            phase_inp % rtm % bt_ch31_3x3_std     = Bt_Ch31_Std_3x3( i,j )
            if ( chan_on_flag_default(27) == 1 ) phase_inp % rtm % Covar_Ch27_Ch31_5x5 = Covar_Ch27_Ch31_5x5( i,j )
            phase_inp % rtm % ref_ch6_clear       = ch(6)%Ref_Toa_Clear( i,j )
            phase_inp % rtm % bt_ch31_atm_sfc     = ch(31)%Bt_Toa_Clear( i,j )
            phase_inp % rtm % bt_ch32_atm_sfc     = ch(32)%Bt_Toa_Clear( i,j )
            phase_inp % rtm % emiss_tropo_ch31    = ch(31)%Emiss_Tropo( i,j )
            
            ! - geo
            phase_inp % geo % sol_zen     = solzen ( i , j )
            
            !- sfc
            phase_inp % sfc % emiss_ch20     = ch(20) % sfc_emiss ( i , j )
            
            call cloud_phase_algo ( phase_inp, phase , ctype)
            
            Cld_Phase ( i, j) = phase
            Cld_Type ( i ,j ) = ctype
            
            deallocate ( phase_inp % rtm % rad_ch31_bb_prof )
            deallocate ( phase_inp % rtm % t_prof )
            deallocate ( phase_inp % rtm % z_prof )
            
            
            !call cloud_type_algo
         end do elem_loop
      end do   line_loop
      
      
      ! - redo for lrc beighborohoods
      line_loop_lrc: do i = 1, num_pix
         elem_loop_lrc: do  j = 1,num_scans_read
            if ( cld_phase ( i,j) == ET_cloud_phase % CLEAR ) cycle
            
            ii = i_lrc(i,j)
            jj = j_lrc(i,j) 
            
            if ( ii < 1 .or.  jj < 1 ) cycle
            if (cld_phase ( i_lrc(i,j),j_lrc(i,j))  == ET_cloud_phase % CLEAR ) cycle
                        
            cld_phase ( i,j) = cld_phase ( i_lrc(i,j),j_lrc(i,j))
            
            !- ice pixels where lrc is water ==> retype to water
            if (cld_type ( i,j) == ET_cloud_type % CIRRUS &
               .or. cld_type ( i,j) == ET_cloud_type % OVERLAP) then
               
               if (  cld_type (ii,jj) >= ET_cloud_type % FIRST_WATER &
                     .and. cld_type (ii,jj) <= ET_cloud_type % LAST_WATER) then
               
                  cld_type (i , j ) = cld_type (ii , jj )
               
               end if                    
            end if   
           
           
            !- water clouds at edge of cirrus
            if (cld_type ( i,j) == ET_cloud_type % FOG &
               .or. cld_type ( i,j) == ET_cloud_type % WATER) then
               
               if (cld_type (ii,jj) == ET_cloud_type % CIRRUS &
                  .or. cld_type ( i,j) == ET_cloud_type % OVERLAP &
                  .or. cld_type ( i,j) == ET_cloud_type % OPAQUE_ICE) then
                  
                  cld_type (i , j ) = ET_cloud_type % CIRRUS   
               end if   
            end if
            
            
            
         end do elem_loop_lrc
      end do   line_loop_lrc
      
      deallocate ( phase_inp % rtm % p_prof )
      
      
      
   end subroutine cloud_type_bridge

end module cloud_type_bridge_module
