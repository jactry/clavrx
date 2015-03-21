! $Id:$
!
!
module nlcomp_hybrid_mod
   
   use dcomp_lut_mod , only : &
       lut_obj , lut_output
   
   use nlcomp_forward_mod, only : &
        
        pixel_vec
        
   public :: vis_channel_cod
   public :: cps_known_cod
   private   

    

contains
   !>
   !!
   !!
   subroutine vis_channel_cod ( &
        rfl_dnb &
      , pixel &
      , alb_sfc &
      , lut_path &
      , cod &
      )
      
      real, intent(in) :: rfl_dnb
      type ( pixel_vec ) , intent ( in) :: pixel
      real , intent(in) :: alb_sfc
      character ( len = 1024 ) , intent ( in ) , optional :: lut_path
      real, intent(out) :: cod
   
      type ( lut_output) :: lut_data
      integer :: i 
      real :: rfl_toa(29)
      integer :: cps_pos
      real :: cod_vec_lut(29)
      real :: cod_wgt
      integer :: phase_num
      
      cps_pos = -999
      cod = -999.
      cod_vec_lut = [(i/10.,i=-6,22,1)]
     
      phase_num = 2
      if ( pixel % is_water_phase ) phase_num = 1
      
      call lut_obj % initialize ( 'VIIRS' , ancil_path = lut_path)
      call lut_obj % set_angles ( pixel % sat_zen , pixel % lun_zen &
            , pixel % lun_rel_azi )   
        
      do i = 1, 29          
         call lut_obj % get_data ( 1, phase_num, cod_vec_lut (i) , 1.0 , lut_data) 
         rfl_toa(i) =  lut_data%refl + lut_data%trn_sol* lut_data%trn_sat * alb_sfc &
            / (1- alb_sfc * lut_data % albsph)
         if ( i > 1 ) then
             if ( rfl_toa(i) .gt. rfl_dnb .and. rfl_toa(i-1) .lt. rfl_dnb ) cps_pos = i 
         end if        
          print*,'COD SEARCH: ',i, phase_num, rfl_toa(i), lut_data%refl
      end do  
 
      
      if ( cps_pos > 1 ) then
         cod_wgt = ( rfl_toa ( cps_pos ) - rfl_dnb ) / (rfl_toa ( cps_pos -1 ) - rfl_toa ( cps_pos  ) )
         cod = cod_vec_lut ( cps_pos ) + cod_wgt * 0.1 
      end if
      
      print*,cps_pos,rfl_dnb,  10**cod
  
   end subroutine vis_channel_cod
   
   
   !>
   !!
   !!
   subroutine cps_known_cod ( &
        obs &
      , cod &
      , rad_clear_toc &
      , pixel &
      , lut_path &
      , emis_cld_ch31, beta_cld_e31_e32 &
      , cps, tsfc )
      
       use clavrx_planck_mod, only: &
         planck_tmp2rad &
         , planck_rad2tmp
      
      implicit none
      real, intent(in) :: obs(4)
      real, intent(in) :: cod
      type ( pixel_vec ) , intent ( inout) :: pixel
      character ( len = 1024 ) , intent ( in ) , optional :: lut_path
      real, intent(in) :: emis_cld_ch31
      real, intent(in) :: beta_cld_e31_e32
      
      real, intent(out) :: cps
      real, intent(in) :: tsfc
      real :: cps_vec_lut(9)
      real :: rad20, rad31, rad32
      real :: bt20, bt31, bt32
      real :: rad_clear_toc (42)
      type ( lut_output) :: lut_data20, lut_data31, lut_data32
      integer :: phase_num
      real :: cod_used, rad20_sfc
      integer :: i , j
      real :: planck_rad20 , planck_rad31 , planck_rad32
      
      
      cps_vec_lut = (/(i,i = 4 , 20 , 2 )/)/10.
 
      planck_rad20 = planck_tmp2rad ( pixel % ctt, 'VIIRS' , 20)
      planck_rad31 = planck_tmp2rad ( pixel % ctt, 'VIIRS' , 31)
      planck_rad32 = planck_tmp2rad ( pixel % ctt, 'VIIRS' , 32)
 
      rad20_sfc = planck_tmp2rad ( tsfc , 'VIIRS' , 20)
      
      print*,'vergleich PLanck sfc20, rad_clear_toc20: ', rad20_sfc , rad_clear_toc(20)
      print*, 'CTT: ',pixel % ctt
      print*, 'rad20,rad31, rad32: ', planck_rad20, planck_rad31, planck_rad32      
      
      bt20 =  planck_rad2tmp ( planck_rad20, 'VIIRS' , 20)
      bt31 =  planck_rad2tmp ( planck_rad31, 'VIIRS' , 31)
      bt32 =  planck_rad2tmp ( planck_rad32, 'VIIRS' , 32) 
      
      print*,'bt: ', bt20, bt31, bt32
      phase_num = 2
      if ( pixel % is_water_phase ) phase_num = 1
      
      
      call lut_obj % initialize ( 'VIIRS' , ancil_path = lut_path)
      call lut_obj % set_angles ( pixel % sat_zen , pixel % lun_zen &
            , pixel % lun_rel_azi )
      
       print*,'obs: ',obs, 'cod: ',cod
      
       
      cod_used = cod
      if ( cod < -100. ) cod_used = 1.02

      
      do i = 1, 9          
         call lut_obj % get_data ( 31, phase_num, cod_used , cps_vec_lut (i) , lut_data31)
         call lut_obj % get_data ( 32, phase_num, cod_used , cps_vec_lut (i) , lut_data32)
         print*,'CPS: ems,beta: ',i,lut_data31 %  ems ,  log(lut_data31 %  ems )/ log(lut_data32 %  ems)
         !rfl_toa(i) =  lut_data%refl + lut_data%trn_sol* lut_data%trn_sat * alb_sfc &
         !   / (1- alb_sfc * lut_data % albsph)
         !if ( i > 1 ) then
         !    if ( rfl_toa(i) .gt. rfl_dnb .and. rfl_toa(i-1) .lt. rfl_dnb ) cps_pos = i 
         !end if        
         ! print*,i, phase_num, rfl_toa(i), lut_data%refl
      end do  
      
      
      do i = 1, 9
         
         call lut_obj % get_data ( 20, phase_num , cod_used, cps_vec_lut (i) , lut_data20) 
         call lut_obj % get_data ( 31, phase_num , cod_used, cps_vec_lut (i) , lut_data31)   
         call lut_obj % get_data ( 32, phase_num , cod_used, cps_vec_lut (i) , lut_data32)
         
         rad20 = lut_data20 %  ems * planck_rad20 + lut_data20 % trn_ems * rad_clear_toc(20)
         rad31 = lut_data31 %  ems * planck_rad31 + lut_data31 % trn_ems * rad_clear_toc(31)
         rad32 = lut_data32 %  ems * planck_rad32 + lut_data32 % trn_ems * rad_clear_toc(32)
         
         bt20 =  planck_rad2tmp ( rad20, 'VIIRS' , 20)
         bt31 =  planck_rad2tmp ( rad31, 'VIIRS' , 31)
         bt32 =  planck_rad2tmp ( rad32, 'VIIRS' , 32)
         
         ! - first thin clouds,We can use emis_acha or beta_acha to be consistent
         ! - later we can do our own ?? Probably not
         
         
         
         print*
        print*, '',  lut_data20 %  ems * planck_rad20 
        print*, 'TRN_EMS * RAD_CLEAR_TOC: ', lut_data20 % trn_ems * rad_clear_toc(20)
        print*,'TRN_EMS: ',lut_data20 % trn_ems
        print*,'EMS: ',lut_data20 %  ems
        
        
         
         print*,i, rad20,bt31-bt32 , bt20 - bt31
         print*,i,obs(2), obs(3),obs(4)
         print*
         
      end do      
      !end do  
   
   end subroutine cps_known_cod


end module nlcomp_hybrid_mod
