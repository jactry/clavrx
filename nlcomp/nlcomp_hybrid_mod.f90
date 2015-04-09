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
      , rad_abv_cld &
      , pixel &
      , lut_path &
      , emis_cld_ch31, beta_cld_e31_e32 &
      , cps, tsfc,rad20,rad31,rad32 )
      
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
      real :: rad_abv_cld (42)
      type ( lut_output) :: lut_data20, lut_data31, lut_data32
      integer :: phase_num
      real :: cod_used, rad20_sfc,rad31_sfc,rad32_sfc
      integer :: i , j
      real :: planck_rad20 , planck_rad31 , planck_rad32
      real :: fwd(3,9)
      integer :: cod_pos(3) = -999
      real :: cps_wgt
      
      
      cps_vec_lut = (/(i,i = 4 , 20 , 2 )/)/10.
 
      planck_rad20 = planck_tmp2rad ( pixel % ctt, 'VIIRS' , 20)
      planck_rad31 = planck_tmp2rad ( pixel % ctt, 'VIIRS' , 31)
      planck_rad32 = planck_tmp2rad ( pixel % ctt, 'VIIRS' , 32)
 
      rad20_sfc = planck_tmp2rad ( tsfc , 'VIIRS' , 20)
      rad31_sfc = planck_tmp2rad ( tsfc , 'VIIRS' , 31)
      rad32_sfc = planck_tmp2rad ( tsfc , 'VIIRS' , 32)
      
      bt20 =  planck_rad2tmp ( planck_rad20, 'VIIRS' , 20)
      bt31 =  planck_rad2tmp ( planck_rad31, 'VIIRS' , 31)
      bt32 =  planck_rad2tmp ( planck_rad32, 'VIIRS' , 32)
      
      
      print*,'vergleich PLanck sfc20, rad_clear_toc20: ', rad20_sfc , rad_clear_toc(20)
      
      print*
      print*,'channel 20:'
      print*,'RAD: ', rad20
      print*,rad20_sfc
      print*,rad_clear_toc(20)
      print*,'cld rad: ', planck_rad20
      print*,'rad abv: ',rad_abv_cld(20)
      print*,'cld BT: ', bt20
      print*
       print*,'channel 31:'
       print*,'RAD: ', rad31
      print*,rad31_sfc
      print*,rad_clear_toc(31)
       print*,'cld rad: ', planck_rad31
      print*,'rad abv: ',rad_abv_cld(31)
       print*,'cld BT: ', bt31
      print*
       print*,'channel 32:'
       print*,'RAD: ', rad32
      print*,rad32_sfc
      print*,rad_clear_toc(32)
       print*,'cld rad: ', planck_rad32
      print*,'rad abv: ',rad_abv_cld(32)
       print*,'cld BT: ', bt32
       
      print*
      print*, 'CTT: ',pixel % ctt
      print*, 'rad20,rad31, rad32: ', planck_rad20, planck_rad31, planck_rad32      
      
 
      
      print*,'bt: ', bt20, bt31, bt32
      phase_num = 2
      if ( pixel % is_water_phase ) phase_num = 1
      
      
      call lut_obj % initialize ( 'VIIRS' , ancil_path = lut_path)
      call lut_obj % set_angles ( pixel % sat_zen , pixel % lun_zen &
            , pixel % lun_rel_azi )
      
       print*,'obs: ',obs
       print*, 'cod: ',cod
      
       
      cod_used = cod
      if ( cod < -100. ) cod_used = 1.02

     
    !  do i = 1, 9          
    !     call lut_obj % get_data ( 31, phase_num, cod_used , cps_vec_lut (i) , lut_data31)
    !     call lut_obj % get_data ( 32, phase_num, cod_used , cps_vec_lut (i) , lut_data32)
    !     print*,'CPS: ems,beta: ',i,lut_data31 %  ems ,  log(lut_data31 %  ems )/ log(lut_data32 %  ems)
         !rfl_toa(i) =  lut_data%refl + lut_data%trn_sol* lut_data%trn_sat * alb_sfc &
         !   / (1- alb_sfc * lut_data % albsph)
         !if ( i > 1 ) then
         !    if ( rfl_toa(i) .gt. rfl_dnb .and. rfl_toa(i-1) .lt. rfl_dnb ) cps_pos = i 
         !end if        
         ! print*,i, phase_num, rfl_toa(i), lut_data%refl
    !  end do  
      
      print*
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
         
         fwd(1,i) = rad20
         fwd(2,i) = bt31-bt32
         fwd(3,i) = bt20 - bt31
         print*,i,'trn ==> ',lut_data20 % trn_ems,lut_data31 % trn_ems,lut_data32 % trn_ems
         print*,i, 'rad FWD==> ', rad20,rad31,rad32
         print*,i,' FWD==> ', rad20,bt31-bt32 , bt20 - bt31
         print*,i,' OBS==> ', obs(2), obs(3),obs(4)
         print*
         
         if ( i > 1 ) then
             if ( fwd(1,i) .gt. obs(2) .and. fwd(1,i-1) .lt. obs(2) .or. fwd(1,i) .lt. obs(2) .and. fwd(1,i-1) .gt. obs(2) ) cod_pos(1) = i 
             if ( fwd(2,i) .gt. obs(3) .and. fwd(2,i-1) .lt. obs(3) .or. fwd(2,i) .lt. obs(3) .and. fwd(2,i-1) .gt. obs(3) ) cod_pos(2) = i 
             if ( fwd(3,i) .gt. obs(4) .and. fwd(3,i-1) .lt. obs(4) .or. fwd(3,i) .lt. obs(4) .and. fwd(3,i-1) .gt. obs(4) ) cod_pos(3) = i 
             
         end if 
         
      end do  
      
      do j = 1 ,3   
         if ( cod_pos(j) > 1 ) then
            cps_wgt = ( fwd ( j, cod_pos(j) ) - obs(j+1) ) / ( fwd ( j, cod_pos(j) - 1 ) -  fwd ( j, cod_pos(j) ))
            cps = cps_vec_lut ( cod_pos(j) ) + cps_wgt * 0.2
            print*,j,10.**cps 
         end if
      end do
      
      
      print*,cod_pos , cps, 10**cps   
      !end do  
   
   end subroutine cps_known_cod


end module nlcomp_hybrid_mod
