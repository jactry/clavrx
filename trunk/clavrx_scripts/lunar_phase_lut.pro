




function lunar_phase_lut , year , month, day , hour 


num_dist_tabvals = 184080L 

default, ancil_path, '/DATA/Ancil_Data/clavrx_ancil_data/dnb_ancils/'
distance_table_file= ancil_path+'DIST_2010-2030_double.bin'
yyyymmddhh = year * 1000000 + month*10000+ day *100 + hour

data_base = dblarr( 4 * num_dist_tabvals )
openr,lun,distance_table_file,/get_lun

readu,lun,data_base

free_lun,lun
data_base = reform(data_base,[ 4 , num_dist_tabvals ] )
idx = value_locate (  data_base[0,*],yyyymmddhh) 


lunar_phase = data_base[1,idx] 


return,lunar_phase



end
