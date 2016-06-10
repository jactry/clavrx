;  this tool creates in include file for clavr-x which 
;  assigns clavr-x internal (andglobal) variables to output in csv file
;

pro make_level2_inc

csv_file= 'clavrx_level2_products.csv'
data = read_csv(csv_file)

file_l2 = 'level2_assign.inc'
openw,10,file_l2
printf,10,'select case(trim(name))'
for i=0,n_elements(data.field01) -1  do begin
   if (data.(1))[i] eq '_filename' then continue
   if (data.(1))[i] eq '_global_attr' then continue
   if trim((data.(2))[i]) eq 'NOT_SET_YET' then continue
   
   out_name = (data.(1))[i]
   printf,10,'case("'+out_name+'")'
   var_dim = (data.(3))[i]
   dtype = (data.(4))[i]
   global_var = (data.(2))[i]
   scaling = (data.(5))[i]
   if global_var eq 'special' then continue
   
   sub = '(Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)'
   if var_dim eq 2 then sub = '(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)'
   data_name = 'data_dim'+string(var_dim,for='(i1)')+'_dtype'$
         +string(dtype,for='(i1)')
   if scaling eq 1 then data_name = 'data_dim' +string(var_dim,for='(i1)')+'_dtype_r4'     
   ; special case where three different output variables come from one clavr-x var
   if out_name eq 'cld_temp_acha_qf' then begin
      global_var = 'acha%oe_quality_flags'
      sub = '(1,:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)'
   endif   
   if out_name eq'cld_emiss_acha_qf' then begin
      global_var = 'acha%oe_quality_flags'
       sub = '(2,:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)'
   endif    
   if out_name eq 'cld_beta_acha_qf' then begin
      global_var = 'acha%oe_quality_flags'
      sub = '(3,:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)'
   endif   
         printf,10,'   if (allocated ( '+global_var+')) then'
         printf,10,'      '+data_name+' = '+global_var+sub
         printf,10,'   end if'
    
    
    
    
endfor

printf,10,'end select'


close,10




end
