;  this tool creates in include file for clavr-x which 
;  assigns clavr-x internal (andglobal) variables to output in csv file
;

pro make_level2_inc, try_alloc=try_alloc

csv_file= 'clavrx_level2_products.csv'
data = read_csv(csv_file)

file_l2 = 'level2_assign.inc'
openw,10,file_l2
printf,10,'select case(trim(name))'
for i=0,n_elements(data.field01) -1  do begin
   if (data.(2))[i] eq '_filename' then continue
   if (data.(2))[i] eq '_global_attr' then continue
   
   printf,10,'case("'+(data.(2))[i]+'")'
   var_dim = (data.(1))[i]
   dtype = (data.(4))[i]
   global_var = (data.(3))[i]
   if global_var eq 'special' then continue
   
   sub = '(Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)'
   if var_dim eq 2 then sub = '(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)'
   data_name = 'data_dim'+string(var_dim,for='(i1)')+'_dtype'$
         +string(dtype,for='(i1)')
   if ~keyword_set(try_alloc) then begin
      
         printf,10,'   if (allocated ( '+global_var+')) then'
         printf,10,'      '+data_name+' = '+global_var+sub
         printf,10,'   end if'
    endif else begin
      printf,10,'allocate('+data_name+' ,source='+global_var+sub+')'
    
    
    endelse     
endfor

printf,10,'end select'


close,10




end
