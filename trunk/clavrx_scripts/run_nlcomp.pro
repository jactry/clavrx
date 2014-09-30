pro run_nlcomp $
	, year = year $
	, doy0 = doy0 $
	, doy1 = doy1 $
	, ancil_path = ancil_path $
	, region = region 

default, year, 2013
default,doy0, 1
default,doy1, 366
default,region,'global'
default, ancil_path, '/data3/Ancil_Data/clavrx_ancil_data/dnb_ancils/'

month = 3


lunar_irrad_file = ancil_path+'dnb_ancils/lunar_irrad_Mean_DNB.bin'
distance_table_file= ancil_path+'dnb_ancils/DIST_2010-2030_double.bin'


case region of
	'global': begin
		lon0 = -180
		lon1 = 180
		lat0 = -90
		lat1 = 90.
	end
	
	'south_amer': begin
			lon0 = -100
			lon1 = -40
			lat0 = -50
			lat1 = 10.
	end
	
	else: begin
		print,'region is wrong!'
    end


endcase


thr  = 60.





for doy = doy0, doy1 do begin
   jd = julday ( 1, 1 , year ) + doy - 1
   caldat, jd , month , day_of_month
   for hour = 0, 23 do begin
      lunar_phase = lunar_phase_lut ( year , month, day_of_month , hour )
      print,year,doy,lunar_phase
      if lunar_phase gt 110 then continue
      for m =0, 59 do begin
	     ; error =''
	     ; catch,error
	     ; if error ne '' then begin
        ;    print, 'error ',d,h,m
	     ;    continue	
	    ;  endif
         
        
         
         timestring = string(year,form='(i4.4)')+string(month,form='(i2.2)') + string(day_of_month,form='(i2.2)') $
               + string(hour,form='(i2.2)') $
               + string(m,form='(i2.2)')+'00'
             
         lunar_zen = read_viirs(timestring , 'LUN_ZEN',/dnb_flag)
         print,timestring+'  number of lunar zenith below 50: ' , total(between ( lunar_zen, 0,thr))
         
         if total(between ( lunar_zen, 0 , thr)) lt 300000 then continue
         dnb = read_viirs(timestring , 'DNB',sol_zen = sol_zen,lon=lon,lat=lat,/dnb_flag)
         n_dnb = total (between ( lunar_zen, 0,thr) and sol_zen gt 100 and dnb gt 1.0e-8 and between(lon,lon0,lon1) and between(lat,lat0,lat1)  )
         print, 'count2: ', n_dnb, total ( sol_zen gt 100 ) , total ( dnb gt 1.0e-8 )
         
         if n_dnb gt 10000 then begin
            spawn,'pwd'
            
            unix_str = './run_viirs_exact.sh  '+string(year,format='(i4.4)')+' '+string(doy,format='(i3.3)')+' '+string(hour,format='(i2.2)') $
                            + ' '+ string ( m,form='(i2.2)') +' 0'
            print,unix_str              
            spawn,unix_str
	         catch,/cancel
         endif
         
     endfor
   endfor
endfor



end
