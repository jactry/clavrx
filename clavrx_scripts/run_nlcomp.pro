pro run_nlcomp

year = 2013
month = 3

ancil_data_dir = '/data3/Ancil_Data/clavrx_ancil_data/dnb_ancils/'
      lunar_irrad_file = ancil_data_dir+'dnb_ancils/lunar_irrad_Mean_DNB.bin'
      distance_table_file= ancil_data_dir+'dnb_ancils/DIST_2010-2030_double.bin'

lon0 = -180
lon1= 180
lat0=-90
lat1=90.
thr = 60.

for d =29, 29 do begin
  for h = 6, 7 do begin
   ;if h gt 8 then continue
     for m =0, 59 do begin
	      error =''
	      catch,error
	      if error ne '' then begin
            print, 'error ',d,h,m
	         continue	
	      endif
         
         timestring = string(year,form='(i4.4)')+string(month,form='(i2.2)') + string(d,form='(i2.2)') $
               + string(h,form='(i2.2)') $
               + string(m,form='(i2.2)')+'00'
             
         lunar_zen = read_viirs(timestring , 'LUN_ZEN',/dnb_flag,lon=lon,lat=lat)
         print,timestring+'  number of lunar zenith below 50: ' , total(between ( lunar_zen, 0,thr))
         if total(between ( lunar_zen, 0 , thr)) lt 300000 then continue
         dnb = read_viirs(timestring , 'DNB',sol_zen = sol_zen,/dnb_flag)
         n_dnb = total (between ( lunar_zen, 0,thr) and sol_zen gt 100 and dnb gt 1.0e-8 and between(lon,lon0,lon1) and between(lat,lat0,lat1)  )
         print, 'count2: ', n_dnb, total ( sol_zen gt 100 ) , total ( dnb gt 1.0e-8 )
         
         if n_dnb gt 10000 then begin
            spawn,'pwd'
            doy = 59 + d
            unix_str = './run_viirs_exact.sh 2013 '+string(doy,format='(i3.3)')+' '+string(h,format='(i2.2)') $
                            + ' '+ string ( m,form='(i2.2)') +' 0'
            print,unix_str              
            spawn,unix_str
	         catch,/cancel
         endif
         
     endfor
   endfor
endfor



end
