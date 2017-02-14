pro aw_image_fortran


!quiet = 1
args = command_line_args()

args=long(args)


print
print,'args: ',args
read_matrix,'aw_image_temp.dat',l
help,l
data=reform(l,args[0],args[1])
read_matrix,'aw_image_temp_lat.dat',l
help,l
lat=reform(l,args[0],args[1])
help,lat
read_matrix,'aw_image_temp_lon.dat',l
lon=reform(l,args[0],args[1])

print,max(data)
!p.multi=0
oo=aw_plot()
outfile = 'aw_image_fortran_'+string(1000*randomu(seed),form='(i3.3)')+(str_sep(systime(),' '))[3]
;oo.startdevice,out=outfile ;,xs=0.5,ys=.7

p_orig = !p


window,xsize=800,ysize=800,/pixmap

fub_image,data,lon,lat,charsize=0.8, void_index =where(data eq -999.),bar_charsize=1.2, $
    format='(f5.0)',min=0,max=max(data),title='test',legend_title = 'CHANNEL REFLECTANCE '
;oo.enddevice
;oo.makejpg

write_jpeg,outfile+'.jpg',tvrd(true=3),true=3

spawn, 'display -resize 1080x824 '+outfile+'.jpg&'

end
