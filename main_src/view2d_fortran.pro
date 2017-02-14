pro view2d_fortran


!quiet = 1
args = command_line_args()

args=long(args)




read_matrix,'view2d_temp.dat',l
help,l
print,args
l=reform(l,args[0],args[1])

!p.multi=0
oo=aw_plot()
outfile = 'view2d_fortran_'+string(1000*randomu(seed),form='(i3.3)')+(str_sep(systime(),' '))[3]
oo.startdevice,out=outfile ;,xs=0.5,ys=.7
view2d,l,/cool,/colo,charsize=0.8, no_data_val =-999.
oo.enddevice
oo.makejpg

spawn, 'display -resize 1280x1024 '+outfile+'.jpg&'

end
