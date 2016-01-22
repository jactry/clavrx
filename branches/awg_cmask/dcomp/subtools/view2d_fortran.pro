pro view2d_fortran


!quiet = 1
args = command_line_args()

args=long(args)


l=make_array([args[0],args[0]],value=0)

read_matrix,'view2d_temp.dat',l
l=reform(l,100,100)

!p.multi=0
oo=aw_plot()
oo.startdevice,out='view2d_fortran_temp' ;,xs=0.5,ys=.7
view2d,l,/cool,/colo,charsize=0.8, no_data_val =-999.
oo.enddevice
oo.makejpg

spawn, 'display -resize 1280x1024 view2d_fortran_temp.jpg'

end
