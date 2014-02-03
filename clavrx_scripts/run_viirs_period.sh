#!/bin/sh
# $Id$
args=("$@")
year=${args[0]} 
doy0=${args[1]} 
doy1=${args[2]} 
dcomp_mode=${args[3]}
shift 4
satname="viirs"

for ddd in $doy0 $doy1
do
	l1b_path=$HOME'/Satellite_Input/'$satname'/'$year'/'$ddd'/'
	out_path=$HOME'/Satellite_Output/'$satname'/'$year'/'$ddd'/'
	mkdir -p $out_path
	./get_gfs.sh $year $ddd
	./write_filelist.sh $l1b_path $out_path GMTC
	./clavrxorb -dcmp_mode $dcomp_mode

done
