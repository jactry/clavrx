#!/bin/sh
# $Id$
args=("$@")
year=${args[0]} 
doy0=${args[1]} 
doy1=${args[2]} 
dcomp_mode=${args[3]}
shift 4
satname="viirs"

region='sam'

case $region in
	"sam")
		lon=-75
		lat=-30;;
	*)
esac

for (( ddd = $doy0; ddd < $doy1; ddd++ ))
do
	l1b_path=$HOME'/Satellite_Input/'$satname'/'$year'/'$ddd'/'$region'/'
	mkdir -p $l1b_path
	out_path=$HOME'/Satellite_Output/'$satname'/'$year'/'$ddd'/'$region'/'
	mkdir -p $out_path
	./get_gfs.sh $year $ddd
	
	./get_viirs_data_regional_day.sh $year $ddd $lat $lon $l1b_path
	./write_filelist.sh $l1b_path $out_path GMTC
	./clavrxorb -dcmp_mode $dcomp_mode
	
done
