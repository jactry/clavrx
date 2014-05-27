#!/bin/sh
# $Id: run_viirs_exact.sh,v 1.2 2013/04/18 19:17:00 awalther Exp $

args=("$@")
year=${args[0]} 
doy=${args[1]} 
hour=${args[2]}
minu=${args[3]}
dcomp_mode=${args[4]}
shift 6

satname="viirs"

echo $HOME


algo_path="main_src/"
cd $algo_path
cd ./../
pwd
options='clavrxorb_default_options'
l1b_path=$HOME'/Satellite_Input/'$satname'/'$year'/'$doy'/'
out_path=$HOME'/Satellite_Output/'$satname'/'$year'/'$doy'/'
[ ! -d $l1b_path ] && mkdir -v -p $l1b_path
[ ! -d $out_path ] &&  mkdir -v -p $out_path


dummy=`grep !nwp $options`

nwp=${dummy:0:1}

case $nwp in
1)echo "NWP option from clavrxorb file : GFS" 
     ./get_gfs.sh $year $doy;;
2)echo "ncep" 
     ./get_sfcr_data.sh $year $doy;;
3)echo "NWP option from clavrxorb file : CFSR" 
     ./get_cfsr.sh $year $doy;;
esac
echo


pwd
 if date -v 1d > /dev/null 2>&1; then
  doy_dum=`expr ${doy} - 1` 

  DAY_OBS=$(./shift_date.sh ${year}0101 +${doy_dum})
  month=${DAY_OBS:4:2}
  day=${DAY_OBS:6:2}
else
  month=$(date -d "01/01/${year} +${doy} days -1 day" "+%m")
  day=$(date -d "01/01/${year} +${doy} days -1 day" "+%d")
  
fi
echo $month, $day
echo "======", $minu,'==='

pwd
echo $year $doy $hour $minu $l1b_path
./get_viirs_data.sh $year $doy $hour $minu $l1b_path
filetype='GMTCO'
./write_filelist.sh $l1b_path $out_path $month$day e$hour$minu $filetype
pwd
./clavrxorb  -dcomp_mode $dcomp_mode
 #rm -fv $l1b_path/*

