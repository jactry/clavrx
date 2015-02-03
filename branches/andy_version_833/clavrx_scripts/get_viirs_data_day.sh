#!/bin/sh
# $Id: get_viirs_data.sh,v 1.7 2013/04/18 16:36:37 awalther Exp $

args=("$@") 
year=${args[0]} 
doy=${args[1]}
path=${args[2]}

 if date -v 1d > /dev/null 2>&1; then
  doy_dum=`expr ${doy} - 1` 

  DAY_OBS=$(./shift_date.sh ${year}0101 +${doy_dum})
  month=${DAY_OBS:4:2}
  day=${DAY_OBS:6:2}
else
  month=$(date -d "01/01/${year} +${doy} days -1 day" "+%m")
  day=$(date -d "01/01/${year} +${doy} days -1 day" "+%d")
  
fi
echo $month
echo $day

cp peate_downloader.sh $path
cd $path
 
sh -c './peate_downloader.sh '$year'-'$month'-'$day'+00:00:00 '$year'-'$month'-'$day'+23:59:59 SVDNB GMTCO GITCO GDNBO SVM01 SVM02 SVM03 SVM04 SVM05 SVM06 SVM07 SVM08 SVM09 SVM10 SVM11 SVM12 SVM13 SVM14 SVM15 SVM16 IICMO SVI01'


exit
