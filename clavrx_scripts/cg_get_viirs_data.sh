#!/usr/bin/env bash

function usage() {

cat<< EOF

CG_GET_VIIRS_DATA

need help? bad luck,  sorry!
ask denis.botambekov@ssec.wisc.edu or andi.walther@ssec.wisc.edu 

EOF

}




BASE_PATH=$HOME"/Satellite_Input/viirs/"
satname='viirs'

while :; do 
   case $1 in
   
   --reg)
      if [ "$2" ]; then
         REG=$2
         shift 2
         continue
      fi
      ;;
   --path)
      if [ "$2" ]; then
         L1B_PATH=$2
         shift 2
         continue
      fi
   ;;   
      *)
      break
   esac
done

n=$#

YEARDOY=$1

YEAR=${YEARDOY:0:4}
DOY=${YEARDOY:4:3}

HOUR0='00'
MINU0='00'


if [ $n -gt 1 ]; then
   TIME0=$2 
   length0=$(echo ${#TIME0}) 
   HOUR0=${TIME0:0:2}
   if [ $length0 -gt 2 ]; then 
      MINU0=${TIME0:2:2}
   fi   
fi 

HOUR1=$HOUR0
MINU1=$MINU0

if [ $n -gt 2 ]; then
   TIME1=$3 
   HOUR1='23'
   MINU1='59'
   length1=$(echo ${#TIME1}) 
   HOUR1=${TIME1:0:2}
   if [ $length1 -gt 2 ]; then 
      MINU1=${TIME1:2:2}
   fi   
fi 



if date -v 1d > /dev/null 2>&1; then
  doy_dum=`expr ${DOY} - 1` 
  DAY_OBS=$(./shift_date.sh ${YEAR}0101 +${DOY})
  MONTH=${DAY_OBS:4:2}
  DAY=${DAY_OBS:6:2}
else
  MONTH=$(date -d "01/01/${YEAR} +${DOY} days -1 day" "+%m")
  DAY=$(date -d "01/01/${YEAR} +${DOY} days -1 day" "+%d") 
fi



mkdir -p $l1b_path


if  [ $REG ]; then
   

   case $REG in
      spc)
         echo "south america pacific"
         ll_lat=-45
         ll_lon=-100
         ur_lat=15
         ur_lon=-60
      ;;
      
      kzn)
         echo "kazachstan"
         ll_lat=20
         ll_lon=90
         ur_lat=50
         ur_lon=130
      ;;   
      *)
      echo "unknown area"
   
   esac
   
else

   echo "global: " 
   
         ll_lat=-80
         ll_lon=-170
         ur_lat=80
         ur_lon=170
         
         
      REG='global'
  
fi

if [ ! $L1B_PATH ];
then

	L1B_PATH=$HOME'/Satellite_Input/'$satname'/'$YEAR'/'$DOY'/'$REG
fi

mkdir -p $L1B_PATH
echo $l1b_path



      sh -c './cg_peate_downloader.sh --path '$l1b_path' --ll '$ll_lat' '$ll_lon' --ur '$ur_lat' '$ur_lon' '$YEAR'-'$MONTH'-'$DAY'+'$HOUR0':'$MINU0':00 '$YEAR'-'$MONTH'-'$DAY'+'$HOUR1':'$MINU1':00 SVDNB GMTCO GITCO GDNBO SVM01 SVM02 SVM03 SVM04 SVM05 SVM06 SVM07 SVM08 SVM09 SVM10 SVM11 SVM12 SVM13 SVM14 SVM15 SVM16 IICMO SVI01'


