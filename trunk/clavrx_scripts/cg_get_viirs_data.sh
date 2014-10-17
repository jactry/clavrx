#!/usr/bin/env bash
#
#  $Id:$
#
function usage() {

cat<< EOF

CG_GET_VIIRS_DATA

need help? bad luck,  sorry!
ask denis.botambekov@ssec.wisc.edu or andi.walther@ssec.wisc.edu 

EOF

}

# calles by check_download
function redo_it () {
#  echo "Check files again"

sch_path=$l1b_path'*h5'
num_real_files=`ls -1 $sch_path | wc -l`
echo "Real $num_real_files"

if [ "$num_real_files" -ne "$num_need_files" ]
then
   echo "Not enough files"
   cd $l1b_path
   bash $down_file
   stat=1
else
   stat=0
fi
return $stat
}

#
#this is needed because peate scripts sometimes doesn't download all files
function check_download() {
args=("$@") 
l1b_path=${args[0]}
down_file=$l1b_path'downloader.sh'
curr_dir=`pwd`

tmp=`grep "files" $down_file`

if [ -f $down_file ] ; then
  num_need_files=`echo $tmp |awk -F " " '{print $2}'`
  echo "Need $num_need_files"
  sed -e 's/-q/-q -nc/g' <$down_file >temp.txt
  cp temp.txt $down_file
  rm temp.txt
  stat=1
fi

while [ $stat -eq 1 ]
do
   redo_it
done

cd $curr_dir
echo "All files are there"

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
echo $YEARDOY
YEAR=${YEARDOY:0:4}
DOY=${YEARDOY:4:3}

echo $YEAR,$DOY' <===='


HOUR0='00'
MINU0='00'
 MINU1='59'

if [ $n -gt 1 ]; then
   TIME0=$2 
   length0=$(echo ${#TIME0}) 
   HOUR0=${TIME0:0:2}
   if [ $length0 -gt 2 ]; then 
      MINU0=${TIME0:2:2}
	  MINU1=$MINU0
   fi   
fi 

HOUR1=$HOUR0


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






if  [ $REG ]; then
   

   case $REG in
      spc)
         echo "south america pacific"
         ll_lat=-45
         ll_lon=-100
         ur_lat=15
         ur_lon=-60
      ;;
	  
	  sal)
	     echo "south atlantic"
         ll_lat=-45
         ll_lon=-30
         ur_lat=15
         ur_lon=15
      ;;	
      
      kaz)
         echo "kazachstan"
         ll_lat=40
         ll_lon=40
         ur_lat=55
         ur_lon=90
      ;;   
	  
	  bal)
         echo "baltic"
         ll_lat=55
         ll_lon=5
         ur_lat=65
         ur_lon=25
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
echo $L1B_PATH



      sh -c './cg_peate_downloader.sh --path '$L1B_PATH' --ll '$ll_lat' '$ll_lon' --ur '$ur_lat' '$ur_lon' '$YEAR'-'$MONTH'-'$DAY'+'$HOUR0':'$MINU0':00 '$YEAR'-'$MONTH'-'$DAY'+'$HOUR1':'$MINU1':00 SVDNB GMTCO GITCO GDNBO SVM01 SVM02 SVM03 SVM04 SVM05 SVM06 SVM07 SVM08 SVM09 SVM10 SVM11 SVM12 SVM13 SVM14 SVM15 SVM16 IICMO SVI01'

check_download $L1B_PATH
