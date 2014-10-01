#!/usr/bin/env bash

function usage() {

cat<< EOF

CG_GET_VIIRS_DATA

need help? bad luck,  sorry!
ask denis.botambekov@ssec.wisc.edu or andi.walther@ssec.wisc.edu 

EOF

}

while :; do 
   case $1 in
   
   --reg)
      if [ "$2" ]; then
         REG=$2
         shift 2
         continue
      fi
      ;;
      *)
      break
   esac
done

n=$#
echo "==> "$n
TIME0=$1

length0=$(echo ${#TIME0})
echo $length

YEAR0=${TIME0:0:4}
MONTH0=${TIME0:4:2}

DAY0='01'
if [ $length0 -gt 6 ]; then 
   DAY0=${TIME0:6:2}
fi   

HOUR0='00'
if [ $length0 -gt 8 ]; then
 HOUR0=${TIME0:8:2}
fi 

MINU0='00'
if [ $length0 -gt 10 ]; then 
   MINU0=${TIME0:10:2}
fi

YEAR1=$YEAR0
MONTH1=$MONTH0
DAY1=$DAY0
if [ $length0 -le 6 ]; then 
   DAY1='31'
fi


HOUR1=$HOUR0
if [ $length0 -le 8 ]; then 
   HOUR1='23'
fi


MINU1=$MINU0
if [ $length0 -le 10 ]; then 
   MINU1='59'
fi

if [ $n -eq 2 ]; then
   TIME1=$2
   length1=$(echo ${#TIME1})
   YEAR1=${TIME1:0:4}
   MONTH1=${TIME1:4:2}
   DAY1='31'
   if [ $length1 -gt 6 ]; then 
      DAY1=${TIME1:6:2}
   fi
      
   HOUR1='23'
   if [ $length1 -gt 8 ]; then 
      HOUR1=${TIME1:8:2}
   fi
      
   MINU1='59'
   if [ $length1 -gt 10 ]; then 
      MINU1=${TIME1:10:2}
   fi   
fi



if  [ $REG ]; then
   

   case $REG in
      spc)
         echo "south america"
         ll_lat=-45
         ll_lon=-100
         ur_lat=15
         ur_lon=-60
      ;;
      *)
      echo "unknown area"
   
   esac
   
   
   
   sh -c './cg_peate_downloader.sh --ll '$ll_lat' '$ll_lon' --ur '$ur_lat' '$ur_lon' '$YEAR0'-'$MONTH0'-'$DAY0'+'$HOUR0':'$MINU0':00 '$YEAR0'-'$MONTH0'-'$DAY1'+'$HOUR1':'$MINU1':00 SVDNB GMTCO GITCO GDNBO SVM01 SVM02 SVM03 SVM04 SVM05 SVM06 SVM07 SVM08 SVM09 SVM10 SVM11 SVM12 SVM13 SVM14 SVM15 SVM16 IICMO SVI01'
   
else
   echo "global: "  
   sh -c './cg_peate_downloader.sh '$YEAR0'-'$MONTH0'-'$DAY0'+'$HOUR0':'$MINU0':00 '$YEAR0'-'$MONTH0'-'$DAY1'+'$HOUR1':'$MINU1':00 SVDNB GMTCO GITCO GDNBO SVM01 SVM02 SVM03 SVM04 SVM05 SVM06 SVM07 SVM08 SVM09 SVM10 SVM11 SVM12 SVM13 SVM14 SVM15 SVM16 IICMO SVI01' 
fi




