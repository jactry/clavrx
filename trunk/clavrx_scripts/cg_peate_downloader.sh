#!/usr/bin/env bash

function usage() {

cat<< EOF


need help? bad luck,  sorry!
ask denis.botambekov@ssec.wisc.edu or andi.walther@ssec.wisc.edu 

EOF

}



while :; do
   case $1 in 
      -h|-\?|--help) 
         usage
         exit
         ;;
      --lon)
         if [ "$2" ]; then
            LON=$2
            echo $LON
            
            if [ "$LON" -lt -180 ] || [ "$LON" -gt 180 ]; 
            then
               echo 'ERROR: Longitude not set or in wrong range! '
               usage
               exit 1
            fi 
            
            shift 2
            continue 
         else
            echo 'ERROR: Must specify a non-empty "--lon lon" argument.' >&2
            exit 1
         fi       
      ;;
       --lat)
         if [ "$2" ]; then
            LAT=$2 
            
            if [ "$LAT" -lt -90 ] || [ "$LAT" -gt 90 ];  
            then
               echo 'ERROR: Latitude not set or in wrong range! '
               usage
               exit 1
            fi
            
            shift 2
            continue 
         else
            echo 'mist specify long'
            exit 1
         fi       
      ;;
      --radius)
         RAD=$2
         
         if [ "$RAD" -le 0 ] || [ "$RAD" -gt 5000 ];  
            then
               echo 'WARNING: Radius not in the correct range; set to 2000 '
               RAD=500
            fi
         
         shift 2
         continue
      ;;
      --ll)
         LL_LAT=$2
         LL_LON=$3
         echo "  ==> ll < ==="
         echo $LL_LAT
         echo $LL_LON
         if [ "$LL_LAT" -le -90 ] || [ "$LL_LAT" -gt 90 ] || [ "$LL_LON" -lt -180 ] || [ "$LL_LON" -gt 180 ];  
            then
               echo 'WARNING: Box not in the correct range; set to 2000 '
               usage
               exit 1
            fi
         
         shift 3
         continue
      ;;
      --ur)
         UR_LAT=$2
         UR_LON=$3
         
         echo "ur "
         echo $UR_LAT
         echo $UR_LON
         
         if [ "$UR_LAT" -le -90 ] || [ "$UR_LAT" -gt 90 ] || [ "$UR_LON" -lt -180 ] || [ "$UR_LON" -gt 180 ];  
            then
               echo 'WARNING: Box not in the correct range; set to 2000 '
               usage
               exit 1
            fi
         
         shift 3
         continue
      ;;
      
      --day)
         DAY='1'
         
      ;;
      
      --night)
         NIGHT='1'
         
      ;;
      
      
      *)
      break
   esac
   shift
done


if ([ $LON ] && [ ! $LAT ]) || ([ ! $LON ] && [ $LAT ]); 
then
  echo 'ERROR: Both, LAT and LON must be set '
  usage
  exit 1
fi

START=$1
END=$2
shift 2
TYPES=$*

echo $TYPES

OUTPUT=wget
URL="http://peate.ssec.wisc.edu/flo/api/find?start=$START&end=$END&output=$OUTPUT"

for ft in $TYPES; do
	URL="$URL&file_type=$ft"
done

if [ $LON ];
   then
   URL="$URL&loc=$LAT,$LON"
   if [ $RAD ];
      then
      URL="$URL&radius=$RAD"
   else
      URL="$URL&radius=500"   
   fi
fi  


if [ $UR_LAT ] && [ $LL_LAT ];
   then
   URL="$URL&ll=$LL_LAT,$LL_LON&ur=$UR_LAT,$UR_LON"
fi

if [ $DAY ];
   then
   URL="$URL&tod=D"
fi 

if [ $NIGHT ];
   then
   URL="$URL&tod=N"
fi 



echo $URL

SCRIPT=downloader.sh

wget -q -O $SCRIPT ${URL}

# make wget look if data exist
old_strg="-q"
new_strg="-q -nc"
sed "s/$old_strg/$new_strg/g"<$SCRIPT >'tmp.txt'
mv 'tmp.txt' $SCRIPT

bash $SCRIPT

echo "sucess"
exit


