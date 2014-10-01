#!/usr/bin/env bash

function usage() {

cat<< EOF


need help? bad luck sorry!

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
         shift 2
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


if ([ $LON ] && [ ! $LAT]) || ([ ! $LON ] && [ $LAT ]); 
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
      URL="$URL&radius=2000"   
   fi
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
bash $SCRIPT

echo "sucess"
exit


