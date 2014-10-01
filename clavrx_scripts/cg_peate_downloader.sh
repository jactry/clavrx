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
            
            shift 2
            continue 
         else
            echo 'mist specify long'
            exit 1
         fi       
      ;;
      --region)
         REG=$2
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
   if [ $REG ];
      then
      URL="$URL&radius=$REG"
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


