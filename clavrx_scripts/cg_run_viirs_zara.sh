#!/bin/sh
# $Id:$


# makes hourly jobs stored in temporary shell scripts 


source /etc/bashrc
module load bundle/basic-1
module load bundle/basic-1 hdf5



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
         BASE_PATH=$2
         shift 2
         continue
      fi
   ;;   
      *)
      break
   esac
done


year=$1
doy_start=$2
doy_end=$3

# definitions 
satname='viirs'



echo "+++++++++++VIIRS PROCESSING ++++++++++++"


# --- Loop days
for (( dd = $doy_start; dd <= doy_end; dd ++ ))
do

   doy=$(printf "%03d " ${dd} )
   month=$(date -d "01/01/${year} +${doy} days -1 day" "+%m")
   day=$(date -d "01/01/${year} +${doy} days -1 day" "+%d")

   echo $doy,$month,$day
done



