#!/bin/bash
# $Id: get_gfs.sh,v 1.4 2012/10/22 16:18:09 awalther Exp $

args=("$@") 
year=${args[0]} 
doy=${args[1]}


#gfs_path='/data/Ancil_Data/gfs/'$year'/'
gfs_path='/data/Ancil_Data/gfs/'
 [ ! -d $gfs_path ] && mkdir -p -v $gfs_path
 
 if date -v 1d > /dev/null 2>&1; then
  doy_dum=`expr ${doy} - 1` 
 
  DAY_OBS=$(./shift_date.sh ${year}0101 +${doy_dum})
  DAY_OBS=${DAY_OBS:2}
else
  
  DAY_OBS=$(date -d "01/01/${year} +${doy} days -1 day" "+%y%m%d")
  echo $(date -d "01/01/${year} +${doy} days -1 day" "+%y%m%d")
fi
#eval $DAY_OBS

 
echo 'check and/or get GFS data from '$DAY_OBS

for hhh in 00 06 12 18
do
  filebase='gfs.'$DAY_OBS$hhh'_F012.hdf'
  file=$gfs_path$filebase
  file_bz=$file.bz2
  
  if [ ! -f $file ]; 
  then
   sh -c 'scp -p -r thor:/data1/Ancil_Data/gfs/hdf/gfs.'$DAY_OBS$hhh'_F012.hdf* '$gfs_path
   [  -f $file_bz ] &&  sh -c 'bunzip2 -v '$file'.bz2' || echo 'all clear'
   
  fi
done

if [ ! -f $file ];
then
  echo " no NWP data are available!!"
  trap 9
  
fi

 if date -v 1d > /dev/null 2>&1; then
  doy_dum=`expr ${doy} - 2` 

  DAY_OBS=$(./shift_date.sh ${year}0101 +${doy_dum})
  DAY_OBS=${DAY_OBS:2}
else
  
  DAY_OBS=$(date -d "01/01/${year} +${doy} days -2 day" "+%y%m%d")
  echo $(date -d "01/01/${year} +${doy} days -2 day" "+%y%m%d")
fi


for hhh in 00 06 12 18
do
  filebase='gfs.'$DAY_OBS$hhh'_F012.hdf'
  file=$gfs_path$filebase
  file_bz=$file.bz2
  
  if [ ! -f $file ]; 
  then
   sh -c 'scp -p -r thor:/data1/Ancil_Data/gfs/hdf/gfs.'$DAY_OBS$hhh'_F012.hdf* '$gfs_path
   [  -f $file_bz ] &&  sh -c 'bunzip2 -v '$file'.bz2' || echo 'all clear'
   
  fi
done
