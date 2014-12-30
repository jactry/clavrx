#!/bin/sh
# $Id: run_viirs.sh,v 1.4 2012/08/02 07:01:51 awalther Exp $
# !!!!! TO RUN ON THE CLUSTER !!!!!!!!!
# qsub -q all.q -S /bin/bash -l matlab=0 -l friendly=1 -p -200 run_viirs.sh
# qsub -q r720.q -l vf=4G -S /bin/bash -l matlab=0 -l friendly=1 -p -200 run_viirs.sh

#args=("$@")

# --- Load hdf5 & intel12 libs
source /etc/bashrc
module load bundle/basic-1
module load bundle/basic-1 hdf5

# --- Set year, date, etc.
year=2014 #${args[0]} 
doy_start=1 #3,99,192,293 #${args[1]} 
doy_end=1 #$doy_start #336 #152 #${args[2]} 
hour0=0
hour1=0 #$hour0 #23  #$hour0
satname="VIIRS"
#region='global'
#region='sahara'
#region='dom_c'
#region='north_pacific'
region='tmp'
work_dir='/fjord/jgs/personal/dbotambekov/patmosx_processing/scripts/'

echo "+++++++++++VIIRS PROCESSING ++++++++++++"

zara_files_path='/home/dbotambekov/src_clavrx_code/zara_run/'
logs_path='/fjord/jgs/personal/dbotambekov/logs/'
options='clavrxorb_default_options_viirs_zara'
filelist='clavrxorb_file_list'
filetype='GMTCO'

# --- Loop days
for (( doy = $doy_start; doy <= doy_end; doy ++ ))
do

  # --- check if doy is less than 10 or 100
  if [ $doy -lt 10 ] ; then
    doy_str=`expr 00$doy`
  elif [ $doy -ge 10 ] && [ $doy -lt 100 ] ; then
    doy_str=`expr 0$doy`
  else
    doy_str=$doy
  fi
  echo 'Pocessing DOY = '$doy_str

  month=$(date -d "01/01/${year} +${doy} days -1 day" "+%m")
  day=$(date -d "01/01/${year} +${doy} days -1 day" "+%d")

# --- Loop through the hours
  for (( hhh = $hour0; hhh <= $hour1; hhh ++ ))
  do
    hhh_str=S$hhh
    hhh_str1=$hhh
    if [ $hhh -lt 10 ] ; then
      hhh_str=S`expr 0$hhh`
      hhh_str1=`expr 0$hhh`
    fi

  l1b_path='/fjord/jgs/personal/dbotambekov/Satellite_Input/'$satname'/'$region'/'$year'/'$doy_str'/'
  out_path='/fjord/jgs/personal/dbotambekov/Satellite_Output/'$satname'/'$region'/'$year'/'$doy_str'/'

   # !!!!!!!!! CREATE A NEW TEMP SCRIPT TO SUBMIT IT TO ZARA
   tmp_script=$work_dir'npp_'$year'_'$doy_str'_'$hhh_str1'_'$region'_patmosx.sh'
   tmp_work_dir=$work_dir'npp_'$year'_'$doy_str'_'$hhh_str1'_'$region
   echo "#!/bin/sh" > $tmp_script
   echo "source /etc/bashrc" >> $tmp_script
   echo "module load bundle/basic-1" >> $tmp_script
   echo "module load bundle/basic-1 hdf5" >> $tmp_script

   echo "echo 'Processing viirs $year $doy_str $hhh_str1 $region'" >> $tmp_script

   echo "echo 'Creating necessary dirs and copying files'" >> $tmp_script
   echo "[ ! -d $tmp_work_dir ] && mkdir -v -p $tmp_work_dir" >> $tmp_script
   echo "[ ! -d $tmp_work_dir/temporary_files ] && mkdir -v -p $tmp_work_dir/temporary_files" >> $tmp_script
   echo "cp $zara_files_path/peate_downloader.sh $tmp_work_dir" >> $tmp_script
   echo "cp $zara_files_path/get_viirs_data_zara.sh $tmp_work_dir" >> $tmp_script
   echo "cp $zara_files_path/sync_viirs_zara.sh $tmp_work_dir" >> $tmp_script
   echo "cp $zara_files_path/write_filelist_zara.sh $tmp_work_dir" >> $tmp_script
   echo "cp $zara_files_path/check_filelist_zara.sh $tmp_work_dir" >> $tmp_script
   echo "cp $zara_files_path/clavrxorb_trunk $tmp_work_dir" >> $tmp_script
   echo "cp $zara_files_path/comp_asc_des_level2b $tmp_work_dir" >> $tmp_script
   echo "cp $zara_files_path/$options $tmp_work_dir" >> $tmp_script
   echo "[ ! -d $l1b_path ] && mkdir -v -p $l1b_path" >> $tmp_script
   echo "[ ! -d $out_path ] && mkdir -v -p $out_path" >> $tmp_script

   echo "echo 'Linking Ancil Data'" >> $tmp_script
   echo "[ ! -d $tmp_work_dir/clavrx_ancil_data ] && mkdir -v -p $tmp_work_dir/clavrx_ancil_data" >> $tmp_script
   echo "ln -s /fjord/jgs/patmosx/Ancil_Data/clavrx_ancil_data/avhrr_data $tmp_work_dir/clavrx_ancil_data" >> $tmp_script
   echo "ln -s /fjord/jgs/patmosx/Ancil_Data/clavrx_ancil_data/bayes $tmp_work_dir/clavrx_ancil_data" >> $tmp_script
   echo "ln -s /fjord/jgs/patmosx/Ancil_Data/clavrx_ancil_data/insolation $tmp_work_dir/clavrx_ancil_data" >> $tmp_script
   echo "ln -s /fjord/jgs/patmosx/Ancil_Data/clavrx_ancil_data/luts $tmp_work_dir/clavrx_ancil_data" >> $tmp_script
   echo "ln -s /fjord/jgs/patmosx/Ancil_Data/clavrx_ancil_data/sfc_data $tmp_work_dir/clavrx_ancil_data" >> $tmp_script
   echo "ln -s /fjord/jgs/patmosx/Ancil_Data/clavrx_ancil_data/dcomp_config $tmp_work_dir/clavrx_ancil_data" >> $tmp_script
   echo "cp -r /fjord/jgs/patmosx/Ancil_Data/clavrx_ancil_data/pfast/ $tmp_work_dir/clavrx_ancil_data" >> $tmp_script

   echo "echo 'Getting l1b data'" >> $tmp_script
   echo "cd $tmp_work_dir" >> $tmp_script
   echo "cp peate_downloader.sh $l1b_path" >> $tmp_script
   echo "./get_viirs_data_zara.sh $year $doy $hhh $l1b_path" >> $tmp_script
   echo "echo 'Making sure all files are there, running sync_viirs_zara.sh'" >> $tmp_script
   echo "./sync_viirs_zara.sh $l1b_path" >> $tmp_script

   echo "echo 'Writing files to the filelist'" >> $tmp_script
   echo "./write_filelist_zara.sh $l1b_path $out_path $filetype d$year$month$day t$hhh_str1" >> $tmp_script
   echo "echo 'Checking files, if already processed delete them from the filelist'" >> $tmp_script
   echo "./check_filelist_zara.sh $filelist $filetype" >> $tmp_script
   echo "echo 'Starting CLAVR-x'" >> $tmp_script
   echo "./clavrxorb_trunk  -default $options -lines_per_seg 400" >> $tmp_script
   echo "echo 'Finished, Deleting All Temp Data'" >> $tmp_script
#   echo "rm -rf $work_dir" >> $tmp_script
#   echo "rm -rf $l1b_path" >> $tmp_script
#   echo "rm $tmp_script" >> $tmp_script

   # --- Submit job to zara
   qsub -q r720.q -l vf=4G -S /bin/bash -l matlab=0 -l friendly=1 -p -200 -o $logs_path -e $logs_path -l h_rt=06:00:00 $tmp_script

   
  done # hours loop

done # days loop

