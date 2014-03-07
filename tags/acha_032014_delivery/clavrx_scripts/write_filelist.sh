#!/bin/sh
# $Id: write_filelist.sh,v 1.4 2012/07/03 22:53:54 awalther Exp $

args=("$@") 
dir_1b_local=${args[0]} 
output_path_local=${args[1]} 




echo 'local level1b path: '$dir_1b_local
echo
mkdir -p output_path_local

dir_1bx=$dir_1b_local
dir_nav_in=$dir_1b_local
dir_nav_out=$output_path_local
dir_cmr=$output_path_local
dir_sst=$output_path_local
dir_cld=$output_path_local
dir_obs=$output_path_local
dir_geo=$output_path_local
dir_rtm=$output_path_local
dir_ash=$output_path_local
dir_level2=$output_path_local
dir_level2b_daily=$output_path_local
dir_level3=$output_path_local
dir_level3_daily=$output_path_local

echo $dir_1b_local > clavrxorb_file_list
echo $dir_1bx >> clavrxorb_file_list
echo $dir_nav_in >> clavrxorb_file_list
echo $dir_nav_out >> clavrxorb_file_list
echo $dir_cmr >> clavrxorb_file_list
echo $dir_sst >> clavrxorb_file_list
echo $dir_cld >> clavrxorb_file_list
echo $dir_obs >> clavrxorb_file_list
echo $dir_geo >> clavrxorb_file_list
echo $dir_rtm >> clavrxorb_file_list
echo $dir_ash >> clavrxorb_file_list
echo $dir_level2 >> clavrxorb_file_list
echo $dir_level3 >> clavrxorb_file_list

case $# in
3)ls $dir_1b_local |  grep $3  >> clavrxorb_file_list;;
4)ls $dir_1b_local |  grep $3 | grep $4  >> clavrxorb_file_list;;
5)ls $dir_1b_local |  grep $3 | grep $4 | grep $5  >> clavrxorb_file_list;;
6)ls $dir_1b_local |  grep $3 | grep $4 | grep $5  >> clavrxorb_file_list;;
esac

echo
echo "files to be processed========>  ",$3,$4,$5
case $# in
3)ls $dir_1b_local |  grep $3 ;;
4)ls $dir_1b_local |  grep $3 | grep $4  ;;
5)ls $dir_1b_local |  grep $3 | grep $4 | grep $5 ;;
esac
echo
