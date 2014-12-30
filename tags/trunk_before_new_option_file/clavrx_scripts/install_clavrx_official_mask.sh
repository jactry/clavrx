#!/bin/bash
# A.Walther 7 Nov 2013
# $Header: https://svn.ssec.wisc.edu/repos/cloud_team_clavrx/trunk/clavrx_scripts/install_clavrx_trunk.sh 36 2014-02-05 22:20:03Z heidinger $
#!/bin/bash
#
#   history 20 January 2014: changed to patmosx branch of clavrx (AW)
#            24 March 2014: for CLAVR-x stable version 

#source ~/.bashrc
#export DISPLAY=:1

set -e

hdf5_path="/opt/hdf5-1.8.8-intel"
hdf4_path="/usr/local/hdf4"

# for convinience: symbol link to local directory
# I encourage everybody to do this on all machines
# ln -s <hdf4_path> ~/lib/hdf4
# ln -s <hdf5_path> ~/lib/hdf5
hdf5_path=$HOME"/lib/hdf5/"
hdf4_path=$HOME"/lib/hdf4/"

stable_version='official_cloud_mask'

path=$stable_version
if [ -n "$1" ]
then
path=$1

fi

echo
echo -e "\033[1mClavrx official cloud mask will be installed in ===> $path \033[0m"
echo


if [ -d "$path" ]; then
   echo -e "\033[1mThis path already exists!!\033[0m"
    read -p "Do you wish to continue and update this path ?? (This may overwrite your local copy ) " yn
   case $yn in
    [Yy]* ) echo -e "\033[1mClavrx official cloud mask branch will be installed in ===> $path \033[0m";;
    [Nn]* ) exit;;
    * ) echo "Please answer Yes or No. ";; 
  esac

fi
echo
echo '...........     svn checkout clavrx official cloud mask version 5.4 .................'
echo

mkdir -p  $path
cd $path

svn checkout -q https://svn.ssec.wisc.edu/repos/cloud_team_clavrx/branches/official_cloud_mask ./


cd dcomp
./configure -hdf5root=$hdf5_path -with-ifort -hdflib=${hdf4_path}/lib


cd ../nlcomp
./configure -hdf5root=$hdf5_path -with-ifort -hdflib=${hdf4_path}/lib

cd ../main_src
cp level2_all_on.inc level2.inc
./configure -hdf5root=$hdf5_path -with-ifort  -hdflib=${hdf4_path}/lib -hdfinc=${hdf4_path}/include -nlcomp_dir=../nlcomp/ -dcomp_dir=../dcomp/ -acha_dir=../cloud_acha/



if make; then 
  printf '\033[32m Clavrx official cloud mask successfully installed %s\033[m\n'
else
   ret=$?
   printf '\033[31m Error !!!! clavrx official cloud mask is not installed  error code $ret %s\033[m\n'
fi

#cp clavrxorb_default_options ../
#cd ../
#ln -s clavrx_bin/* ./
#./clavrxorb

exit


