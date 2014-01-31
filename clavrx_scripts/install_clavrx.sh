#!/bin/bash
# A.Walther 7 Nov 2013
# $Header$
#!/bin/bash
#
#   history 20 January 2014: changed to patmosx branch of clavrx (AW)

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

path='clavrx_test'
if [ -n "$1" ]
then
path=$1

fi

echo
echo -e "\033[1mClavrx DEVELOPMENT branch will be installed in ===> $path \033[0m"
echo


if [ -d "$path" ]; then
   echo -e "\033[1mThis path already exists!!\033[0m"
    read -p "Do you wish to continue and update this path ?? (This may overwrite your local copy ) " yn
   case $yn in
    [Yy]* ) echo -e "\033[1mClavrx trunk branch will be installed in ===> $path \033[0m";;
    [Nn]* ) exit;;
    * ) echo "Please answer Yes or No. ";; 
  esac

fi
echo
echo '...........     cvs checkout/update programs .................'
echo

mkdir -p  $path
cd $path


svn checkout -q https://svn.ssec.wisc.edu/repos/cloud_team_clavrx/branches/development clavrx_src_development

svn checkout -q https://svn.ssec.wisc.edu/repos/cloud_team_dcomp/trunk cloud_team_dcomp
svn checkout -q https://svn.ssec.wisc.edu/repos/cloud_team_dcomp/trunk cloud_team_nlcomp



cd cloud_team_dcomp
./configure -hdf5root=$hdf5_path -with-ifort -hdflib=${hdf4_path}/lib


cd ../cloud_team_nlcomp
./configure -hdf5root=$hdf5_path -with-ifort -hdflib=${hdf4_path}/lib

cd ../clavrx_src_development/main_src
cp level2_all_on.inc level2.inc
./configure -hdf5root=$hdf5_path -with-ifort  -hdflib=${hdf4_path}/lib -hdfinc=${hdf4_path}/include -nlcomp_dir=../../cloud_team_nlcomp/ -dcomp_dir=../../cloud_team_dcomp/ -acha_dir=../cloud_acha/



if make; then 
  printf '\033[32m Clavrx development branch successfully installed %s\033[m\n'
else
   ret=$?
   printf '\033[31m Error !!!! clavrx development branch is not installed  error code $ret %s\033[m\n'
fi

#cp clavrxorb_default_options ../
#cd ../
#ln -s clavrx_bin/* ./
#./clavrxorb

exit


