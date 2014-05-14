How to install and run PATMOS-x clavrxorb
last update 20130719 C.C. Molling

Instructions valid for...
akh_clavrx_src distribution branch akh_clavrx_src_development



***Note*** If you have previously been using the version of the Makefile
that puts the executables in one directory up from akh_clavrx_src, 
please remove the old executables before compiling this code.  New 
executables will go into ../clavrx_bin/.
-------------------------------------------------------------

Some of these instructions are a repeat of information found at
 https://groups.ssec.wisc.edu/groups/clavr-x/clavr-x-delivery-documentation/installing-clavr-x-on-your-own-machine

Clavrxorb
=========

You will need to first check out from the CVS server the clavrx code, 
the ACHA code, and the DCOMP code.  You may also need grib2hdf if you will
be converting your own NWP files.

>cd my_patmosx_src_dir
>cvs checkout -r akh_clavrx_src_development akh_clavrx_src
>cvs checkout cloud_team_acha
>cvs checkout cloud_team_dcomp
>cvs checkout cloud_team_nlcomp
>cvs checkout akh_cloud_type_repo
>cvs checkout jpssrr_nbvcm
>cvs checkout grib2hdf

Before compiling, you use the configure command to create the file config.mk
which contains all the relevant paths and flags needed by the Makefile.

1) Enter the patmosx src directory

>cd my_patmosx_src_dir/akh_clavrx_src

2) Configure the code.  You must specify the compiler, the hdf library and 
include file locations (if not /usr/local/hdf/lib/ and /usr/local/hdf/include/), 
the acha src directory (if not ../cloud_team_acha/), the dcomp src directory
(if not ../cloud_team_dcomp/) and the cloud type directory (if not
../akh_cloud_type_repo/).  If you wish to run clavrxorb with VIIRS, you
must specify the -hdf5root directory.  Omit -hdf5root in the configure command 
if you will not be running on VIIRS data.

>./configure --with-ifort -hdflib=(HDF4 lib location)
   -hdfinc=(HDF4 include files location) -acha_dir=(ACHA directory)
   -hdf5root=(dif containing HDF5 lib and include)

Above, the location of your HDF4 library is indicated with -hdflib
and the location of your HDF5 root directory is indicated with -hdf5root.
Note that configure currently requires you to have hdf5's lib and
include as subdirectories of the same hdf5root directory (if used).

Example: configure the code to use the ifort compiler...

>./configure --with-ifort -hdflib=/home/heidinger/hdf4.2-gcc/lib
   -hdfinc=/opt/hdf4-intel/include/ -hdf5root=/opt/hdf5-1.8.8-intel/

WARNING!  At this time the configure process does not configure the code
inside aw_dcomp, so you must do this separately.  In cloud_team_dcomp/ issue the command

> ./configure --help

to get the help page for the configure.  Then, similarly to the configure 
above, issue the configure command with the correct options and library 
locations.


2) Before compiling, in akh_clavrx_src copy over a level2*.inc to level2.inc, 
and edit this level2.inc to turn on and off desired output variables.
Then compile using the command...

> make

The make process will compile everything, create a directory called ../clavrx_bin/ and
place all the executables in ../clavrx_bin/.

3) To run the clavrxorb code...

Below, (srcdir) is the directory in which this source code exits
(the full or relative path for akh_clavrx_src)
and (workingdir) is a scratch directory where you can run the
program and create diagnostic output.  (outputdir) is where you want 
the permanent output to be written.

>cd (workingdir)
>mkdir temporary_files diag_output (outputdir)
>cp (srcdir)/clavrxorb_default_options my_clavrxorb_default_options

Edit my_clavrxorb_default_options to reflect stuff you want and 
where input files are.

Then either copy to (workingdir) the file clavrxorb_file_list and edit, or use a
script like run1day.bash or patmosx_goes_daily.pl (in 
patmosx_tools CVS dir).

If not using a script, use the command...

>(srcdir)/../clavrx_bin/clavrxorb -filelist clavrxorb_file_list -default
my_clavrxorb_default_options



Grib2hdf
========

Please follow the instructions found at
 https://groups.ssec.wisc.edu/groups/clavr-x/clavr-x-delivery-documentation/installing-clavr-x-on-your-own-machine
