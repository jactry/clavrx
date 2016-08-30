Dark Composite Code for GOES Area Files

last upadte 20141016 C.C.Molling

This code computes the second darkest pixel at each pixel location of the reference 
channel 1 area file.  This second darkest pixel is used by the SASRAB algorithm to compute
surface total insolation and diffuse insolation.  This version allows you to have all the 
ch1 area files in either one directory or yyyy_doy directories.  These files can be either
uncompresssed or compressed (.gz or .bz2).  The output can be the historical .dat binary 
format or area file format.  Output will be either uncompressed or compressed (gzip).

Your executable clavrxorb may not read in area format at this time, so check to make sure
which output format you need!!!  Versions before October 16, 2014 definitely need the .dat 
format.

1) Look inside make_gsip_dark_ch1 to see if you need to change compiler 
or flags.  Edit as needed.

2) Make the program

make

which will create 'dark_ch1'.  Eventually, this step will be done with the install script.

See dark_ch1_start_file.example for the input file format.

3) At this time the program, 'batch_dark_ch1', a batch script to run the program on
multiple files, is not made because it's code does not match the input for dark_ch1.
If you do correct the code batch_dark_ch1.f90, you should change the makefile so it
is compiled with "make".

Create an input file called batch_dark_ch1_input in order to use batch_dark_ch1.

4) At this time, you need to edit the 
1) Look inside make_gsip_dark_ch1 to see if you need to change compiler 
or flags.  Edit as needed.

2) Make the program

make

which will create 'dark_ch1'.  Eventually, this step will be done with the install script.

See dark_ch1_start_file.example for the input file format.

3) At this time the program, 'batch_dark_ch1', a batch script to run the program on
multiple files, is not made because it's code does not match the input for dark_ch1.
If you do correct the code batch_dark_ch1.f90, you should change the makefile so it
is compiled with "make".

Create an input file called batch_dark_ch1_input in order to use batch_dark_ch1.

4) At this time, you need to edit the 
1) Look inside make_gsip_dark_ch1 to see if you need to change compiler 
or flags.  Edit as needed.

2) Make the program

make

which will create 'dark_ch1'.  Eventually, this step will be done with the install script.

See dark_ch1_start_file.example for the input file format.

3) At this time the program, 'batch_dark_ch1', a batch script to run the program on
multiple files, is not made because it's code does not match the input for dark_ch1.
If you do correct the code batch_dark_ch1.f90, you should change the makefile so it
is compiled with "make".

Create an input file called batch_dark_ch1_input in order to use batch_dark_ch1.

4) At this time, you need to edit the file dark_ch1_gsip.f90 and create new domain
logic each time you have a new area coverage (e.g., you make a florida coverage).
Then recompile.
