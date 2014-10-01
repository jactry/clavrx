#!/usr/bin/env bash
# $Id:$


function usage() {

cat<< EOF

cg_run_viirs
need help? bad luck,  sorry!
ask denis.botambekov@ssec.wisc.edu 
or andi.walther@ssec.wisc.edu 

EOF

}

satname="viirs"

START=$1
END=$2




YEAR0=${START:0:4}
DOY0=${START:5:3}


YEAR1=${END:0:4}
DOY1=${END:5:3}




PWD_LOC=${PWD}
echo $PWD_LOC
for (( ddd = $DOY0; ddd < $DOY1; ddd++ ))
do
	l1b_path=$HOME'/Satellite_Input/'$satname'/'$YEAR0'/'$ddd'/'
	mkdir -p $l1b_path
	out_path=$HOME'/Satellite_Output/'$satname'/'$YEAR0'/'$ddd'/'
	mkdir -p $out_path
	./get_gfs.sh $YEAR0 $ddd
	for (( hhh = 0; hhh < 24; hhh++ ))
	do
		cp cg_peate_downloader.sh $l1b_path
      cd $l1b_path
      echo $l1b_path
		./cg_get_viirs_data.sh --reg spc 20140312
      cd $PWD_LOC
      echo ${PWD}
		./write_filelist.sh $l1b_path $out_path GMTC
		#./clavrxorb -dcmp_mode $dcomp_mode
	done

done





usage
exit
