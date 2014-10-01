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
		sh -c './cg_peate_downloader.sh '$YEAR0'-'$MONTH0'-'$DAY0'+'$hhh':00:00 '$YEAR0'-'$MONTH0'-'$DAY0'+'$hhh':00:00 SVDNB GMTCO GITCO GDNBO SVM01 SVM02 SVM03 SVM04 SVM05 SVM06 SVM07 SVM08 SVM09 SVM10 SVM11 SVM12 SVM13 SVM14 SVM15 SVM16 IICMO SVI01'
      cd $PWD_LOC
      echo ${PWD}
		./write_filelist.sh $l1b_path $out_path GMTC
		#./clavrxorb -dcmp_mode $dcomp_mode
	done

done





usage
exit
