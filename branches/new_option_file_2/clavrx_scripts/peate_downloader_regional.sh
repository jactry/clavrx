#!/usr/bin/env bash

function usage () {
  cat <<EOF

Atmosphere PEATE Downloader Script
----------------------------------
Usage: $0 [-t (curl|wget|txt)] <start> <end> <file_type> [<file_type> ...]

  -t          Output type, (default cur) If curl or wget, the downloaded 
              script will be excecuted after download.
  -x          Use xargs to download 4 files sumultaneously (curl, wget only).
  start, end  Query window. YYYY-MM-DD[+HH:MM:SS]
  file_type   Ingest System file type name (RCRIS-RNSCA, SVI01, etc...)

Need help?  peate.admin@ssec.wisc.edu

EOF
}

XARGS="0"
OUTPUT=curl
while getopts xt: o
do 
case "$o" in
  t) OUTPUT=$OPTARG;;
  x) XARGS="1";;
  *)
esac
done
shift `expr $OPTIND - 1`
LON=$1
LAT=$2
START=$3
END=$4
shift 4
TYPES=$*

echo "lldlld"$LON$LAT

# Make sure start, end, and at least one file type are provided
if [[ -z $START ]] || [[ -z $END ]] || [[ -z $TYPES ]]; then
  usage
  exit 1
fi

URL="http://peate.ssec.wisc.edu/flo/api/find?start=$START&end=$END&output=$OUTPUT&xargs=$XARGS"
for ft in $TYPES; do
  URL="$URL&file_type=$ft"
 
done

 URL="$URL&loc=$LAT,$LON"
  URL="$URL&radius=2000"
echo "$URL"

SCRIPT=downloader.sh
if [ $OUTPUT = "txt" ]; then
  SCRIPT=downloader.txt
fi

wget -q -O $SCRIPT ${URL}

if [ "$?" -ne "0" ]; then
  echo "Retreiving file list failed"
  rm $SCRIPT &> /dev/null
  exit 1
else
  echo "Downloaded $SCRIPT"
fi

if [ $OUTPUT = "txt" ]; then
  exit 0
fi

echo "Running $SCRIPT"
bash $SCRIPT

if [ ! $? ]; then
  echo "Warning: $SCRIPT returned non-zero status"
fi

