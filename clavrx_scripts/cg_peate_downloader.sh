#!/usr/bin/env bash

function usage() {

cat<< EOF


need help? sorry!

EOF

}

START=$1
END=$2
shift 2
TYPES=$*

OUTPUT=wget
URL="http://peate.ssec.wisc.edu/flo/api/find?start=$START&end=$END&output=$OUTPUT"

for ft in $TYPES; do
	URL="$URL&file_type=$ft"
done

echo $URL

SCRIPT=downloader.sh

wget -q -O $SCRIPT ${URL}

usage
exit 1
