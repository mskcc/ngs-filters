#!/bin/bash

if [ "$#" -lt "2" ]; then
    echo "usage: do_Normal_Cohort_Test.sh INPUT.maf NC_FILL.maf"
    exit 1
fi

INPUT=$1
FILL=$2
shift 2

TAG=_$(uuidgen)
#echo TAG=$TAG

../filter_cohort_normals.R  -m $INPUT -f $FILL -o _${TAG}_OUTPUT.maf $*
egrep -v "^#" _${TAG}_OUTPUT.maf | cut -f120 | sort | uniq -c > _${TAG}_TARGET.out

diff _Normal_Cohort_Test__TARGET _${TAG}_TARGET.out >_${TAG}_DIFF
nDiff=$(wc -l _${TAG}_DIFF | awk '{print $1}')

if [ "$nDiff" != "0" ]; then
    echo TEST $(basename $0) FAILED
    echo
    cat _${TAG}_DIFF
    exit 1
else
    echo TEST $(basename $0) PASSED
fi