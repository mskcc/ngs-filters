#!/bin/bash

SDIR="$( cd "$( dirname "$0" )" && pwd )"
SVERSION=$(git --git-dir=$SDIR/.git --work-tree=$SDIR describe --always --long)

usage() {
    echo "filterWrapper.sh FILTER_SCRIPT IN_MAF OUT_MAF [Additional parameters for Filter]"
    echo
    echo "    Run the filter (FILTER_SCRIPT) on the input maf and"
    echo "    write a properly headered output maf"
    echo
    echo "    If the filter script require more parameters they"
    echo "    can be passed after the output maf name"
    echo ""
    echo $SVERSION
    echo $SDIR
    exit
}

if [ "$#" -lt "3" ]; then
    usage
fi

FILTER=$1
MAFIN=$2
MAFOUT=$3

