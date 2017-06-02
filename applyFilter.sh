#!/bin/bash

SDIR="$( cd "$( dirname "$0" )" && pwd )"
SVERSION=$(git --git-dir=$SDIR/.git --work-tree=$SDIR describe --always --long)
VTAG="$(basename $SDIR)/$(basename $0) VERSION=$SVERSION"

usage() {
    echo "applyFilter.sh FILTER_SCRIPT IN_MAF OUT_MAF [Additional parameters for Filter]"
    echo
    echo "    Run the filter (FILTER_SCRIPT) on the input maf and"
    echo "    write a properly headered output maf"
    echo
    echo "    If the filter script require more parameters they"
    echo "    can be passed after the output maf name"
    echo ""
    echo "  "$VTAG
    echo
    exit
}

if [ "$#" -lt "3" ]; then
    usage
fi

FILTER=$1
MAFIN=$2
MAFOUT=$3

shift 3

# Resolve filter to full path
FILTER=$SDIR/$(basename $FILTER)
if [ ! -e $FILTER ]; then
    echo "ERROR: Non-existant filter" $(basename $FILTER)
    exit 1
fi

#
# Check if MAF has proper version header
# If not add it.
#

HEADER=$(head -1 $MAFIN)
if [[ ! "$HEADER" =~ ^# ]]; then
    echo "#version 2.4" > $MAFOUT
else
    egrep "^#" $MAFIN > $MAFOUT
fi

# Add version tag

echo "#$VTAG FILTER=$(basename $FILTER)" >>$MAFOUT

TMPMAF=$(uuidgen).maf


Rscript --vanilla  $FILTER -m $MAFIN -o $TMPMAF $*
EXIT_CODE=$?

if [ "$EXIT_CODE" != "0" ]; then
    echo "ERROR IN R-script"
    exit 1
fi


cat $TMPMAF >>$MAFOUT
rm $TMPMAF
