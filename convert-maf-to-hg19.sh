#!/bin/bash

if [ "$1" == "-h" -o $# -eq 0 ]; then
  echo "Usage: `basename $0` /path/to/file.maf"
  exit 0
fi

MAF=$1
OUTNAME=$(basename ${MAF%.*})

grep -v "#" $MAF | head -1 > ${OUTNAME}_hg19.maf
grep -v "#" $MAF | tail -n+2 | awk 'BEGIN{OFS = FS = "\t"} {$5 = "chr"$5; print}' >> ${OUTNAME}_hg19.maf
sed -i 's/chrMT/chrM/g' ${OUTNAME}_hg19.maf
