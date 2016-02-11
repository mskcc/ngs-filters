#!/usr/bin/env python

##########################################################################################
##########################################################################################
# MSKCC CMO
descr = 'Perform fillout operation on MAF file using GetBaseCountsMultiSample'
##########################################################################################
##########################################################################################

import argparse
import subprocess
import os
import string

parser = argparse.ArgumentParser(description = descr, formatter_class = argparse.RawTextHelpFormatter)
parser.add_argument('-m', '--maf', help = 'MAF file on which to filllout', required = True)
parser.add_argument('-b', '--bams', help = 'BAM files to fillout with', required = True, nargs='+')
parser.add_argument('-g', '--genome', help = 'Reference assembly of BAM files, e.g. hg19/grch37/b37', required = True)
parser.add_argument('-o', '--output', help = 'Prefix for output file', required = False)
args = parser.parse_args()

maf = args.maf
bams = args.bams
genome = args.genome
if args.output is None:
	output = os.path.splitext(os.path.basename(maf))[0]+'.fillout'
else:
	output = args.output

### Path to GetBaseCountsMultiSample
gbcmPath = '/home/socci/Code/Zeng/GetBaseCountsMultiSample/GetBaseCountsMultiSample'

### Set genome path
if genome == 'hg19':
	genomePath = '/ifs/depot/assemblies/H.sapiens/hg19/hg19.fasta'
if genome == 'grch37':
	genomePath = '/ifs/depot/assemblies/H.sapiens/GRCh37/grch37.fasta'
if genome == 'b37':
	genomePath = '/ifs/depot/assemblies/H.sapiens/b37/b37.fasta'

### Parse BAM files into string
bamString = []
for bam in bams:
	bamString.append('--bam '+os.path.splitext(os.path.basename(bam))[0]+':'+bam)
bamString = string.join(bamString)

### Call GetBaseCountsMultiSample
gbcmCall = gbcmPath+' --thread 20 --filter_improper_pair 0 --fasta %s --maf %s --output %s %s' % (genomePath, maf, output, bamString)
subprocess.call(gbcmCall, shell = True)