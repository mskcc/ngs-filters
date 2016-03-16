#!/usr/bin/env python

##########################################################################################
# MSKCC CMO
descr = 'Perform fillout operation on MAF file using GetBaseCountsMultiSample'
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
parser.add_argument('-n', '--n_threads', help = 'Multithread', default = 10, required = False)
args = parser.parse_args()

maf = args.maf
bams = args.bams
genome = args.genome.lower()
n = args.n_threads
if args.output is None:
	output = os.path.splitext(os.path.basename(maf))[0]+'.fillout'
else:
	output = args.output

### Path to GetBaseCountsMultiSample
gbcmPath = '/home/socci/Code/Zeng/GetBaseCountsMultiSample/GetBaseCountsMultiSample'

### Set genome path
if genome == 'hg19':
	genomePath = '/common/data/assemblies/H.sapiens/hg19/hg19.fasta'
if genome == 'grch37':
	genomePath = '/ifs/depot/assemblies/H.sapiens/GRCh37/gr37.fasta'
if genome == 'b37':
	genomePath = '/ifs/depot/assemblies/H.sapiens/b37/b37.fasta'
if genome == 'b37_dmp':
	genomePath = '/ifs/depot/assemblies/H.sapiens/b37_dmp/b37.fasta'

### Parse BAM files into string
bamString = []
for bam in bams:
	bamString.append('--bam '+os.path.splitext(os.path.basename(bam))[0]+':'+bam)
bamString = string.join(bamString)

### Check if MAF has right genome
mafGenome = subprocess.check_output('grep -v ^# '+maf+' | tail -1 | cut -f4', shell = True)
print 'Using '+genome
print 'MAF genome seems to be '+mafGenome.strip().lower()
if mafGenome.strip().lower() is not genome:
	print 'Genome build different from that in MAF file, might fail'

### Check if genome in BAM header 
for bam in bams:
	try:
		out = subprocess.check_output('samtools view -H '+bam+' | grep '+genome, shell = True)
	except subprocess.CalledProcessError:
		print 'Genome in '+bam+' does not agree with input genome'

### Make a temporary simplified MAF
uniqRscript = os.path.dirname(os.path.realpath(__file__))+'/maf_uniq_tags.R'
uniqRCall = uniqRscript+' '+maf+' > ___tmp.maf'
subprocess.call(uniqRCall, shell = True)

### Call GetBaseCountsMultiSample
gbcmCall = gbcmPath+' --thread %s --filter_improper_pair 0 --fasta %s --maf ___tmp.maf --output %s %s' % (n, genomePath, output, bamString)
subprocess.call(gbcmCall, shell = True)

### Remove temporary MAF
subprocess.call('rm -f ___tmp.maf', shell = True)
