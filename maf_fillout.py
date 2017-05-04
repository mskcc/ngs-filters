#!/usr/bin/env python
##########################################################################################
# MSKCC CMO
descr = 'Perform fillout operation on MAF file using GetBaseCountsMultiSample'
##########################################################################################

import argparse
import subprocess
import os
import string
import uuid
import cmo

parser = argparse.ArgumentParser(description = descr, formatter_class = argparse.RawTextHelpFormatter)
parser.add_argument('-m', '--maf', help = 'MAF file on which to filllout', required = True)
parser.add_argument('-b', '--bams', help = 'BAM files to fillout with', required = True, nargs='+')
parser.add_argument('-g', '--genome', help = 'Reference assembly of BAM files, e.g. hg19/grch37/b37', required = True, choices = cmo.util.genomes.keys())
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
#should be cmo.util.programs['GetBaseCountsMultiSample]['1.0.0'], but zeng hasn't packaged this for release yet
gbcmPath = '/home/socci/Code/Zeng/GetBaseCountsMultiSample/GetBaseCountsMultiSample'

### Set genome path
genomePath = cmo.util.genomes[args.genome]['fasta']

### Parse BAM files into string
bamString = []
for bam in bams:
	bamString.append('--bam '+os.path.splitext(os.path.basename(bam))[0]+':'+bam)
bamString = string.join(bamString)

### Check if MAF has right genome
mafGenome = subprocess.check_output('grep -v ^# '+maf+' | tail -1 | cut -f4', shell = True)
print 'Using '+genome
print 'MAF genome seems to be '+mafGenome.strip().lower()
if mafGenome.strip().lower() is not genome.lower():
	print 'Genome build different from that in MAF file, might fail'

### Check if genome in BAM header 
for bam in bams:
	try:
		out = subprocess.check_output('samtools view -H '+bam+' | grep '+genome, shell = True)
	except subprocess.CalledProcessError:
		print 'Genome in '+bam+' does not agree with input genome'

### Make a temporary simplified MAF
tmpMaf = uuid.uuid4().hex+'_tmp.maf'
uniqRscript = os.path.dirname(os.path.realpath(__file__))+'/maf_uniq_tags.R'
uniqRCall = uniqRscript+' '+maf+' > '+tmpMaf
subprocess.call(uniqRCall, shell = True)

### Call GetBaseCountsMultiSample
gbcmCall = gbcmPath+' --thread %s --filter_improper_pair 0 --fasta %s --maf %s --output %s %s' % (n, genomePath, tmpMaf, output, bamString)
subprocess.call(gbcmCall, shell = True)

### Remove temporary MAF
subprocess.call('rm -f '+tmpMaf, shell = True)
