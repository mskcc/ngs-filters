#!/usr/bin/env python

##########################################################################################
# MSKCC CMO
descr = 'Wrapper for fpfilter.pl\.Technical filtering of SNV calls from NGS experiments.\
		\nRequires bam-readcount, bgzip, tabix and vcf-annotate in PATH.\
		\nfpfilter.pl available at https://github.com/ckandoth/variant-filter.\
		\nfpfilter.pl criteria and bam-readcount arguments are subject to changes depending\
		on experiment type.'
##########################################################################################

import argparse
import subprocess
import os
import string

parser = argparse.ArgumentParser(description = descr, formatter_class = argparse.RawTextHelpFormatter)
parser.add_argument('-v', '--vcf', help = 'VCF file with SNVs', required = True)
parser.add_argument('-b', '--bam', help = 'Tumor sample BAM file', required = True)
parser.add_argument('-g', '--genome', help = 'Reference assembly of BAM files, e.g. hg19/grch37/b37', required = True, choices=cmo.util.genomes.keys())
parser.add_argument('-f', '--fppath', help = 'Full path to fpfilter.pl', required = False)
args = parser.parse_args()

### Set up environment
vcf = args.vcf
bam = args.bam
genome = args.genome
prefix = os.path.splitext(os.path.basename(vcf))[0]
if args.fppath is None:
	fpfilter = 'fpfilter.pl'
else:
	fpfilter = args.fppath

genomePath = cmo.util.genomes[args.genome]['fasta']

outPath = prefix+'_fpfilter'
if not os.path.exists(outPath):
	os.makedirs(outPath)

### Generate list of SNVs
print('Making list of SNVs')
varOut = outPath+'/'+prefix+'.var'
varCmd = "perl -ane 'print join(\"\t\",@F[0,1,1]).\"\n\" unless(m/^#/)' %s  > %s" % (vcf, varOut)
subprocess.call(varCmd, shell = True)

### Run bam-readcount on on variants in sample BAM
print('Counting reads')
readcountOut = outPath+'/'+prefix+'.readcount'
readcountCmd = 'bam-readcount -q1 -b15 -w1 -f %s -l %s %s > %s' % (genomePath, varOut, bam, readcountOut)
subprocess.call(readcountCmd, shell = True)

### Run fpfilter.pl
print('Running fpfilter.pl')
fpfilterOut = outPath+'/'+prefix+'.fpfilter'
fpfilterCmd = 'perl %s --var-file %s --readcount-file %s --output-file %s' % (fpfilter, vcf, readcountOut, fpfilterOut) 
subprocess.call(fpfilterCmd, shell = True)

### Prepare for vcf-annotate
print('Indexing filter output')
tabixOut = fpfilterOut+'.gz'
tabixCmd = 'bgzip -c %s > %s && tabix -p vcf %s' % (fpfilterOut, tabixOut, tabixOut)
subprocess.call(tabixCmd, shell = True)

### Annotate VCF 
print('Annotating VCF')
annotateOut = prefix+'_filter.vcf'
annotateCmd = 'cat %s | vcf-annotate --annotations %s --columns CHROM,POS,REF,ALT,-,-,-,FILTER > %s' % (vcf, tabixOut, annotateOut)
subprocess.call(annotateCmd, shell = True)

### Remove temporary folder
subprocess.call('rm -rf '+outPath, shell = True)
