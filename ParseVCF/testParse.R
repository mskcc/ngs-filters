#!/usr/bin/env Rscript

args=commandArgs(trailing=F)
TAG="--file="
path_idx=grep(TAG,args)
SDIR=dirname(substr(args[path_idx],nchar(TAG)+1,nchar(args[path_idx])))
if(length(SDIR)==0) SDIR=getwd()

source(file.path(SDIR,"parseVCF.R"))

if( ! interactive() ) {

  pkgs = c('data.table', 'argparse', 'reshape2', 'dplyr', 'tidyr', 'stringr')
  junk <- lapply(pkgs, function(p){suppressPackageStartupMessages(require(p, character.only = T))})
  rm(junk)

  parser=ArgumentParser()
  parser$add_argument('-f', '--fillout', type='character', help='GetBaseCountsMultiSample output')
  parser$add_argument('-o', '--outfile', type='character', help='Output file', default = 'stdout')
  args=parser$parse_args()

  fillout <- suppressWarnings(fread(args$fillout, showProgress = F))
  outfile <- args$outfile

  fillMaf = parse_vcf(fillout)
  maf.out = fillMaf

  if (outfile == 'stdout') {
    write.table(maf.out, stdout(), sep = "\t", col.names = T, row.names = F, quote = F, na="")
  } else {
    write.table(maf.out, outfile, sep = "\t", col.names = T, row.names = F, quote = F, na="")
  }

}