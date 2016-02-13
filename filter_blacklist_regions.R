#!/usr/bin/env Rscript

##########################################################################################
##########################################################################################
# MSKCC CMO
# Annotate DAC blacklisted regions obtained from:
# https://www.encodeproject.org/annotations/ENCSR636HFF/
##########################################################################################
##########################################################################################

annotate_maf <- function(maf, blacklist) {
  
  setkey(blacklist, Chromosome, Start_Position, End_Position)
  
  fo <- foverlaps(maf[, .(Chromosome, Start_Position, End_Position)], 
                  blacklist, 
                  type="any")
  
  maf.annotated <- maf[, blacklist_region := fo$Info]
  
  return(maf.annotated)
  
}

if( ! interactive() ) {
  
  pkgs = c('data.table', 'argparse')
  junk <- lapply(pkgs, require, character.only = T)
  rm(junk)
  
  parser=ArgumentParser()
  parser$add_argument('-m', '--maf', type='character', help='SOMATIC_FACETS.vep.maf file')
  parser$add_argument('-b', '--blacklist', type='character', help='DAC Blacklisted Regions')
  parser$add_argument('-o', '--outfile', type='character', help='Output file')
  args=parser$parse_args()
  
  maf <- fread(args$maf)
  blacklist <- fread(args$blacklist)
  outfile <- args$outfile
  
  blacklist[, c('V5', 'V6') := NULL]
  setnames(blacklist, c("Chromosome", "Start_Position", "End_Position", "Info"))
  blacklist[, Chromosome := gsub("chr", "", Chromosome)]
  
  maf.out <- annotate_maf(maf, blacklist)
  write.table(maf.out, outfile, sep = "\t", 
              col.names = T, row.names = F,
              quote = F)
  
}
