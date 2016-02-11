#!/usr/bin/env Rscript

##########################################################################################
##########################################################################################
# MSKCC CMO
# Annotate common_variants identified by maf2maf
##########################################################################################
##########################################################################################

annotate_maf <- function(maf, flagged) {
  
  # check chrM
  if (nrow(maf) == nrow(flagged) & length(setdiff(unique(maf$Chromosome), unique(flagged$Chromosome))) == 0) {
  
      maf.annotated <- maf[, common_variant := flagged$FILTER]
  
  } else { # chrM dropped by maf2maf
    
    # index
    flagged[, TAG := stringr::str_c('chr', Chromosome,
                                    ':', Start_Position,
                                    '-', End_Position,
                                    ':', Reference_Allele,
                                    ':', Tumor_Seq_Allele2)]
    
    common_variants <- flagged[FILTER == "common_variant"]$TAG
    maf[, common_variant := "."]
    maf.annotated <- maf[, common_variant := TAG %in% common_variants]
  }
  
  return(maf.annotated)
  
}

if( ! interactive() ) {
  
  pkgs = c('data.table', 'argparse')
  junk <- lapply(pkgs, require, character.only = T)
  rm(junk)
  
  parser=ArgumentParser()
  parser$add_argument('-m', '--maf', type='character', help='SOMATIC_FACETS.vep.maf file')
  parser$add_argument('-f', '--flagged', type='character', help='vep.flagged.maf file')
  parser$add_argument('-o', '--outfile', type='character', help='Output file')
  args=parser$parse_args()
  
  maf <- fread(args$maf)
  flagged <- fread(args$flagged)
  outfile <- args$outfile
  
  maf.out <- annotate_maf(maf, flagged)
  write.table(maf.out, outfile, sep = "\t", 
              col.names = T, row.names = F,
              quote = F)
  
}