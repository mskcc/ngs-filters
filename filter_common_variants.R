#!/usr/bin/env Rscript

##########################################################################################
# MSKCC CMO
# Annotate common_variants identified by maf2maf
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
  junk <- lapply(pkgs, function(p){suppressPackageStartupMessages(require(p, character.only = T))})
  rm(junk)

  parser=ArgumentParser()
  parser$add_argument('-m', '--maf', type='character', help='SOMATIC_FACETS.vep.maf file', default = 'stdin')
  parser$add_argument('-f', '--flagged', type='character', help='vep.flagged.maf file')
  parser$add_argument('-o', '--outfile', type='character', help='Output file', default = 'stdout')
  args=parser$parse_args()

  if (args$maf == 'stdin') { maf = suppressWarnings(fread('cat /dev/stdin', showProgress = F))
  } else { maf <- suppressWarnings(fread(args$maf, showProgress = F)) }
  flagged <- suppressWarnings(fread(args$flagged, showProgress = F))
  outfile <- args$outfile

  maf.out <- annotate_maf(maf, flagged)
  if (outfile == 'stdout') { write.table(maf.out, stdout(), sep = "\t", col.names = T, row.names = F, quote = F)
  } else { write.table(maf.out, outfile, sep = "\t", col.names = T, row.names = F, quote = F) }
}