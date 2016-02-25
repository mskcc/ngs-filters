#!/usr/bin/env Rscript

##########################################################################################
# MSKCC CMO
# Identify loci with 3 or more alternate reads in a FILLOUT .maf and
# annotate somatic variant calls accordingly
##########################################################################################

annotate_maf <- function(maf, fillout,
                         alt.reads = 3,
                         normal.regex = NA) {

  # select normal samples
  if (!is.na(normal.regex)) {
    fillout <- fillout[Tumor_Sample_Barcode %like% normal]
  } else {
    fillout <- fillout[Tumor_Sample_Barcode %like% "*N$"]
  }

  # identify loci with 3+ alternate reads in any normal sample
  fillout <- fillout[t_alt_count >= alt.reads]

  # Add TAG to MAF
  if (!('TAG' %in% names(maf))) {
    maf[, TAG := stringr::str_c('chr', Chromosome,
                                ':', Start_Position,
                                '-', End_Position,
                                ':', Reference_Allele,
                                ':', Tumor_Seq_Allele2)]
  }

  # index
  fillout[, TAG := stringr::str_c('chr', Chromosome,
                                  ':', Start_Position,
                                  '-', End_Position,
                                  ':', Reference_Allele,
                                  ':', Tumor_Seq_Allele2)]

  if (!('FILTER' %in% names(maf))) maf$FILTER = '.'
  fillout.blacklist <- unique(fillout$TAG)
  maf.annotated <- maf[, cohort_normal := TAG %in% fillout.blacklist]
  maf.annotated <- maf[, FILTER := ifelse(FILTER == '.' & cohort_normal == TRUE, 'cohort_normal',
                                          ifelse(FILTER != '.' & cohort_normal == TRUE,
                                                 paste0(FILTER, ',cohort_normal'), FILTER))]

  return(maf.annotated)

}

if( ! interactive() ) {

  pkgs = c('data.table', 'argparse')
  junk <- lapply(pkgs, function(p){suppressPackageStartupMessages(require(p, character.only = T))})
  rm(junk)

  parser=ArgumentParser()
  parser$add_argument('-m', '--maf', type='character', help='SOMATIC_FACETS.vep.maf file', default = 'stdin')
  parser$add_argument('-f', '--fillout', type='character', help='FILLOUT.vep.maf file')
  parser$add_argument('-n', '--reads', type='character', default=3, help='Alternate read threshold')
  parser$add_argument('-o', '--outfile', type='character', help='Output file', default = 'stdout')
  args=parser$parse_args()

  if (args$maf == 'stdin') { maf = suppressWarnings(fread('cat /dev/stdin', showProgress = F))
  } else { maf <- suppressWarnings(fread(args$maf, showProgress = F)) }
  fillout <- suppressWarnings(fread(args$fillout, showProgress = F))
  alt.reads <- args$reads
  outfile <- args$outfile

  maf.out <- annotate_maf(maf, fillout)
  if (outfile == 'stdout') { write.table(maf.out, stdout(), sep = "\t", col.names = T, row.names = F, quote = F)
  } else { write.table(maf.out, outfile, sep = "\t", col.names = T, row.names = F, quote = F) }
}