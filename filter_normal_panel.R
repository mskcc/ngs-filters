#!/usr/bin/env Rscript

##########################################################################################
# MSKCC CMO
# Identify variant present in 3 or more alternate reads in at least 3 samples in
# normal panel and annotate somatic variant calls accordingly
##########################################################################################

annotate_maf <- function(maf, fillout,
                         normal.count=3) {

  # identify loci with 3+ alternate reads in any normal sample
  fillout <- fillout[normal_count >= normal.count]

  # Add TAG to MAF
  if (!('TAG' %in% names(maf))) {
    maf[, TAG := str_c(Chromosome,
                       ':', Start_Position,
                       '-', End_Position,
                       ':', Reference_Allele,
                       ':', Tumor_Seq_Allele2)]
  }

  if (!('FILTER' %in% names(maf))) maf$FILTER = '.'
  normal_panel.blacklist <- unique(fillout$TAG)
  maf.annotated <- maf[, normal_panel := TAG %in% normal_panel.blacklist]
  maf.annotated <- maf[, FILTER := ifelse(FILTER == '.' & normal_panel == TRUE, 'normal_panel',
                                          ifelse(FILTER != '.' & normal_panel == TRUE,
                                                 paste0(FILTER, ',normal_panel'), FILTER))]

  return(maf.annotated)

}

parse_fillout <- function(fillout) {

  # Convert GetBaseCountsMultiSample output
  fillout = melt(fillout, id.vars = colnames(fillout)[1:34], variable.name = 'Tumor_Sample_Barcode') %>%
    separate(value, into = c('n_depth','n_ref_count','n_alt_count','n_var_freq'), sep = ';') %>%
    mutate(n_depth = str_extract(n_depth, regex('[0-9].*'))) %>%
    mutate(n_ref_count = str_extract(n_ref_count, regex('[0-9].*'))) %>%
    mutate(n_alt_count = str_extract(n_alt_count, regex('[0-9].*'))) %>%
    mutate(n_var_freq = str_extract(n_var_freq, regex('[0-9].*'))) %>%
    mutate(TAG = str_c(Chrom, ':', Start, '-', Start, ':', Ref, ':', Alt))

  # Note, variant might be present mutliple times if occuring in more than one sample, fix this at the fillout step
  # by de-duping the MAF?
  fillout = mutate(fillout, tmp_id = str_c(Tumor_Sample_Barcode, Chrom, Start, Ref, Alt, Gene))
  fillout = fillout[!duplicated(fillout$tmp_id),]

  # Calculate frequencies and return
  group_by(fillout, TAG) %>%
    summarize(normal_count = sum(n_alt_count>=3))
}

if( ! interactive() ) {

  pkgs = c('data.table', 'argparse', 'reshape2', 'dplyr', 'tidyr', 'stringr')
  junk <- lapply(pkgs, function(p){suppressPackageStartupMessages(require(p, character.only = T))})
  rm(junk)

  parser=ArgumentParser()
  parser$add_argument('-m', '--maf', type='character', help='SOMATIC_FACETS.vep.maf file', default = 'stdin')
  parser$add_argument('-f', '--fillout', type='character', help='GetBaseCountsMultiSample output')
  parser$add_argument('-n', '--normal_count', type='character', default=3, help='Normal panel count threshold')
  parser$add_argument('-o', '--outfile', type='character', help='Output file', default = 'stdout')
  args=parser$parse_args()

  maf <- suppressWarnings(fread(args$maf, showProgress = F))
  fillout <- suppressWarnings(fread(args$fillout, showProgress = F))
  alt.reads <- args$normals
  outfile <- args$outfile

  parsed_fillout = parse_fillout(fillout)

  maf.out <- annotate_maf(maf, parsed_fillout)
  if (outfile == 'stdout') { write.table(maf.out, stdout(), sep = "\t", col.names = T, row.names = F, quote = F, na="")
  } else { write.table(maf.out, outfile, sep = "\t", col.names = T, row.names = F, quote = F, na="") }
}
