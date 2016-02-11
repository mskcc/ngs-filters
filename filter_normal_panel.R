#!/usr/bin/env Rscript

##########################################################################################
##########################################################################################
# MSKCC CMO
# Identify variant present in 3 or more alternate reads in at least 3 samples in
# normal panel and annotate somatic variant calls accordingly
##########################################################################################
##########################################################################################

annotate_maf <- function(maf, fillout, 
                         normal.count=3) {
  
  # select normal samples
  if (!is.na(normal.regex)) {
    fillout <- fillout[Tumor_Sample_Barcode %like% normal]
  } else {
    fillout <- fillout[Tumor_Sample_Barcode %like% "*N$"]
  }
  
  # identify loci with 3+ alternate reads in any normal sample
  fillout <- fillout[normal_count >= normal.count]
  
  # index
  fillout[, TAG := stringr::str_c('chr', Chromosome,
                                  ':', Start_Position,
                                  '-', End_Position,
                                  ':', Reference_Allele,
                                  ':', Tumor_Seq_Allele2)]
  
  normal_panel.blacklist <- unique(fillout$TAG)
  maf.annotated <- maf[, normal_panel := TAG %in% normal_panel.blacklist]
  
  return(maf.annotated)
  
}

parse_fillout <- function(fillout) {

  # Convert GetBaseCountsMultiSample output
  fillout = melt(fillout, id.vars = colnames(fillout)[1:34], variable.name = 'Tumor_Sample_Barcode') %>%
    separate(value, into = c('t_depth','t_ref_count','t_alt_count','t_var_freq'), sep = ';') %>%
    mutate(t_depth = str_extract(t_depth, regex('[0-9].*')) %>%
    mutate(t_ref_count = str_extract(t_ref_count, regex('[0-9].*')) %>%
    mutate(t_alt_count = str_extract(t_alt_count, regex('[0-9].*')) %>%
    mutate(t_var_freq = str_extract(t_var_freq, regex('[0-9].*')) %>%
    mutate(TAG = str_c('chr', Chromosome, ':', Start, '-', Start, ':', 'Ref', ':', 'Alt')) 

  # Note, variant might be present mutliple times if occuring in more than one sample // fix this at the fillout step
  # by de-duping the MAF?
  fillout = mutate(fillout, tmp_id = str_c(Tumor_Sample_Barcode, Chrom, Start, Ref, Alt, Gene))
  fillout = fillout[!duplicated(fillout$tmp_id),]

  # Calculate frequencies and return
  group_by(fillout, TAG) %>%
    summarize(normal_count = sum(t_alt_count>=3), avg_depth = mean(t_depth), avg_alt_count = mean(t_alt_count))

}

if( ! interactive() ) {
  
  pkgs = c('data.table', 'argparse', 'reshape2', 'dplyr', 'tidyr')
  junk <- lapply(pkgs, require, character.only = T)
  rm(junk)
  
  parser=ArgumentParser()
  parser$add_argument('-m', '--maf', type='character', help='SOMATIC_FACETS.vep.maf file')
  parser$add_argument('-f', '--fillout', type='character', help='GetBaseCountsMultiSample output')
  parser$add_argument('-n', '--normal_count', type='character', default=3, help='Normal panel count threshold')
  parser$add_argument('-o', '--outfile', type='character', help='Output file')
  args=parser$parse_args()
  
  maf <- fread(args$maf)
  fillout <- fread(args$fillout)
  alt.reads <- args$normals
  outfile <- args$outfile
  
  parsed_fillout = parse_fillout(fillout)

  maf.out <- annotate_maf(maf, parsed_fillout)
  write.table(maf.out, outfile, sep = "\t", 
              col.names = T, row.names = F,
              quote = F)
  
}