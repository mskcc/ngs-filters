#!/usr/bin/env Rscript

##########################################################################################
# MSKCC CMO
# Identify loci with 3 or more alternate reads in a FILLOUT .maf and
# annotate somatic variant calls accordingly
##########################################################################################

annotate_maf <- function(maf, fillout, normal.samples, alt.reads = 3) {

  # select normal samples
  fillout <- fillout[Tumor_Sample_Barcode %in% normal.samples]

  # identify loci with 3+ alternate reads in any normal sample
  fillout <- fillout[t_alt_count >= alt.reads]

  if(nrow(fillout)==0){
    cat("**WARNING** cohort normal fillout empty. Check normal samples or normal.regex=r/",normal.regex,"/\n")
  }

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
  parser$add_argument('-n', '--reads', type='double', default=3, help='Alternate read threshold')
  parser$add_argument('-o', '--outfile', type='character', help='Output file', default = 'stdout')
  parser$add_argument('-r', '--regex', type='character', help='RegEx for normals', default = '')
  parser$add_argument('-N', '--normalSamplesFile', type='character', help='File with list of normal samples', default='')
  args=parser$parse_args()

  if (args$maf == 'stdin') { maf = suppressWarnings(fread('cat /dev/stdin', showProgress = F))
  } else { maf <- suppressWarnings(fread(args$maf, colClasses=c(Chromosome="character"), showProgress = F)) }
  fillout <- suppressWarnings(fread(args$fillout, colClasses=c(Chromosome="character"),showProgress = F))
  alt.reads <- args$reads
  outfile <- args$outfile
  normal.regex <- args$regex
  normalSamplesFile <- args$normalSamplesFile

  # get normal samples

  if(normalSamplesFile!="") {

    normal.samples=scan(normalSamplesFile,"")

  } else {

    fillout.samples=unique(fillout$Tumor_Sample_Barcode)
    if (normal.regex=="") {
      normal.regex="*N$"
    }
    normal.samples=fillout.samples[grepl(normal.regex,fillout.samples)]
  }

  maf.out <- annotate_maf(maf, fillout, normal.samples)
  if (outfile == 'stdout') { write.table(maf.out, stdout(), na="", sep = "\t", col.names = T, row.names = F, quote = F, na="")
  } else { write.table(maf.out, outfile, na="", sep = "\t", col.names = T, row.names = F, quote = F, na="") }
}
