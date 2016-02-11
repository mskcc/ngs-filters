#!/usr/bin/env Rscript

##########################################################################################
##########################################################################################
# Alex Penson
# MSKCC CMO
# Identify samples with FFPE artifacts
# or filter the artifacts from a maf file.
##########################################################################################
##########################################################################################

write.maf <- function (...) {
  write.table(..., quote = F, col.names = T, row.names = F,
              sep = "\t")
}

#' Annotate maf with Stratton Plot bin
add_Mut_Tri <- function(maf){
  if(! "TriNuc" %in% names(maf)){
    if("Ref_Tri" %in% names(maf)){
      maf[, TriNuc := Ref_Tri]
    } else {
      stop("must have either Ref_Tri or TriNuc column")
    }
  }
  maf <- maf[Variant_Type == "SNP"]
  maf[, TriNuc_CT := ifelse(Reference_Allele %in% c("G", "A"),
                            as.character(
                              Biostrings::reverseComplement(
                                Biostrings::DNAStringSet(TriNuc)
                              )
                            ),
                            TriNuc)
      ]
  maf[, Tumor_Seq_Allele2_CT := ifelse(Reference_Allele %in% c("G", "A"),
                                       as.character(
                                         Biostrings::reverseComplement(
                                           Biostrings::DNAStringSet(Tumor_Seq_Allele2)
                                         )
                                       ),
                                       Tumor_Seq_Allele2)
      ]

  maf[, Mut_Tri := paste0(substr(TriNuc_CT, 1, 2),
                          Tumor_Seq_Allele2_CT,
                          substr(TriNuc_CT, 3, 3))]
  maf
}


#' For each sample, this function returns
#' the fraction of low allele fraction mutations that have NCG > NTG context
#' Useful to identify samples with FFPE artifacts.
identify_artifacts <- function(maf, threshold = 0.1){
  m <- merge(
    maf[t_var_freq < threshold, list(frac_lo = mean(Mut_Tri %like% "CTG$")),
        keyby = Tumor_Sample_Barcode],
    maf[t_var_freq > threshold, list(frac_hi = mean(Mut_Tri %like% "CTG$")),
        keyby = Tumor_Sample_Barcode])
  m[, frac_diff := frac_lo - frac_hi]
  m <- m[order(frac_diff)]
  m
}

#' Remove low allele fraction mutations with the FFPE context
filter_artifacts <- function(maf, threshold = 0.1){
    maf <- maf[t_var_freq > threshold | ! Mut_Tri %like% "CTG$"]
    maf
}


if( ! interactive() ) {

  pkgs = c('data.table', 'argparse')
  junk <- lapply(pkgs, function(p){suppressPackageStartupMessages(require(p, character.only = T))})
  rm(junk)

  parser=ArgumentParser()
  parser$add_argument('-m', '--maf', type='character',
                      help='SOMATIC_FACETS.vep.maf file')
  parser$add_argument('-t', '--threshold', type='double', default = 10,
                      help='allele fraction threshold (%)')
  parser$add_argument('--identify', action="store_true", default = FALSE,
                      help='Identify samples with FFPE artifacts')
  parser$add_argument('--filter', action="store_true", default = FALSE,
                      help='Filter FFPE artifacts')
  args=parser$parse_args()

  threshold <- as.numeric(args$threshold) / 100

  if(args$identify & args$filter){
    stop("Cannot identify and filter at the same time")
  }

  maf <- suppressWarnings(fread(args$maf, showProgress = F))
  if(args$identify) {
    maf <- add_Mut_Tri(maf)
    m <- identify_artifacts(maf, threshold = threshold)
    write.maf(format(m, digits=2), stdout())
  } else if(args$filter) {
    maf <- filter_artifacts(maf)
    write.maf(maf.out, stdout())
  }

}
