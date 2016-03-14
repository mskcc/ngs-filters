#!/usr/bin/env Rscript

##########################################################################################
# MSKCC CMO
# Identify samples with FFPE artifacts or filter the artifacts from a maf file.
##########################################################################################

# Annotate maf with Stratton Plot bin
add_mut_tri <- function(maf) {

    if (!"TriNuc" %in% names(maf)) {
        if ("Ref_Tri" %in% names(maf)) {
            maf[, TriNuc := Ref_Tri]
        } else {
        stop("must have either Ref_Tri or TriNuc column")
        }
    }

  if (!'t_var_freq' %in% names(maf)) maf[, t_var_freq := t_alt_count/(t_alt_count+t_ref_count)]

  maf[, c('TriNuc_CT', 'Tumor_Seq_Allele2_CT', 'Mut_Tri') := 'X']
  maf[Variant_Type == "SNP" & !is.na(Reference_Allele) & !is.na(t_var_freq) & !is.na(TriNuc),
      TriNuc_CT := ifelse(Reference_Allele %in% c("G", "A"),
                          as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(TriNuc))),
                          TriNuc)]

  maf[Variant_Type == "SNP" & !is.na(Reference_Allele) & !is.na(t_var_freq) & !is.na(Tumor_Seq_Allele2),
      Tumor_Seq_Allele2_CT := ifelse(Reference_Allele %in% c("G", "A"),
                                     as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(Tumor_Seq_Allele2))),
                                     Tumor_Seq_Allele2)]

  maf[Variant_Type == "SNP" & !is.na(TriNuc_CT) & !is.na(Tumor_Seq_Allele2_CT),
      Mut_Tri := paste0(substr(TriNuc_CT, 1, 2),
                        Tumor_Seq_Allele2_CT,
                        substr(TriNuc_CT, 3, 3))]

  maf
}


# For each sample, this function returns
# the fraction of low allele fraction mutations that have NCG > NTG context
# Useful to identify samples with FFPE artifacts.
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

# Remove low allele fraction mutations with the C>T context
filter_artifacts <- function(maf, threshold = 0.1) {

    if (!('FILTER' %in% names(maf))) maf$FILTER = '.'
    maf.annotated = maf[, ffpe_artifact := t_var_freq <= threshold & Mut_Tri %like% ".CT."]
    maf.annotated <- maf[, FILTER := ifelse(FILTER == '.' & ffpe_artifact == TRUE, 'ffpe_artifact',
                                            ifelse(FILTER != '.' & ffpe_artifact == TRUE,
                                                   paste0(FILTER, ',ffpe_artifact'), FILTER))]
    maf.annotated
}


if (!interactive()) {

    pkgs = c('data.table', 'argparse', 'Biostrings')
    junk <- lapply(pkgs, function(p){suppressPackageStartupMessages(require(p, character.only = T))})
    rm(junk)

    parser=ArgumentParser()
    parser$add_argument('-m', '--maf', type = 'character', help = 'SOMATIC_FACETS.vep.maf file', default = 'stdin')
    parser$add_argument('-t', '--threshold', type = 'double', default = 10, help = 'allele fraction threshold (pct.)')
    parser$add_argument('-i', '--identify', action = "store_true", default = FALSE, help = 'Identify samples with FFPE artifacts')
    parser$add_argument('-f', '--filter', action = "store_true", default = TRUE, help = 'Filter FFPE artifacts')
    parser$add_argument('-o', '--outfile', type = 'character', help = 'Output file', default = 'stdout')

    args=parser$parse_args()
    outfile = args$outfile

    threshold <- as.numeric(args$threshold) / 100

    if(args$identify & args$filter) stop("Cannot identify and filter at the same time")

    if (args$maf == 'stdin') { maf = suppressWarnings(fread('cat /dev/stdin', showProgress = F))
    } else { maf <- suppressWarnings(fread(args$maf, showProgress = F)) }

    if(args$identify) {
        maf <- add_mut_tri(maf)
        m <- identify_artifacts(maf, threshold = threshold)
        write.table(format(m, digits=2), stdout(), col.names = T, row.names = F, sep = '\t', quote = F)
    } else if(args$filter) {
        maf <- add_mut_tri(maf)
        maf.out <- filter_artifacts(maf)
        if (outfile == 'stdout') { write.table(maf.out, stdout(), sep = "\t", col.names = T, row.names = F, quote = F)
        } else { write.table(maf.out, outfile, sep = "\t", col.names = T, row.names = F, quote = F) }
    }
}
