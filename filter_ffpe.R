#!/usr/bin/env Rscript

##########################################################################################
# MSKCC CMO
# Identify samples with FFPE artifacts or filter the artifacts from a maf file.
##########################################################################################

strReverse <- function(x) sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")
revc <- function(x) strReverse(chartr('ACGT', 'TGCA', x))

# Annotate maf with Stratton Plot bin
add_mut_tri <- function(maf) {

  if (!"TriNuc" %in% names(maf)) {
	  if (!"flanking_bps" %in% names(maf)) {
		  stop("must have TriNuc/flanking_bps column")
	  }
  }

  ### check for t_var_freq
  if (!'t_var_freq' %in% names(maf))
    maf[!t_ref_count %in% c(NA,'.') & !t_alt_count %in% c(NA,'.'),
  t_var_freq := as.numeric(t_alt_count)/(as.numeric(t_alt_count)+as.numeric(t_ref_count))]

  ### reverse complement Ref_Tri if ref is either G or A
  if("TriNuc" %in% names(maf)){
  maf[Variant_Type == "SNP",
      Ref_Tri := ifelse(Reference_Allele %in% c('G', 'A'),
                        revc(TriNuc),
                        TriNuc)]
	}else{
		maf[Variant_Type == "SNP",
				Ref_Tri := ifelse(Reference_Allele %in% c('G', 'A'),
						revc(flanking_bps),
						flanking_bps)]
	}
  ### reverse complement Tumor_Seq_Allele2 if ref is either G or A
  Tumor_Seq_Allele2_CT <- maf[Variant_Type == "SNP",
                              ifelse(Reference_Allele %in% c('G', 'A'),
                                     revc(Tumor_Seq_Allele2),
                                     Tumor_Seq_Allele2)]

  ### combine Ref_Tri and Tumor_Seq_Allele2
  ### (with conditional reverse compliment)
  maf[Variant_Type == "SNP",
      Mut_Tri := paste0(substr(Ref_Tri, 1, 2),
                        Tumor_Seq_Allele2_CT,
                        substr(Ref_Tri, 3, 3))]

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

    pkgs = c('data.table', 'argparse')
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
