#!/usr/bin/env Rscript

##########################################################################################
# MSKCC CMO
# Annotate low-confidence mutation calls
##########################################################################################

annotate_maf <- function(maf) {

    # Add TAG to MAF
    if (!('TAG' %in% names(maf))) {
        maf[, TAG := str_c('chr', Chromosome,
                           ':', Start_Position,
                           '-', End_Position,
                           ':', Reference_Allele,
                           ':', Tumor_Seq_Allele2)]
    }

    if (!('FILTER' %in% names(maf))) maf$FILTER = '.'
    maf.annotated <- maf[, low_confidence := (ifelse(n_alt_count =='.', FALSE, as.numeric(n_alt_count) > 1) | as.numeric(t_depth) < 20 | as.numeric(t_alt_count) <= 3)]
    maf.annotated <- maf[, FILTER := ifelse(FILTER == '.' & low_confidence == TRUE, 'low_confidence',
                                            ifelse(FILTER != '.' & low_confidence == TRUE,
                                                   paste0(FILTER, ',low_confidence'), FILTER))]
    return(maf.annotated)
}

if (!interactive()) {

    pkgs = c('data.table', 'argparse', 'stringr')
    junk <- lapply(pkgs, function(p){suppressPackageStartupMessages(require(p, character.only = T))})
    rm(junk)

    parser=ArgumentParser()
    parser$add_argument('-m', '--maf', type = 'character', help = 'SOMATIC_FACETS.vep.maf file', default = 'stdin')
    parser$add_argument('-o', '--outfile', type = 'character', help = 'Output file', default = 'stdout')
    args=parser$parse_args()

    if (args$maf == 'stdin') { maf = suppressWarnings(fread('cat /dev/stdin', colClasses=c(Chromosome="character"), showProgress = F))
    } else { maf <- suppressWarnings(fread(args$maf, colClasses=c(Chromosome="character"), showProgress = F)) }
    outfile <- args$outfile

    maf.out <- annotate_maf(maf)
    if (outfile == 'stdout') { write.table(maf.out, stdout(), na="", sep = "\t", col.names = T, row.names = F, quote = F)
    } else { write.table(maf.out, outfile, na="", sep = "\t", col.names = T, row.names = F, quote = F) }
}
