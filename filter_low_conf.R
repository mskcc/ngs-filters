#!/usr/bin/env Rscript

##########################################################################################
##########################################################################################
# MSKCC CMO
# Annotate common_variants identified by maf2maf
##########################################################################################
##########################################################################################

annotate_maf <- function(maf, flagged) {
    # Add TAG to MAF
    if (!('TAG' %in% names(maf))) {
        maf[, TAG := str_c('chr', Chromosome,
                           ':', Start_Position,
                           '-', End_Position,
                           ':', Reference_Allele,
                           ':', Tumor_Seq_Allele2)]
    }

    maf.annotated <- maf[, low_confidence := (n_alt_count > 1 | t_depth < 20 | t_alt_count <= 3)]
    maf.annotated <- maf[, FILTER := ifelse(FILTER == '.' & low_confidence == TRUE, 'low_confidence',
                                            ifelse(FILTER != '.' & low_confidence == TRUE,
                                                   paste0(FILTER, ',low_confidence'), FILTER))]

    return(maf.annotated)
}

if( ! interactive() ) {

    pkgs = c('data.table', 'argparse')
    junk <- lapply(pkgs, function(p){suppressPackageStartupMessages(require(p, character.only = T))})
    rm(junk)

    parser=ArgumentParser()
    parser$add_argument('-m', '--maf', type='character', help='SOMATIC_FACETS.vep.maf file')
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