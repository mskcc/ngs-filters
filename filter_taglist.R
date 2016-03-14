#!/usr/bin/env Rscript

##########################################################################################
# MSKCC CMO
# Annotate MAF with arbitrary list of TAGs
##########################################################################################

annotate_maf <- function(maf, flagged = NULL) {

    if (!('TAG' %in% names(maf))) {
        maf[, TAG := stringr::str_c('chr', Chromosome,
                                    ':', Start_Position,
                                    '-', End_Position,
                                    ':', Reference_Allele,
                                    ':', Tumor_Seq_Allele2)]
    }

    if (!is.null(flagged)) {
        flagged[, TAG := stringr::str_c('chr', Chromosome,
                                        ':', Start_Position,
                                        '-', End_Position,
                                        ':', Reference_Allele,
                                        ':', Tumor_Seq_Allele2)]

        common_variants <- flagged[FILTER == "common_variant"]$TAG
        maf[, common_variant := "."]

    } else {
        common_variants <- maf[ExAC_AF <= .0004]$TAG
    }

    maf.annotated <- maf[, common_variant := TAG %in% common_variants]
    maf.annotated <- maf[, FILTER := ifelse(FILTER == '.' & common_variant == TRUE, 'common_variant',
                                            ifelse(FILTER != '.' & common_variant == TRUE,
                                                   paste0(FILTER, ',common_variant'), FILTER))]

    return(maf.annotated)

}

if( ! interactive() ) {

    pkgs = c('data.table', 'argparse')
    junk <- lapply(pkgs, function(p){suppressPackageStartupMessages(require(p, character.only = T))})
    rm(junk)

    parser=ArgumentParser()
    parser$add_argument('-m', '--maf', type='character', help='SOMATIC_FACETS.vep.maf file', default = 'stdin')
    parser$add_argument('-f', '--taglist', type='character', help='Headerless list of mutation TAGs', default = NULL)
    parser$add_argument('-o', '--outfile', type='character', help='Output file', default = 'stdout')
    args=parser$parse_args()

    if (args$maf == 'stdin') { maf = suppressWarnings(fread('cat /dev/stdin', showProgress = F))
    } else { maf <- suppressWarnings(fread(args$maf, showProgress = F)) }

    ### Check if user provided a tagged MAF, else if there's already an ExAC_AF column in MAF
    if (!is.null(args$taglist)) { stop('Needs list of mutation TAGs.\n') }
    taglist = fread(args$taglist, header = 'TAG')

    maf.out <- annotate_maf(maf, taglist)

    outfile <- args$outfile
    if (outfile == 'stdout') { write.table(maf.out, stdout(), sep = "\t", col.names = T, row.names = F, quote = F)
    } else { write.table(maf.out, outfile, sep = "\t", col.names = T, row.names = F, quote = F) }
}