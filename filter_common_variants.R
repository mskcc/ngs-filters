#!/opt/common/CentOS_6-dev/R/R-3.1.3/lib64/R/bin/Rscript

##########################################################################################
# MSKCC CMO
# Annotate common_variants identified by maf2maf
##########################################################################################

annotate_maf <- function(maf, flagged = NULL) {

    # check chrM
    # if (nrow(maf) == nrow(flagged) & length(setdiff(unique(maf$Chromosome), unique(flagged$Chromosome))) == 0) {

    #     maf.annotated <- maf[, common_variant := flagged$FILTER]

    # } else { # chrM dropped by maf2maf

    # index

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
    parser$add_argument('-f', '--flagged', type='character', help='vep.flagged.maf file', default = NULL)
    parser$add_argument('-o', '--outfile', type='character', help='Output file', default = 'stdout')
    args=parser$parse_args()

    if (args$maf == 'stdin') { maf = suppressWarnings(fread('cat /dev/stdin', showProgress = F))
    } else { maf <- suppressWarnings(fread(args$maf, showProgress = F)) }

    ### Check if user provided a tagged MAF, else if there's already an ExAC_AF column in MAF
    if (!is.null(args$flagged)) {
        flagged <- suppressWarnings(fread(args$flagged, showProgress = F))
        maf.out <- annotate_maf(maf, flagged)
    } else if ('ExAC_AF' %in% names(maf)) {
        maf.out = annotate_maf(maf)
    } else { stop('Needs ExAC_AF in either input MAF or supplementary flagged MAF.') }

    outfile <- args$outfile
    if (outfile == 'stdout') { write.table(maf.out, stdout(), sep = "\t", col.names = T, row.names = F, quote = F, na="")
    } else { write.table(maf.out, outfile, sep = "\t", col.names = T, row.names = F, quote = F, na="") }
}