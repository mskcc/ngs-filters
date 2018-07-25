#!/usr/bin/env Rscript

##########################################################################################
# MSKCC CMO
# Identify variants with at least 3 reads of support in at least 5 (default) samples in a
# panel of curated normals, and annotate these calls with FILTER tag "normal_panel"
##########################################################################################

annotate_maf <- function(maf, fillout, normal.count=8) {

    # Subset to loci with read support in at least 5 samples in the normal panel
    fillout <- fillout[fillout$threshold > normal.count,]

    # Add TAG to MAF
    # if (!('TAG' %in% names(maf))) {
    #     maf[, TAG := stringr::str_c('chr', Chromosome,
    #                     ':', Start_Position,
    #                     '-', End_Position,
    #                     ':', Reference_Allele,
    #                     ':', Tumor_Seq_Allele2)]
    # }
    # maf[, tmp_id := stringr::str_c('chr', Chromosome,
    #                 ':', Start_Position,
    #                 '-', End_Position,
    #                 ':', Reference_Allele,
    #                 ':', Tumor_Seq_Allele1,
    #                 ':', Tumor_Sample_Barcode)]


    if (!('FILTER' %in% names(maf))) maf$FILTER = '.'
    # maf[is.na(maf)] <- '.'

    normal_panel.blacklist <- unique(fillout$tmp_id)
    maf.annotated <- maf[, normal_panel := tmp_id %in% normal_panel.blacklist]
    maf.annotated <- maf[, FILTER := ifelse(normal_panel == TRUE, ifelse((FILTER == '' | FILTER == '.' | FILTER == 'PASS' | is.na(FILTER) ), 'normal_panel', paste0(FILTER, ';normal_panel')), FILTER)]
    return(maf.annotated)
}


# parse_fillout <- function(fillout) {

#     # Convert GetBaseCountsMultiSample output
#     fillout = melt(fillout, id.vars = colnames(fillout)[1:34], variable.name = 'Tumor_Sample_Barcode') %>%
#             separate(value, into = c('n_depth','n_ref_count','n_alt_count','n_var_freq'), sep = ';') %>%
#             mutate(n_depth = str_extract(n_depth, regex('[0-9].*'))) %>%
#             mutate(n_ref_count = str_extract(n_ref_count, regex('[0-9].*'))) %>%
#             mutate(n_alt_count = str_extract(n_alt_count, regex('[0-9].*'))) %>%
#             mutate(n_var_freq = str_extract(n_var_freq, regex('[0-9].*'))) %>%
#             mutate(TAG = stringr::str_c('chr', Chrom, ':', Start, '-', Start, ':', Ref, ':', Alt))

#     # Note, variant might be present multiple times if occuring in more than one sample, fix this
#     # at the fillout step by de-duping the MAF
#     fillout = mutate(fillout, tmp_id = stringr::str_c(Tumor_Sample_Barcode, Chrom, Start, Ref, Alt, Gene))
#     fillout = fillout[!duplicated(fillout$tmp_id),]

#     # Calculate frequencies and return
#     return(group_by(fillout, TAG) %>% summarize(normal_count = sum(n_alt_count>=3)))
# }

parse_fillout_maf <- function(maf, fillout) {
    # index
    fillout[, TAG := stringr::str_c('chr', Chromosome,
                    ':', Start_Position,
                    '-', End_Position,
                    ':', Reference_Allele,
                    ':', Tumor_Seq_Allele1)]
    fillout[, tmp_id := stringr::str_c('chr', Chromosome,
                    ':', Start_Position,
                    '-', End_Position,
                    ':', Reference_Allele,
                    ':', Tumor_Seq_Allele1,
                    ':', Tumor_Sample_Barcode)]
    fillout = fillout[!duplicated(fillout$tmp_id),]

    if (!('TAG' %in% names(maf))) {
        maf[, TAG := stringr::str_c('chr', Chromosome,
                        ':', Start_Position,
                        '-', End_Position,
                        ':', Reference_Allele,
                        ':', Tumor_Seq_Allele2)]
    }
    maf[, tmp_id := stringr::str_c('chr', Chromosome,
                    ':', Start_Position,
                    '-', End_Position,
                    ':', Reference_Allele,
                    ':', Tumor_Seq_Allele1,
                    ':', Tumor_Sample_Barcode)]
    # Calculate frequencies and return
    maf$vaf <- maf$t_alt_count / maf$t_depth
    maf$ptvf <- maf$vaf / 10
    maf.shortlist<-select(maf,TAG,tmp_id,ptvf)

    normpanel<-select(fillout,TAG,t_variant_frequency,t_alt_count)
    normpanel <- normpanel[normpanel$t_alt_count >= 3,]
    fulljoin.maf<-full_join(maf.shortlist,normpanel,by='TAG')
    fulljoin.maf$threshold <- fulljoin.maf$t_variant_frequency > fulljoin.maf$ptvf
    
    return(group_by(fulljoin.maf,tmp_id) %>% summarize(threshold=sum(threshold)))
}

if( ! interactive() ) {

    pkgs = c('data.table', 'argparse', 'reshape2', 'dplyr', 'tidyr', 'stringr')
    junk <- lapply(pkgs, function(p){suppressPackageStartupMessages(require(p, character.only = T))})
    rm(junk)

    parser=ArgumentParser()
    parser$add_argument('-m', '--maf', type='character', help='SOMATIC_FACETS.vep.maf file', default = 'stdin')
    parser$add_argument('-f', '--fillout', type='character', help='GetBaseCountsMultiSample output')
    parser$add_argument('-fo', '--fillout_format', type='double', help='GetBaseCountsMultiSample output format. MAF(1), Tab-delimited with VCF coordinates (2:default)', default=2)
    parser$add_argument('-n', '--normal_count', type='double', default=8, help='Minimum number of normal panel samples with 3+ supporting reads')
    parser$add_argument('-o', '--outfile', type='character', help='Output file', default = 'stdout')
    args=parser$parse_args()

    maf <- suppressWarnings(fread(args$maf, colClasses=c(Chromosome="character"), showProgress = F))
    fillout <- suppressWarnings(fread(args$fillout, colClasses=c(Chromosome="character"), showProgress = F))
    fillout.format<-args$fillout_format
    normal.count <- args$normal_count
    outfile <- args$outfile

    # if(fillout.format == 2) {
    #     parsed_fillout = parse_fillout(fillout)
    # }
    # else {
        parsed_fillout = parse_fillout_maf(maf,fillout)
    # }
    maf.out <- annotate_maf(maf, parsed_fillout, normal.count)
    if (outfile == 'stdout') {
        write.table(maf.out, stdout(), na="", sep = "\t", col.names = T, row.names = F, quote = F)
    }
    else {
        write.table(maf.out, outfile, na="", sep = "\t", col.names = T, row.names = F, quote = F)
    }
}
