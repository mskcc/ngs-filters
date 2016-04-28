parse_vcf <- function(vcf){

    format = strsplit(vcf$FORMAT[1],":")[[1]]


    # Unroll VCF samples columns into MAF rows and fix #CHROM column name

    maf = melt(vcf, id.vars = colnames(fillout)[1:9], variable.name = 'Tumor_Sample_Barcode') %>%
        mutate(CHROM=`#CHROM`)

    # convert events to MAF format convention
    #   INS X XN ==> - N
    #   DEL XN X ==> N - (also need to fix position)

    # TODO
    # N.B. *BUG* This code does not work for more complex events like
    #   XNM XN
    #   XN  XNM

    # First compute type of event by sign of change in length from REF to ALT
    # +1 == INS
    # -1 == DEL
    #  0 == SNP

    maf=mutate(maf,mSignDelta=sign(str_length(ALT)-str_length(REF))) %>%

        # Remove the leading base (convert to '-' if len==1)

        mutate(REF=ifelse(mSignDelta==0,REF,ifelse(mSignDelta>0,"-",str_sub(REF,2,-1)))) %>%
        mutate(ALT=ifelse(mSignDelta==0,ALT,ifelse(mSignDelta<0,"-",str_sub(ALT,2,-1)))) %>%

        # For deletions need to fix POS

        mutate(POS=ifelse(mSignDelta < 0, POS + 1, POS)) %>%

        # Parse format column of VCF

        separate(value, into = format, sep = ':') %>%

        # New format for index TAG. Although we are in MAF coordinates
        # only use start point (POS). Too lazy to compute end point
        # and it is not needed

        mutate(TAG = str_c(CHROM, ':', POS, ':', REF, ':', ALT))

    # TODO
    # Bug in current fillout code when given a MAF for input
    # which causes lines to be duplicated if the same event
    # appears multiple times. This should be fixed in fillout
    # However STAG might be usefull

    maf = mutate(maf, STAG = str_c(TAG,":",Tumor_Sample_Barcode))
    maf = maf[!duplicated(maf$STAG),]

}

