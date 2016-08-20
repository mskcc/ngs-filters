#!/opt/common/CentOS_6-dev/R/R-3.1.3/lib64/R/bin/Rscript

suppressWarnings(library(data.table))
suppressWarnings(library(stringr))
args = commandArgs(TRUE)

if (length(args) < 1) stop('Usage: maf_uniq_tags.R input.maf')

maf = suppressWarnings(fread(args[1], showProgress = F))

maf[, TAG := str_c('chr', Chromosome,
					':', Start_Position,
                    '-', End_Position,
                    ':', Reference_Allele,
                    ':', Tumor_Seq_Allele2)]

maf = maf[!duplicated(TAG)]

write.table(maf, stdout(), sep = '\t', col.names = T, row.names = F, quote = F)
