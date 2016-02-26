#!/usr/bin/env Rscript

suppressWarnings(library(data.table))
suppressWarnings(library(stringr))
args = commandArgs(TRUE)

maf = suppressWarnings(fread(args[1], showProgress = F))

maf[, TAG := str_c('chr', Chromosome,
					':', Start_Position,
                    '-', End_Position,
                    ':', Reference_Allele,
                    ':', Tumor_Seq_Allele2)]

maf = maf[!duplicated(TAG)]

write.table(maf, stdout(), sep = '\t', col.names = T, row.names = F, quote = F)