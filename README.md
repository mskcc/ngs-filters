# Modular filters for flagging of putatative false-positive mutation calls in whole-exome sequencing
![Alt Text](http://i.giphy.com/14bJDgZJb8SI4E.gif)

## Use
These scripts are intended to be used to add annotation to a MAF whether a given variant is a possible false positive. All take `stdin` and can write to `stdout` and are standalone with two exceptions, for which a fillout operation needs to be run. Filter flags are added to the `FILTER` column, in a comma-separated manner. This filters almost exclusively operate on SNVs. Additionally, this repo contains a wrapper for running a [VCF-based false-positive filter](https://github.com/ckandoth/variant-filter) which populates the FILTER field of a VCF file, which can be retained if conversion to MAF is carried out with [vcf2maf](https://github.com/mskcc/vcf2maf).

## Filters
* Common variants
A variant is considered common if its minor allele frequency in [ExAC](http://exac.broadinstitute.org/) exceeds 0.0004. This filter needs an `ExAC_AF` column which easiest is can be added to a MAF by running [maf2maf](https://github.com/mskcc/vcf2maf), which now also annotates the `FILTER` column. This hopefully will render this filter script obsolete. With the `-f` flag this filter will annotate a maf with information from another MAF.
```bash
./filter_common_variants.R -m input.maf -o output.maf
``` 

* Low-confidence calls
A variant is considered a low-confidence call if it fulfills `n_alt_count > 1 | t_depth < 20 | t_alt_count <= 3`. Interpretation and use of this filter depends on the nature of the sequencing experiment.
```bash
./filter_low_conf.R -m input.maf -o output.maf
```

* Presence in study normals
Flags a variant if it is supported by 3 reads or more in any of the normals sequenced in the same study. The cut-off for supporting reads can be set with the `-n` flag. See instructions below for how to generate a fillout file.
```bash
./filter_cohort_normals.R -m input.maf -o output.maf -f study.fillout
```

* Presence in pool of normals
Similarily to the previous filter, a variant is flagged by this filter if it is supported by 3 reads or more in at least 3 samples in a pool of normals. See instructions below for how to generate a fillout file.
```bash
./filter_norma_panel.R -m input.maf -o output.maf -f pon.fillout
```

* Presence in FFPE pool
Flags a variant if it is supported by 3 reads or more in a fillout against an FFPE pool. The cut-off for supporting reads can be set with the `-n` flag. See instructions below for how to generate a fillout file.
```bash
./filter_ffpe_pool.R -m input.maf -o output.maf -f ffpe.fillout
```

* FFPE artifact
Flags a variant if it looks like an FFPE artifact, i.e. occurs at low VAF and is a C>T substitution. This script also can help identifying samples suffering from FFPE artifacts by using the `-i` flag.
```bash
./filter_ffpe.R -m input.maf -o output.maf -i
 ### or
./filter_ffpe.R -m input.maf -o output.maf
```

* Low-mappability ("blacklisted") regions
Filter variants located in regions of to which sequencing reads are hard to map, as defined by [ENCODE](encodeproject.org/annotations/ENCSR636HFF/) and [RepeatMasker](http://www.repeatmasker.org/species/hg.html). See `data/source.txt` for details on the files used for this annotation. 
```bash
./filter_blacklist_regions.R -m input.maf -o output.maf
```

***
## Fillout wrapper
This script wraps `GetBaseCountsMultiSample` on luna and can be used to generate fillout files (i.e. allele counts for variants in input MAF) across a set of BAM files. The `-n` flag can be used to run multithreaded. The genome of the MAF and the BAMs needs to be consistent and specified with the `-g` flag, which knows where the assemblies for GRCh37, hg19, b37, and b37_dmp are located on luna. The script `convert-maf-to-hg19.sh` can be used to fake an hg19 MAF.
```bash
./maf_fillout.py -m input.maf -b file1.bam file2.bam [..] -g genome -n threads -o output.fillout
```

## fpfilter.pl wrapper
This script wraps `fpfilter.pl` from [variant-filter](https://github.com/ckandoth/variant-filter). Filter parameters in `fpfilter.pl` might be ajusted according to the nature of the sequencing experiment. Temporary files generated are removed upon completion. Like the fillout wrapper, this script knows where GRCh37, hg19, b37, and b37_dmp are located on luna.
```bash
./run-fpfilter.py -v input.vcf -b tumor.bam -g genome -f path/to/fpfilter.pl
```
