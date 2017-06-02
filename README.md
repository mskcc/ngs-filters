# Modular filters for flagging of putative false-positive mutation calls in whole-exome sequencing

## Use
These scripts are intended to be used to add annotation to a MAF whether a given variant is a possible false positive. All take `stdin` and can write to `stdout` and are standalone with two exceptions, for which a fillout operation needs to be run. Filter flags are added to the `FILTER` column, in a comma-separated manner. This filters almost exclusively operate on SNVs. Additionally, this repo contains a wrapper for running a [VCF-based false-positive filter](https://github.com/ckandoth/variant-filter) which populates the FILTER field of a VCF file, which can be retained if conversion to MAF is carried out with [vcf2maf](https://github.com/mskcc/vcf2maf).

## Generic script to apply filters: applyFilter.sh

This script is a wrapper which will run any of the R based filters in this repository. The output MAF is annotated with headers to indicate which filter was used and which version of the repository.

Usage:
```bash
	applyFilter.sh FILTER_NAME INPUT_MAF OUTPUT_MAF [Additional Parameters]
```

example:

```bash
	applyFilter.sh filter_blacklist_regions.R \
		Proj_1234_CMO_MAF.txt filteredMAF.txt
```

The first lines of the output MAF will look as follows:

```
#version 2.4
#wes-filters/applyFilter.sh VERSION=v1.0.1-2-g4d3694b FILTER=filter_blacklist_regions.R
```

## Script to run all wes-filter scripts 

This script currently runs the following scripts in given order using *applyFilter.sh*:
- tag_hotspots
- filter_blacklist_region
- filter_low_conf
- filter_ffpe
- filter_dmp
- filter_ffpe_pool (if fillout maf for ffpe sample is given)
- filter_normal_panel (if fillout maf for standard normal sample is given)
- filter_cohort_normal (if fillout maf for cohort normal sample is given)

Usage:
```
usage: run_wes-filters.py [options]

This tool helps to tag hotspot events

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         make lots of noise
  -m SomeID.maf, --input-maf SomeID.maf
                        Input maf file which needs to be tagged
  -o SomeID.maf, --output-maf SomeID.maf
                        Output maf file name
  -outdir /somepath/output, --outDir /somepath/output
                        Full Path to the output dir.
  -npmaf /somepath/to/normalpanel.maf, --normal-panel-maf /somepath/to/normalpanel.maf
                        Path to fillout maf file of panel of standard normals
  -fpmaf /somepath/to/ffpe_pool.maf, --ffpe_pool_maf /somepath/to/ffpe_pool.maf
                        Path to fillout maf file for FFPE artifacts
  -ncmaf /somepath/to/normalcohort.maf, --normal-cohort-maf /somepath/to/normalcohort.maf
                        Path to fillout maf file of cohort normals
  -nsf /somepath/to/normalcohort.list, --normalSamplesFile /somepath/to/normalcohort.list
                        File with list of normal samples
  -hsp SomeID.txt, --input-hotspot SomeID.txt
                        Input txt file which has hotspots
                        
```

example:
```
python /home/shahr2/git/wes-filters/run_wes-filters.py -m output.maf -o output_wes.maf -npmaf /ifs/work/prism/shahr2/cmo_fill/output_fill.maf -f/ifs/work/prism/shahr2/cmo_fill/output_FFPE.maf -hsp /home/shahr2/git/hotspot-whitelist/v2/hotspot-list-union-v1-v2.txt -v
```

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
./filter_normal_panel.R -m input.maf -o output.maf -f pon.fillout
```

* Presence in FFPE pool
Flags a variant if it is supported by 3 reads or more in a fillout against an FFPE pool. The cut-off for supporting reads can be set with the `-n` flag. See instructions below for how to generate a fillout file.
```bash
./filter_ffpe_pool.R -m input.maf -o output.maf -f ffpe.fillout
```

* FFPE artifact
Flags a variant if it looks like an FFPE artifact, i.e. occurs at low VAF and is a C>T substitution. This script also can help identifying samples suffering from FFPE artifacts by using the `-i` flag.
```bash
./filter_ffpe.R -m input.maf -o output.maf
 ### or
./filter_ffpe.R -m input.maf -i
```

* Low-mappability ("blacklisted") regions
Filter variants located in regions of to which sequencing reads are hard to map, as defined by [ENCODE](encodeproject.org/annotations/ENCSR636HFF/) and [RepeatMasker](http://www.repeatmasker.org/species/hg.html). See `data/source.txt` for details on the files used for this annotation.
```bash
./filter_blacklist_regions.R -m input.maf -o output.maf
```

* Hotspots ("whitelisted") sites
Tag variants located in sites define as hotspots by [hotspot-whitelist](https://github.com/mskcc/hotspot-whitelist)
```python
./tag_hotspots.py -m input.maf -itxt hotspot-list-union-v1-v2.txt -o output.maf
```

* Add DMP Filter Tag
Flags a variant if it dose not pass allele count threshold set at DMP i.e for hotspot:AD>=8;DP>=20;VF=>0.02 & for non-hotspots D>=10;DP>=20;VF=>0.05 occurs at low VAF and is a C>T substitution. This script also can help identifying samples suffering from FFPE artifacts by using the `-i` flag.
```bash
./filter_dmp.R -m input.maf -o output.maf
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
