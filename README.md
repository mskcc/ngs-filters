# wes-filters

## Rules
* Never remove mutations from a maf
* Filters (in the GATK sense of not removing mutations from a maf) may add many columns but the last should be TRUE/FALSE (KEEP/REJECT?) best-guess keep decision.

## Filter modules for WES:
  1. Filter common variants according to ExAC/EVS/1000Genomes
  2. Filter against panel of WES normals
  3. Filter against cohort-specific panel of normals
  4. Filter FFPE artifacts (includes deep-sequenced FFPE pool)

![Alt Text](http://i.giphy.com/14bJDgZJb8SI4E.gif)
