#!/opt/common/CentOS_6-dev/bin/current/python

#######################################
### Input:
### MAF file
### list of normal BAMs
### list of tumor BAMs

#######################################
### 1. Check if run through maf2maf:
### if not, do so
### once done, continue to next filter
#######################################

#######################################
### 2. Filter against panel of WES normals
### Check if MAF fillout exists, if not, do so
### once done, parse fillout and compile stats called variants in normals
### then proceed
#######################################

#######################################
### 3. Filter against cohort-specific panel of normals
### Check if MAF fillout exists, if not, do so
### once done, parse fillout and compile stats called variants in normals
### then proceed
#######################################

#######################################
### 4. FFPE filter
#######################################
