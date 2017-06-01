#!/usr/bin/python
'''
@Description : This tool helps to run all scripts in wes-filter
@Created :  05/30/2017
@Updated : 05/30/2017
@author : Ronak H Shah

NORMALCOHORTBAMS="/ifs/res/share/pwg/NormalCohort/SetA/CuratedBAMsSetA"
FFPEPOOLDIR="/ifs/res/share/soccin/Case_201601/Proj_06049_Pool/r_001"
FFPEPOOLBAM=os.path.join(FFPEPOOLDIR,"/alignments/Proj_06049_Pool_indelRealigned_recal_s_UD_ffpepool1_N.bam")

'''
from __future__ import division
import argparse
import sys
import os
import time
import logging
import tempfile
import subprocess

logging.basicConfig(
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=logging.DEBUG)
logger = logging.getLogger('run_wes-filters')
try:
    import coloredlogs
    coloredlogs.install(level='DEBUG')
except ImportError:
    logger.warning("run_wes-filters: coloredlogs is not installed, please install it if you wish to see color in logs on standard out.")
    pass
try:
    import pandas as pd
except ImportError:
    logger.fatal("run_wes-filters: pandas is not installed, please install pandas as it is required to run the process.")
    sys.exit(1)

def main():
   parser = argparse.ArgumentParser(prog='run_wes-filters.py', description=' This tool helps to tag hotspot events', usage='%(prog)s [options]')
   parser.add_argument("-v", "--verbose", action="store_true", dest="verbose", help="make lots of noise")
   parser.add_argument("-m", "--input-maf", action="store", dest="inputMaf", required=True, type=str, metavar='SomeID.maf', help="Input maf file which needs to be tagged")
   parser.add_argument("-o","--output-maf", action="store", dest="outputMaf", required=True, type=str, metavar='SomeID.maf', help="Output maf file name")
   parser.add_argument("-outdir", "--outDir", action="store", dest="outdir", required=False, type=str, metavar='/somepath/output', help="Full Path to the output dir.")
   parser.add_argument("-npmaf","--normal-panel-maf",action="store", dest="NormalPanelMaf", required=False, type=str, metavar='/somepath/to/normalpanel.maf', help="Path to fillout maf file of panel of standard normals")
   parser.add_argument("-fpmaf", "--ffpe_pool_maf", action="store", dest="FFPEPoolMaf", required=False, type=str, metavar='/somepath/to/ffpe_pool.maf', help="Path to fillout maf file for FFPE artifacts")
   parser.add_argument("-ncmaf","--normal-cohort-maf",action="store", dest="NormalCohortMaf", required=False, type=str, metavar='/somepath/to/normalcohort.maf', help="Path to fillout maf file of cohort normals")
   parser.add_argument("-hsp", "--input-hotspot", action="store", dest="inputHSP", required=True, type=str, metavar='SomeID.txt', help="Input txt file which has hotspots")
   
   args = parser.parse_args()
   if(args.verbose):
       logger.info("run_wes-filters: Started the run for wes-filters")
   (finalmaf) = run_wes_filters(args)
   if(args.verbose):
       logger.info("run_wes-filters: Output is written in %s", finalmaf)
       logger.info("run_wes-filters: Finished the run for wes-filters.")


def run_wes_filters(args):
    tmpdir = tempfile.mkdtemp()
    wes_filter_bin = os.path.dirname(os.path.realpath(__file__))
    tag_hotspot = os.path.join(wes_filter_bin,"tag_hotspots.py")
    apply_filter = os.path.join(wes_filter_bin,"applyFilter.sh")
    
    #tag_hotspots
    if(args.verbose):
        logger.info("run_wes-filters: Tagging Hotspots")
    tempMaf0 = os.path.join(tmpdir,"maf0.maf")
    cmd =  tag_hotspot + " -v -m " + args.inputMaf  + " -itxt " + args.inputHSP + " -o " + tempMaf0
    if(args.verbose):
        logger.info("run_wes-filters: Running, %s",cmd)
    subprocess.call(cmd, shell = True)
    
    #black_list_region
    if(args.verbose):
        logger.info("run_wes-filters: Applying filter_blacklist_regions")
    tempMaf1 = os.path.join(tmpdir,"maf1.maf")
    cmd =  apply_filter + " filter_blacklist_regions.R " + tempMaf0  + " " + tempMaf1
    if(args.verbose):
        logger.info("run_wes-filters: Running, %s",cmd)
    subprocess.call(cmd, shell = True)
    
    #filter_low_conf
    if(args.verbose):
        logger.info("run_wes-filters: Applying filter_low_conf")
    tempMaf2 = os.path.join(tmpdir,"maf2.maf")
    cmd =  apply_filter + " filter_low_conf.R " + tempMaf1  + " " + tempMaf2
    if(args.verbose):
        logger.info("run_wes-filters: Running, %s",cmd)
    subprocess.call(cmd, shell = True)
    
    #filter_ffpe
    if(args.verbose):
        logger.info("run_wes-filters: Applying filter_ffpe")
    tempMaf3 = os.path.join(tmpdir,"maf3.maf")
    cmd =  apply_filter + " filter_ffpe.R " + tempMaf2  + " " + tempMaf3 
    if(args.verbose):
        logger.info("run_wes-filters: Running, %s",cmd)
    subprocess.call(cmd, shell = True)
    
    #filter_dmp
    if(args.verbose):
        logger.info("run_wes-filters: Applying filter_dmp")
    tempMaf4 = os.path.join(tmpdir,"maf4.maf")
    cmd =  apply_filter + " filter_dmp.R " + tempMaf3  + " " + tempMaf4
    if(args.verbose):
        logger.info("run_wes-filters: Running, %s",cmd)
    subprocess.call(cmd, shell = True)
    
    #filter_ffpe_pool
    tempMaf5 = os.path.join(tmpdir,"maf5.maf")
    if(args.FFPEPoolMaf):
        if(args.verbose):
            logger.info("run_wes-filters: Applying filter_ffpe_pool")
        md =  apply_filter + " filter_ffpe_pool.R " + tempMaf4  + " " + tempMaf5 + " -f " + args.FFPEPoolMaf + " -fo 1"
        if(args.verbose):
            logger.info("run_wes-filters: Running, %s",cmd)
        subprocess.call(cmd, shell = True)
    if(os.path.isfile(tempMaf4)):
        tempMaf5 = tempMaf5
    else:
        tempMaf5 = tempMaf4
    
    #filter_normal_panel
    tempMaf6 = os.path.join(tmpdir,"maf6.maf")
    if(args.NormalPanelMaf):
        if(args.verbose):
            logger.info("run_wes-filters: Applying filter_normal_panel")
        md =  apply_filter + " filter_normal_panel.R " + tempMaf5  + " " + tempMaf6 + " -f " + args.NormalPanelMaf + " -fo 1"
        if(args.verbose):
            logger.info("run_wes-filters: Running, %s",cmd)
        subprocess.call(cmd, shell = True)
    if(os.path.isfile(tempMaf6)):
        tempMaf6 = tempMaf6
    else:
        tempMaf6 = tempMaf5
    
    #filter_cohort_normals
    tempMaf7 = os.path.join(tmpdir,"maf7.maf")
    if(args.NormalCohortMaf):
        if(args.verbose):
            logger.info("run_wes-filters: Applying filter_cohort_normals")
        md =  apply_filter + " filter_cohort_normals.R " + tempMaf5  + " " + tempMaf6 + " -f " + args.NormalCohortMaf
        if(args.verbose):
            logger.info("run_wes-filters: Running, %s",cmd)
        subprocess.call(cmd, shell = True)
    if(os.path.isfile(tempMaf6)):
        tempMaf7 = tempMaf7
    else:
        tempMaf7 = tempMaf6
    
    #move to output name
    if(args.verbose):
        logger.info("run_wes-filters: Moving and cleaning temp directory")
    if(args.outdir):
        outfile = os.path.join(args.outdir + args.outputMaf)
    else:
        outfile = os.path.join(os.getcwd() + args.outputMaf)
    cmd = "mv" + tempMaf7 + " " + outfile
    if(args.verbose):
        logger.info("run_wes-filters: Running, %s",cmd)
    subprocess.call(cmd, shell = True)
    subprocess.call("rm -rf " + tmpdir, shell = True)
    return(outfile)

    
if __name__ == "__main__":
    start_time = time.time()  
    main()
    end_time = time.time()
    totaltime = end_time - start_time
    logging.info("run_wes-filters: Elapsed time was %g seconds", totaltime)
    sys.exit(0)
