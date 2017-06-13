#!/usr/bin/python
'''
@Description : This tool helps to run all scripts in ngs-filter
@Created :  05/30/2017
@Updated : 06/01/2017
@author : Ronak H Shah

'''

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
    logger.warning("run_ngs-filters: coloredlogs is not installed, please install it if you wish to see color in logs on standard out.")
    pass

def main():
   parser = argparse.ArgumentParser(prog='run_ngs-filters.py', description=' This tool helps to tag hotspot events', usage='%(prog)s [options]')
   parser.add_argument("-v", "--verbose", action="store_true", dest="verbose", help="make lots of noise")
   parser.add_argument("-m", "--input-maf", action="store", dest="inputMaf", required=True, type=str, metavar='SomeID.maf', help="Input maf file which needs to be tagged")
   parser.add_argument("-o","--output-maf", action="store", dest="outputMaf", required=True, type=str, metavar='SomeID.maf', help="Output maf file name")
   parser.add_argument("-outdir", "--outDir", action="store", dest="outdir", required=False, type=str, metavar='/somepath/output', help="Full Path to the output dir.")
   parser.add_argument("-npmaf","--normal-panel-maf",action="store", dest="NormalPanelMaf", required=False, type=str, metavar='/somepath/to/normalpanel.maf', help="Path to fillout maf file of panel of standard normals")
   parser.add_argument("-fpmaf", "--ffpe_pool_maf", action="store", dest="FFPEPoolMaf", required=False, type=str, metavar='/somepath/to/ffpe_pool.maf', help="Path to fillout maf file for FFPE artifacts")
   parser.add_argument("-ncmaf","--normal-cohort-maf",action="store", dest="NormalCohortMaf", required=False, type=str, metavar='/somepath/to/normalcohort.maf', help="Path to fillout maf file of cohort normals")
   parser.add_argument('-nsf', '--normalSamplesFile', action="store", dest="NormalCohortSamples", required=False, type=str,metavar='/somepath/to/normalcohort.list', help='File with list of normal samples')
   parser.add_argument("-hsp", "--input-hotspot", action="store", dest="inputHSP", required=False, type=str, metavar='SomeID.txt', help="Input txt file which has hotspots")
   
   args = parser.parse_args()
   this_dir, this_filename = os.path.split(__file__)
   if(args.NormalCohortMaf):
       if(args.NormalCohortSamples):
           pass
       else:
           logger.critical("run_ngs-filters: Value for --normal-cohort-maf is give but --normalSamplesFile is empty")
           sys.exit(1)
   else:
        pass
   if(args.inputHSP):
       pass
   else:
       args.inputHSP = os.path.join(this_dir,"data","hotspot-list-union-v1-v2.txt")
   if(args.verbose):
       logger.info("run_ngs-filters: Started the run for wes-filters")
   (finalmaf) = run_wes_filters(args)
   if(args.verbose):
       logger.info("run_ngs-filters: Output is written in %s", finalmaf)
       logger.info("run_ngs-filters: Finished the run for wes-filters.")


def run_wes_filters(args):
    tmpdir = tempfile.mkdtemp()
    wes_filter_bin = os.path.dirname(os.path.realpath(__file__))
    tag_hotspot = os.path.join(wes_filter_bin,"tag_hotspots.py")
    apply_filter = os.path.join(wes_filter_bin,"applyFilter.sh")
    
    #tag_hotspots
    if(args.verbose):
        logger.info("run_ngs-filters: Tagging Hotspots")
    tempMaf0 = os.path.join(tmpdir,"maf0.maf")
    cmd =  "python " + tag_hotspot + " -v -m " + args.inputMaf  + " -itxt " + args.inputHSP + " -o " + tempMaf0
    if(args.verbose):
        logger.info("run_ngs-filters: Running, %s",cmd)
    subprocess.call(cmd, shell = True)
    
    #black_list_region
    if(args.verbose):
        logger.info("run_ngs-filters: Applying filter_blacklist_regions")
    tempMaf1 = os.path.join(tmpdir,"maf1.maf")
    cmd =  apply_filter + " " + os.path.join(wes_filter_bin,"filter_blacklist_regions.R") + " " + tempMaf0  + " " + tempMaf1
    if(args.verbose):
        logger.info("run_ngs-filters: Running, %s",cmd)
    subprocess.call(cmd, shell = True)
    
    #filter_low_conf
    if(args.verbose):
        logger.info("run_ngs-filters: Applying filter_low_conf")
    tempMaf2 = os.path.join(tmpdir,"maf2.maf")
    cmd =  apply_filter + " " + os.path.join(wes_filter_bin, "filter_low_conf.R") + " " + tempMaf1  + " " + tempMaf2
    if(args.verbose):
        logger.info("run_ngs-filters: Running, %s",cmd)
    subprocess.call(cmd, shell = True)
    
    #filter_ffpe
    if(args.verbose):
        logger.info("run_ngs-filters: Applying filter_ffpe")
    tempMaf3 = os.path.join(tmpdir,"maf3.maf")
    cmd =  apply_filter + " " + os.path.join(wes_filter_bin,"filter_ffpe.R") + " " + tempMaf2  + " " + tempMaf3 
    if(args.verbose):
        logger.info("run_ngs-filters: Running, %s",cmd)
    subprocess.call(cmd, shell = True)
    
    #filter_dmp
    if(args.verbose):
        logger.info("run_ngs-filters: Applying filter_dmp")
    tempMaf4 = os.path.join(tmpdir,"maf4.maf")
    cmd =  apply_filter + " " + os.path.join(wes_filter_bin,"filter_dmp.R") + " " + tempMaf3  + " " + tempMaf4
    if(args.verbose):
        logger.info("run_ngs-filters: Running, %s",cmd)
    subprocess.call(cmd, shell = True)
    
    #filter_ffpe_pool
    tempMaf5 = os.path.join(tmpdir,"maf5.maf")
    if(args.FFPEPoolMaf):
        if(args.verbose):
            logger.info("run_ngs-filters: Applying filter_ffpe_pool")
        cmd =  apply_filter + " " + os.path.join(wes_filter_bin,"filter_ffpe_pool.R") + " " + tempMaf4  + " " + tempMaf5 + " -f " + args.FFPEPoolMaf + " -fo 1"
        if(args.verbose):
            logger.info("run_ngs-filters: Running, %s",cmd)
        subprocess.call(cmd, shell = True)
    if(os.path.isfile(tempMaf5)):
        tempMaf5 = tempMaf5
    else:
        tempMaf5 = tempMaf4
    
    #filter_normal_panel
    tempMaf6 = os.path.join(tmpdir,"maf6.maf")
    if(args.NormalPanelMaf):
        if(args.verbose):
            logger.info("run_ngs-filters: Applying filter_normal_panel")
        cmd =  apply_filter + " " + os.path.join(wes_filter_bin,"filter_normal_panel.R") + " " + tempMaf5  + " " + tempMaf6 + " -f " + args.NormalPanelMaf + " -fo 1"
        if(args.verbose):
            logger.info("run_ngs-filters: Running, %s",cmd)
        subprocess.call(cmd, shell = True)
    if(os.path.isfile(tempMaf6)):
        tempMaf6 = tempMaf6
    else:
        tempMaf6 = tempMaf5
    
    #filter_cohort_normals
    tempMaf7 = os.path.join(tmpdir,"maf7.maf")
    if(args.NormalCohortMaf):
        if(args.verbose):
            logger.info("run_ngs-filters: Applying filter_cohort_normals")
        cmd =  apply_filter + " " + os.path.join(wes_filter_bin,"filter_cohort_normals.R") + " " + tempMaf5  + " " + tempMaf6 + " -f " + args.NormalCohortMaf + " -N " + args.NormalCohortSamples
        if(args.verbose):
            logger.info("run_ngs-filters: Running, %s",cmd)
        subprocess.call(cmd, shell = True)
    if(os.path.isfile(tempMaf7)):
        tempMaf7 = tempMaf7
    else:
        tempMaf7 = tempMaf6
    
    #move to output name
    if(args.verbose):
        logger.info("run_ngs-filters: Moving and cleaning temp directory")
    if(args.outdir):
        outfile = os.path.join(args.outdir , args.outputMaf)
    else:
        outfile = os.path.join(os.getcwd() , args.outputMaf)
    cmd = "mv " + tempMaf7 + " " + outfile
    if(args.verbose):
        logger.info("run_ngs-filters: Running, %s",cmd)
    subprocess.call(cmd, shell = True)
    if(args.verbose):
        logger.info("run_ngs-filters: Running, rm -rf %s",tmpdir)
    subprocess.call("rm -rf " + tmpdir, shell = True)
    return(outfile)

    
if __name__ == "__main__":
    start_time = time.time()  
    main()
    end_time = time.time()
    totaltime = end_time - start_time
    logging.info("run_ngs-filters: Elapsed time was %g seconds", totaltime)
    sys.exit(0)
