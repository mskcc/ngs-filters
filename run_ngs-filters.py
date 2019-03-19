#!/usr/bin/python
'''
@Description : This tool helps to run all scripts in ngs-filter
@Created :  05/30/2017
@Updated : 10/11/2017
@author : Ronak H Shah, Cyriac Kandoth

'''

import argparse
import sys
import os
import shutil
import time
import logging
import tempfile
import subprocess
import getpass

logging.basicConfig(
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    datefmt='%m/%d/%Y %I:%M:%S %p',
    level=logging.DEBUG)
logger = logging.getLogger('run_wes-filters')


def main():
    parser = argparse.ArgumentParser(prog='run_ngs-filters.py', description=' This tool helps to tag hotspot events', usage='%(prog)s [options]')
    parser.add_argument("-v", "--verbose", action="store_true", dest="verbose", help="make lots of noise")
    parser.add_argument("-m", "--input-maf", action="store", dest="inputMaf", required=True, type=str, metavar='SomeID.maf', help="Input maf file which needs to be tagged")
    parser.add_argument("-o", "--output-maf", action="store", dest="outputMaf", required=True, type=str, metavar='SomeID.maf', help="Output maf file")
    parser.add_argument("-npmaf", "--normal-panel-maf", action="store", dest="NormalPanelMaf", required=False, type=str, metavar='/somepath/to/normalpanel.maf', help="Path to fillout maf file of panel of standard normals")
    parser.add_argument("-hsp", "--input-hotspot", action="store", dest="inputHSP", required=False, type=str, metavar='SomeID.txt', help="Input txt file which has hotspots")
    args = parser.parse_args()
    this_dir, this_filename = os.path.split(__file__)

    if(args.inputHSP):
        pass
    else:
        args.inputHSP = os.path.join(this_dir, "data", "hotspot-list-union-v1-v2.txt")

    if(args.verbose):
        logger.info("run_ngs-filters: Started the run for wes-filters")
    (finalmaf) = run_wes_filters(args)
    if(args.verbose):
        logger.info("run_ngs-filters: Output is written in %s", finalmaf)
        logger.info("run_ngs-filters: Finished the run for wes-filters.")


def run_wes_filters(args):
    wes_filter_bin = os.path.dirname(os.path.realpath(__file__))
    tag_filters = os.path.join(wes_filter_bin, "tag_filters.py")
    apply_filter = os.path.join(wes_filter_bin, "applyFilter.sh")

    # Create two MAFs that we'll alternate between, as input/output to each filter
    tmpdir = tempfile.mkdtemp()
    tempMaf0 = os.path.join(tmpdir, "maf0.maf")
    tempMaf1 = os.path.join(tmpdir, "maf1.maf")

    #tag_filters
    if(args.verbose):
        logger.info("run_ngs-filters: Tagging Hotspots, Applying nad3 filter, Fix MuTect artifacts")
    cmd = "python " + tag_filters + " -v -m " + args.inputMaf + " -itxt " + args.inputHSP + " -o " + tempMaf0
    if(args.verbose):
        logger.info("run_ngs-filters: Running, %s", cmd)
    subprocess.call(cmd, shell=True)

    #black_list_region
    if(args.verbose):
        logger.info("run_ngs-filters: Applying filter_blacklist_regions")
    cmd = apply_filter + " " + os.path.join(wes_filter_bin, "filter_blacklist_regions.R") + " " + tempMaf0 + " " + tempMaf1
    if(args.verbose):
        logger.info("run_ngs-filters: Running, %s", cmd)
    subprocess.call(cmd, shell=True)
    os.remove(tempMaf0)

    #filter_normal_panel
    if(args.NormalPanelMaf):
        if(args.verbose):
            logger.info("run_ngs-filters: Applying filter_normal_panel")
        cmd = apply_filter + " " + os.path.join(wes_filter_bin, "filter_normal_panel.R") + " " + tempMaf1 + " " + tempMaf0 + " -f " + args.NormalPanelMaf + " -fo 1 "
        if(args.verbose):
            logger.info("run_ngs-filters: Running, %s", cmd)
        subprocess.call(cmd, shell=True)
        os.remove(tempMaf1)

    # Move final MAF to user-defined final output location
    if(args.verbose):
        logger.info("run_ngs-filters: Moving output-maf and cleaning temp directory")
    #General case to make sure the file exists before removing. Assumes that we clean up correctly in prior steps
    if os.path.exists(tempMaf1) and os.path.exists(tempMaf0):
        sys.exit('Error: Both temp maf files exist. Probably because of race condition.')
    elif os.path.exists(tempMaf1):
        shutil.move(tempMaf1, args.outputMaf)
    elif os.path.exists(tempMaf0):
        shutil.move(tempMaf0, args.outputMaf)
    else:
        sys.exit('Error: Temp maf files were malformed')
    shutil.rmtree(tmpdir)
    return(args.outputMaf)


if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    totaltime = end_time - start_time
    logging.info("run_ngs-filters: Elapsed time was %g seconds", totaltime)
    sys.exit(0)
