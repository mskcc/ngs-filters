#!/usr/bin/python
'''
@Description : This tool helps to run all scripts in ngs-filter
@Created :  05/30/2017
@Updated : 10/11/2017
@author : Ronak H Shah, Cyriac Kandoth

'''

import argparse, sys, os, shutil, time, logging, tempfile, subprocess, getpass

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
    parser.add_argument("-npmaf", "--normal-panel-maf",action="store", dest="NormalPanelMaf", required=False, type=str, metavar='/somepath/to/normalpanel.maf', help="Path to fillout maf file of panel of standard normals")
    parser.add_argument("-ncmaf", "--normal-cohort-maf",action="store", dest="NormalCohortMaf", required=False, type=str, metavar='/somepath/to/normalcohort.maf', help="Path to fillout maf file of cohort normals")
    parser.add_argument('-nsf', "--normalSamplesFile", action="store", dest="NormalCohortSamples", required=False, type=str, metavar='/somepath/to/normalcohort.list', help='File with list of normal samples')
    parser.add_argument("-hsp", "--input-hotspot", action="store", dest="inputHSP", required=False, type=str, metavar='SomeID.txt', help="Input txt file which has hotspots")

    args = parser.parse_args()
    this_dir, this_filename = os.path.split(__file__)
    if(args.NormalCohortMaf):
        if(args.NormalCohortSamples):
            pass
        else:
            logger.critical("run_ngs-filters: Value for --normal-cohort-maf is given but --normalSamplesFile is empty")
            sys.exit(1)
    else:
         pass
    if(args.inputHSP):
        pass
    else:
        args.inputHSP = os.path.join(this_dir,"data","hotspot-list-union-v1-v2.txt")

    # Create a TMPDIR at /scratch/<username>, if it doesn't already exist
    tmp_root = "/scratch/" + getpass.getuser()
    if not os.path.exists(tmp_root):
        os.makedirs(tmp_root)
    os.environ["TMPDIR"] = tmp_root

    if(args.verbose):
        logger.info("run_ngs-filters: Started the run for wes-filters")
    (finalmaf) = run_wes_filters(args)
    if(args.verbose):
        logger.info("run_ngs-filters: Output is written in %s", finalmaf)
        logger.info("run_ngs-filters: Finished the run for wes-filters.")

def run_wes_filters(args):
    wes_filter_bin = os.path.dirname(os.path.realpath(__file__))
    tag_hotspot = os.path.join(wes_filter_bin,"tag_hotspots.py")
    apply_filter = os.path.join(wes_filter_bin,"applyFilter.sh")

    # Create two MAFs that we'll alternate between, as input/output to each filter
    tmpdir = tempfile.mkdtemp()
    tempMaf0 = os.path.join(tmpdir,"maf0.maf")
    tempMaf1 = os.path.join(tmpdir,"maf1.maf")

    #tag_hotspots
    if(args.verbose):
        logger.info("run_ngs-filters: Tagging Hotspots")
    cmd =  "python " + tag_hotspot + " -v -m " + args.inputMaf  + " -itxt " + args.inputHSP + " -o " + tempMaf0
    if(args.verbose):
        logger.info("run_ngs-filters: Running, %s",cmd)
    subprocess.call(cmd, shell = True)

    #black_list_region
    if(args.verbose):
        logger.info("run_ngs-filters: Applying filter_blacklist_regions")
    cmd =  apply_filter + " " + os.path.join(wes_filter_bin,"filter_blacklist_regions.R") + " " + tempMaf0  + " " + tempMaf1
    if(args.verbose):
        logger.info("run_ngs-filters: Running, %s",cmd)
    subprocess.call(cmd, shell = True)
    os.remove(tempMaf0)

    #filter_vaf
    if(args.verbose):
        logger.info("run_ngs-filters: Applying filter_vaf")
    cmd =  apply_filter + " " + os.path.join(wes_filter_bin,"filter_vaf.R") + " " + tempMaf1  + " " + tempMaf0
    if(args.verbose):
        logger.info("run_ngs-filters: Running, %s",cmd)
    subprocess.call(cmd, shell = True)
    os.remove(tempMaf1)

    #filter_normal_panel
    if(args.NormalPanelMaf):
        if(args.verbose):
            logger.info("run_ngs-filters: Applying filter_normal_panel")
        cmd =  apply_filter + " " + os.path.join(wes_filter_bin,"filter_normal_panel.R") + " " + tempMaf0  + " " + tempMaf1 + " -f " + args.NormalPanelMaf + " -fo 1 "
        if(args.verbose):
            logger.info("run_ngs-filters: Running, %s",cmd)
        subprocess.call(cmd, shell = True)
        os.remove(tempMaf0)
    else:
        # In case this filter was skipped, swap around the input/output MAFs
        tempMaf0, tempMaf1 = tempMaf1, tempMaf0

    #filter_cohort_normals
    if(args.NormalCohortMaf):
        if(args.verbose):
            logger.info("run_ngs-filters: Applying filter_cohort_normals")
        cmd =  apply_filter + " " + os.path.join(wes_filter_bin,"filter_cohort_normals.R") + " " + tempMaf1  + " " + tempMaf0 + " -f " + args.NormalCohortMaf + " -N " + args.NormalCohortSamples
        if(args.verbose):
            logger.info("run_ngs-filters: Running, %s",cmd)
        subprocess.call(cmd, shell = True)
        os.remove(tempMaf1)
    else:
        # In case this filter was skipped, swap around the input/output MAFs
        tempMaf0, tempMaf1 = tempMaf1, tempMaf0

    # Move final MAF to user-defined final output location
    if(args.verbose):
        logger.info("run_ngs-filters: Moving output-maf and cleaning temp directory")
    shutil.move(tempMaf0, args.outputMaf)
    shutil.rmtree(tmpdir)
    return(args.outputMaf)

if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    totaltime = end_time - start_time
    logging.info("run_ngs-filters: Elapsed time was %g seconds", totaltime)
    sys.exit(0)
