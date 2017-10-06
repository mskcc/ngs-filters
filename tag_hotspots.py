#!/usr/bin/python
'''
@Description : This tool helps to tag hotspot events
@Created :  05/10/2017
@Updated : 10/06/2017
@author : Ronak H Shah, Cyriac Kandoth

'''
from __future__ import division
import argparse, sys, os, time, logging, csv, re

logging.basicConfig(
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=logging.DEBUG)
logger = logging.getLogger('tag_hotspots')
try:
    import coloredlogs
    coloredlogs.install(level='DEBUG')
except ImportError:
    logger.warning("tag_hotspots: coloredlogs is not installed, please install it if you wish to see color in logs on standard out.")
    pass

def main():
    parser = argparse.ArgumentParser(prog='tag_hotspots.py', description=' This tool helps to tag hotspot events', usage='%(prog)s [options]')
    parser.add_argument("-v", "--verbose", action="store_true", dest="verbose", help="make lots of noise")
    parser.add_argument("-m", "--input-maf", action="store", dest="inputMaf", required=True, type=str, metavar='SomeID.maf', help="Input maf file which needs to be tagged")
    parser.add_argument("-itxt", "--input-hotspot", action="store", dest="inputTxt", required=True, type=str, metavar='SomeID.txt', help="Input txt file which has hotspots")
    parser.add_argument("-o","--output-maf", action="store", dest="outputMaf", required=True, type=str, metavar='SomeID.maf', help="Output maf file name")
    parser.add_argument("-outdir", "--outDir", action="store", dest="outdir", required=False, type=str, metavar='/somepath/output', help="Full Path to the output dir.")

    args = parser.parse_args()
    if(args.verbose):
        logger.info("tag_hotspots: Started the run for tagging hotspots")

    # Load hotspots into a dict for easy lookup by chr:pos:ref:alt, and store AA position changed
    hotspot = dict()
    with open(args.inputTxt, 'rb') as infile:
        reader = csv.DictReader(infile, delimiter='\t')
        for row in reader:
            key = ':'.join([row['Chromosome'], row['Start_Position'], row['Reference_Allele'], row['Tumor_Seq_Allele2']])
            # Extract the amino-acid position from the HGVSp code, to store in the dict
            aa_pos = re.match( r'^p\.\D+(\d+)', row['HGVSp_Short'])
            if aa_pos:
                hotspot[key] = aa_pos.group(1)

    # Parse through input MAF, and create a new one with an extra column tagging hotspots
    with open(args.inputMaf, 'rb') as infile, open(args.outputMaf, 'wb') as outfile:
        # ::NOTE:: Comment lines are tossed, though they may need to be retained in some use cases
        reader = csv.DictReader((row for row in infile if not row.startswith('#')), delimiter='\t')
        writer = csv.DictWriter(outfile, delimiter='\t', lineterminator='\n', fieldnames=reader.fieldnames+["hotspot_whitelist"])
        writer.writeheader()
        for row in reader:
            row['hotspot_whitelist'] = "FALSE"
            key = ':'.join([row['Chromosome'], row['Start_Position'], row['Reference_Allele'], row['Tumor_Seq_Allele2']])
            if key in hotspot:
                row['hotspot_whitelist'] = "TRUE"
            writer.writerow(row)

    if(args.verbose):
        logger.info("tag_hotspots: Finished the run for tagging hotspots.")

if __name__ == "__main__":
    start_time = time.time()  
    main()
    end_time = time.time()
    totaltime = end_time - start_time
    logging.info("tag_hotspots: Elapsed time was %g seconds", totaltime)
    sys.exit(0)
