#!/usr/bin/python
'''
@Description : This tool helps to tag hotspot events, apply nad3 filter, fix mutect artifacts
@Created :  05/10/2017
@Updated : 03/19/2019
@author : Ronak H Shah, Cyriac Kandoth, Tim Song

'''
from __future__ import division
import argparse, sys, os, time, logging, csv, re

logging.basicConfig(
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=logging.DEBUG)
logger = logging.getLogger('tag_filters')

def main():
    parser = argparse.ArgumentParser(prog='tag_filters.py', description=' This tool helps to tag hotspot events', usage='%(prog)s [options]')
    parser.add_argument("-v", "--verbose", action="store_true", dest="verbose", help="make lots of noise")
    parser.add_argument("-m", "--input-maf", action="store", dest="inputMaf", required=True, type=str, metavar='SomeID.maf', help="Input maf file which needs to be tagged")
    parser.add_argument("-itxt", "--input-hotspot", action="store", dest="inputTxt", required=True, type=str, metavar='SomeID.txt', help="Input txt file which has hotspots")
    parser.add_argument("-o","--output-maf", action="store", dest="outputMaf", required=True, type=str, metavar='SomeID.maf', help="Output maf file name")
    parser.add_argument("-outdir", "--outDir", action="store", dest="outdir", required=False, type=str, metavar='/somepath/output', help="Full Path to the output dir.")

    args = parser.parse_args()
    if(args.verbose):
        logger.info("tag_filters: Started the run for tagging Hotspots, applying nad3 filter, fixing MuTect artifacts")

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
        header = reader.fieldnames
        writer = csv.DictWriter(outfile, delimiter='\t', lineterminator='\n', fieldnames=reader.fieldnames+["hotspot_whitelist"])
        writer.writeheader()
        for row in reader:
            row['hotspot_whitelist'] = "FALSE"
            key = ':'.join([row['Chromosome'], row['Start_Position'], row['Reference_Allele'], row['Tumor_Seq_Allele2']])
            if key in hotspot:
                row['hotspot_whitelist'] = "TRUE"

            else: #Non hotspots
                if 'fillout_n_alt' in header: #Start tagging nad3
                    try:
                        if int(row['fillout_n_alt']) > 3:
                            if row['FILTER'] == 'PASS' or row['FILTER'] == '':
                                row['FILTER'] = 'nad3'
                            else:
                                row['FILTER'] = row['FILTER'] + ';nad3'
                    except ValueError: #Handles cases if fillout_n_alt is empty or string
                        pass
                if 'VSB' in header and 'set' in header and 'fillout_t_alt' in header:
                    if 'mutect' not in row['set'].lower() and str(row['VSB']) == '1' and int(row['fillout_t_alt']) > 10 and (int(row['fillout_t_forward_alt'])==0 or int(row['fillout_t_reverse_alt']) == 0 ):
                        if row['FILTER'] == 'PASS' or row['FILTER'] == '':
                            row['FILTER'] = 'asb'
                        else:
                            row['FILTER'] = row['FILTER'] + ';asb'

            if 'set' in header and 'FAILURE_REASON' in header:  # Fix artifact of bcftools annotate upstream
                if 'mutect' in row['set'].lower():
                    if row['FAILURE_REASON'] == '':
                        pass
                    elif row['FAILURE_REASON'] != 'None':
                        setlist = row['set'].split(',')
                        newsetlist = []
                        for caller in setlist:
                            if caller.lower() != 'mutect':
                                newsetlist.append(caller)
                        row['set'] = ','.join(newsetlist)


            writer.writerow(row)

    if(args.verbose):
        logger.info("tag_filters: Finished the run for tagging filters.")

if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    totaltime = end_time - start_time
    logging.info("tag_filters: Elapsed time was %g seconds", totaltime)
    sys.exit(0)
