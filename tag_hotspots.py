#!/usr/bin/python
'''
@Description : This tool helps to tag hotspot events
@Created :  05/10/2017
@Updated : 05/10/2017
@author : Ronak H Shah

'''
from __future__ import division
import argparse
import sys
import os
import time
import logging

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
try:
    import pandas as pd
except ImportError:
    logger.fatal("tag_hotspots: pandas is not installed, please install pandas as it is required to run the process.")
    sys.exit(1)

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
   (mafDF) =read_maf(args)
   (hotspotDF) =read_hotspots(args)
   (taggedDF) = tag_hotspots(args, mafDF, hotspotDF)
   write_output(args,taggedDF)
   if(args.verbose):
       logger.info("tag_hotspots: Finished the run for tagging hotspots.")

def read_maf(args):
    dataDF = pd.read_table(args.inputMaf, comment="#", low_memory=False)
    return(dataDF)

def read_hotspots(args):
    dataDF = pd.read_table(args.inputTxt, comment="#", low_memory=False)
    return(dataDF)


def tag_hotspots(args,mafDF,hotspotsDF):
    mafDF_copy = mafDF.copy()
    mafDF_copy["hotspot_whitelist"] = None
    for i_index, i_row in mafDF.iterrows():
        m_chr = i_row.loc['Chromosome']
        m_start = i_row.loc['Start_Position']
        m_end = i_row.loc['End_Position']
        m_type = i_row.loc['TYPE']
        m_vt = i_row.loc['Variant_Type']
        m_ref = (str(i_row.loc['Reference_Allele'])).rstrip()
        m_alt = (str(i_row.loc['Tumor_Seq_Allele2'])).rstrip()
        mafDF_copy.set_value(i_index,"hotspot_whitelist","FALSE")
        is_SNV = False
        if(len(m_ref) == len(m_alt)):
            is_SNV = True
        m_hgvs_p_short = (str(i_row.loc['HGVSp_Short'])).rstrip()
        hotspot_subDF = hotspotsDF.loc[hotspotsDF['Chromosome'] == m_chr]
        for j_index,j_row in hotspot_subDF.iterrows():
            j_chr = j_row.loc['Chromosome']
            j_start = j_row.loc['Start_Position']
            j_end = j_row.loc['End_Position']
            j_hgvs_p_short = (str(i_row.loc['HGVSp_Short'])).rstrip()
            if(is_SNV):
                if(j_start == m_start and j_hgvs_p_short == m_hgvs_p_short):
                    mafDF_copy.set_value(i_index,"hotspot_whitelist","TRUE")
                    break
                else:
                    continue
            else:
                if(j_start == m_start):
                    mafDF_copy.set_value(i_index,"hotspot_whitelist","TRUE")
                    break
                else:
                    continue
    return(mafDF_copy)

def write_output(args,output_DF):
    if(args.outdir):
        outFile = os.path.join(args.outdir,args.outputMaf)
    else:
        outFile = os.path.join(os.getcwd(),args.outputMaf)
    output_DF.to_csv(outFile, sep='\t', index=False)
    
if __name__ == "__main__":
    start_time = time.time()  
    main()
    end_time = time.time()
    totaltime = end_time - start_time
    logging.info("tag_hotspots: Elapsed time was %g seconds", totaltime)
    sys.exit(0)
