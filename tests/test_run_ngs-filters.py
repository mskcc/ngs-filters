'''
@Description : This tool helps test ngs-filters
@Created :  06/01/2017
@Updated : 03/27/2018
@author : Ronak H Shah, Cyriac Kandoth

'''

import filecmp
import os
import subprocess
from subprocess import Popen
import shlex
import nose
import logging, tempfile, shutil, sys

new_dir = tempfile.mkdtemp()
this_dir = os.path.realpath(__file__)
this_dir = "/".join(this_dir.split("/")[0:-2])

def setup_module(): 
    inputFileMaf = os.path.join(this_dir, "data", "sample_input.maf")
    outFileMaf = os.path.join(new_dir, "sample_output.maf")
    hotspotFile = os.path.join(this_dir, "data", "hotspot-list-union-v1-v2.txt")
    NormalPanelFill = os.path.join(this_dir, "data", "sample_input_fill.maf")
    scriptFile = os.path.join(this_dir, "run_ngs-filters.py")
    cmd = "python " + scriptFile + " -v -m " + inputFileMaf + " -o " + outFileMaf + " -npmaf " + NormalPanelFill + " -hsp " + hotspotFile
    args = shlex.split(cmd)
    if(os.path.isfile(outFileMaf)):
        os.remove(outFileMaf)
        
    try:
        proc = Popen(args)
        proc.wait()
        retcode = proc.returncode
        if(retcode >= 0):
            pass
    except:
        e = sys.exc_info()[0]
        logging.info("Running of python command: %s \n has failed. The exception produced is %s Thus we will exit",cmd,e)
        sys.exit(1)
             
def teardown_module():
     #shutil.rmtree(new_dir)
     pass

def test_maf_fileSimilarity():
    outFileMaf = os.path.join(new_dir, "sample_output.maf")
    cmpFileMaf = os.path.join(this_dir, "data", "sample_output.maf")
    
    cmd = "sed -i -e '/^#/ d' " + outFileMaf 
    subprocess.call(cmd, shell=True)
    nose.tools.ok_(filecmp.cmp(outFileMaf, cmpFileMaf), msg="The current result text file and the original result text file for run_ngs-filters are not the same") 

if __name__ == '__main__':
    nose.main()
