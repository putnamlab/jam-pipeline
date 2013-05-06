
# (c) 2012-2013 Rice University & Nicholas H. Putnam
#
# This file is part of jam-pipeline
#
# This work is licensed under the Creative Commons Attribution 3.0 Unported License. To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/.


import unittest
import os
#import random
import subprocess
import hashlib 
import re
import glob


import gzip 

class TestJamMmScan(unittest.TestCase):

    def setUp(self):
        self.jam_root = os.environ.get('JAM_ROOT')
        self.cksum={}
        f = open(os.path.join( self.jam_root,"checksums.txt"))
        while True:
            l=f.readline()
            if not l:
                break
            c=l.strip().split()
            self.cksum[c[0]]=c[1]
        f.close()

#    def x_test_make_bv_commands(self):
#        wc = subprocess.call(["kmerpipe.py","--drop","--gspec","Ltest","--debug","--force","--serial"])
        
    def test_mmscan(self):

        seqdir = os.path.join( self.jam_root,"projects","Limulus_testpolyphemus", "sequence","FastaMasked" )
        files = glob.glob(seqdir+"/*.fam.gz")
        print "files:",seqdir,files
        kmers_dir = os.path.join( self.jam_root,"projects","Limulus_testpolyphemus", os.environ['JAM_ANALYSIS_DIR'], "kmers")
        mmscan_dir = os.path.join( self.jam_root,"projects","Limulus_testpolyphemus", os.environ['JAM_ANALYSIS_DIR'], "mmscan")

        try:
            os.mkdir(mmscan_dir)
        except:
            print "MmScan dir already there"

        for f in files:
            print f
            m = re.search("([^\/]*\.fam\.gz)",f)
            of = re.sub("fam.gz","Mmscan.gz",m.group(1))
            of = os.path.join( mmscan_dir, of)
            
            cmd = "GenomeMmScan -o 23  -i <( cat %s/MmTable.11slice5.txt %s/snpmers-filt.txt ) -a -H 100000 -s -S 11:5 %s 2> /dev/null | gzip > %s" % (kmers_dir,kmers_dir,f,of)
            print cmd
            subprocess.call(["bash","-c",cmd])


        files = glob.glob(mmscan_dir+"/*f.Mmscan.gz")
        for f in files:
            cmd= "matchpairs.pl %s" % (f)
            print cmd
            subprocess.call(["bash","-c",cmd])
            
        files = glob.glob(mmscan_dir+"/*mated.Mmscan.gz")
        for f in files:
            outfile = re.sub("Mmscan.gz$","map3.gz",f)
            cmd = "GenomeLinkContigs  -d oasi -o 23 -c %s/contigs.txt -n 90000 -H 1700000 %s 2> /dev/null | gzip > %s" % (kmers_dir,f,outfile)
            print cmd
            subprocess.call(["bash","-c",cmd])
                

        return(True)

if __name__ == '__main__':
    #unittest.main()
    suite = unittest.TestLoader().loadTestsFromTestCase(TestJamMmScan)
#    unittest.TextTestRunner(verbosity=2).run(suite)
    r=unittest.TextTestRunner(verbosity=2).run(suite)
    if not r.wasSuccessful():
        exit(1)
    else:
        exit(0)
