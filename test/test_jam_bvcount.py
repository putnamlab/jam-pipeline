
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



import gzip 

class TestJamBV(unittest.TestCase):

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
        
    def test_bvcount(self):

        cmd = "DriveGenomeBVcount.py --serial -C gbv_commands.txt -a %s %s/projects/Limulus_testpolyphemus"%(os.environ['JAM_ANALYSIS_DIR'],os.environ['JAM_ROOT'])
        print cmd
        wc = subprocess.call(cmd.split()) 


        cf=open("gbv_commands.txt")
        while True:
            l = cf.readline()
            if not l:
                break
            m = re.search("(GenomeBVcount\.\d+-\d+\.out)",l)
            ofn = m.group(1)
            m = re.search(" > (.*GenomeBVcount\.\d+-\d+\.out)",l)
            ofn_fullpath = m.group(1)
            print ofn_fullpath
            c=l.strip().split()
            print ofn,
            wc = subprocess.call(["bash","-c",l.strip()])
#            outfilename = c[-1].split(os.sep)[-1]
            fh = open(ofn_fullpath,"rb")
            d = fh.readlines()
            d.sort()
            fh.close()
            cksum = hashlib.sha1("".join(d)).hexdigest()
            print ofn,self.cksum.get(ofn),cksum
            self.assertEqual(self.cksum.get(ofn),cksum)


if __name__ == '__main__':
    #unittest.main()
    suite = unittest.TestLoader().loadTestsFromTestCase(TestJamBV)
#    unittest.TextTestRunner(verbosity=2).run(suite)
    r=unittest.TextTestRunner(verbosity=2).run(suite)
    if not r.wasSuccessful():
        exit(1)
    else:
        exit(0)
