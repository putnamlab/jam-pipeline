import unittest
import os
#import random
import subprocess
import hashlib 
import re



cksums={
"GenomeBVcount.17-0.out":"995b4b83a08aea9ea02d01460cf73c5ec399e9ef",
 "GenomeBVcount.17-1.out":"e442e9651fc10d81345bff5df1e31cb349ee025c",
 "GenomeBVcount.17-10.out":"2a6836072f91828ef09e3c2da2cecb211c69fd6e",
 "GenomeBVcount.17-11.out":"138648c94b55d8d2b39ac0d210c8ab17f2c96def",
 "GenomeBVcount.17-12.out":"56b3fc78800c148972289d90f033b7a2bafbd17c",
 "GenomeBVcount.17-13.out":"964984a807d2b223ffc9b058bf72863caa5bdf44",
 "GenomeBVcount.17-14.out":"a6cf56cc4703b8c4593da44a662f12caea176c12",
 "GenomeBVcount.17-15.out":"fd9efa05e6acc09c17190832a231d727f8c3141e",
 "GenomeBVcount.17-16.out":"0a606a83b21cc5126aecb8a5374870e2ba683067",
 "GenomeBVcount.17-2.out":"ae75598036aa46b570ea9198e3fa3d1776e99519",
 "GenomeBVcount.17-3.out":"52f6fa628ad431746d63875b841f21a72c4711c8",
 "GenomeBVcount.17-4.out":"e5b06c4eb74fdb4d8356afda2ce4799798fa4565",
 "GenomeBVcount.17-5.out":"78ac953c351226fddbeea19732df5d18b444f086",
 "GenomeBVcount.17-6.out":"1c59e4806aaf120a96c6e18f77286cd5e7919a5c",
 "GenomeBVcount.17-7.out":"110bc42cc253c19771813d9f647ea2d9cc57c1be",
 "GenomeBVcount.17-8.out":"a7c3c3d8f4f391bc99bd2c25cb51e75b164b89f4",
 "GenomeBVcount.17-9.out":"3f9431e739c098c3935b2e4d1d5cc7446d2c0ab4"
}


#" > /Users/nputnam/JAMtestSequence/projects/Limulus_testpolyphemus/JAM-2013.04.22/kmers/GenomeBVcount.17-0.out 

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

    def test_make_trim_commands(self):
#        wc = subprocess.call(["kmerpipe.py","--drop","--gspec","Ltest","--debug","--force","--serial"])
        cmd = "DriveGenomeBVcount.py --serial -C gbv_commands.txt %s/projects/Limulus_testpolyphemus"%(os.environ['JAM_ROOT'])
        wc = subprocess.call(cmd.split()) 
        
    def test_trim(self):
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
    unittest.TextTestRunner(verbosity=2).run(suite)
