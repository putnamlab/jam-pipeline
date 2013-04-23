import unittest
import os
#import random
import subprocess
import hashlib 
import re



import gzip 

class TestJamSNPmers(unittest.TestCase):

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

    def test_run_snpmer_pairing(self):
#        wc = subprocess.call(["kmerpipe.py","--drop","--gspec","Ltest","--debug","--force","--serial"])
        f=open("gbv_commands.txt")
        lines = f.readlines()
        f.close()
        bvfiles=[]
        kmers_dir=""
        for l in lines:
            m=re.search("> (\S+) ",l)
            bvfiles.append(m.group(1))
            #print m.group(1)
        m=re.search("^(.*)\/[^/]+",bvfiles[0])
        kmers_dir = m.group(1)
        cmd = "cat " + " ".join(bvfiles) + " | GenomeMmTable -o 23 -H 100000000 > %s/snpmers.txt 2> %s/snpmers.err" % (kmers_dir,kmers_dir)
        print cmd
        wc = subprocess.call(["bash","-c",cmd]) 

        cmd = "cat %s/snpmers.txt | awk '$5==3 || $5==21 || $5==12' > %s/snpmers-filt.txt"  % (kmers_dir,kmers_dir)
        print cmd
        wc = subprocess.call(["bash","-c",cmd]) 

        fh = open("%s/snpmers-filt.txt"%kmers_dir,"rb")
        d = fh.readlines()
        d.sort()
        fh.close()
        cksum = hashlib.sha1("".join(d)).hexdigest()
        print "snpmers-filt.txt",self.cksum.get("snpmers-filt.txt"),cksum
        self.assertEqual(self.cksum.get("snpmers-filt.txt"),cksum)



if __name__ == '__main__':
    #unittest.main()
    suite = unittest.TestLoader().loadTestsFromTestCase(TestJamSNPmers)
#    unittest.TextTestRunner(verbosity=2).run(suite)
    r=unittest.TextTestRunner(verbosity=2).run(suite)
    if not r.wasSuccessful():
        exit(1)
    else:
        exit(0)
