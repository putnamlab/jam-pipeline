import unittest
import os
#import random
import subprocess
import hashlib 



import gzip 

class TestJamTrim(unittest.TestCase):

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
        wc = subprocess.call(["kmerpipe.py","--drop","--gspec","Ltest","--debug","--force","--serial"])
        
    def test_trim(self):
        cf=open("Ltest.cmds")
        while True:
            l = cf.readline()
            if not l:
                break
            c=l.strip().split()
            print l.strip()
            wc = subprocess.call(["bash","-c",l.strip()])
            outfilename = c[-1].split(os.sep)[-1]
            fh = gzip.open(c[-1],"rb")
            d = fh.read()
            fh.close()
            print outfilename ,self.cksum.get(outfilename),hashlib.sha1(d).hexdigest()
            self.assertEqual(self.cksum.get(outfilename),hashlib.sha1(d).hexdigest())


if __name__ == '__main__':
    #unittest.main()
    suite = unittest.TestLoader().loadTestsFromTestCase(TestJamTrim)
#    unittest.TextTestRunner(verbosity=2).run(suite)
    r=unittest.TextTestRunner(verbosity=2).run(suite)
    if not r.wasSuccessful():
        exit(1)
    else:
        exit(0)
