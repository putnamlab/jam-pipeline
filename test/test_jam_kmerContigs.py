import unittest
import os
#import random
import subprocess
import hashlib 
import re



import gzip 

class TestJamContigs(unittest.TestCase):

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

    def test_kmer_edges(self):

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
#GenomeMmContigs  -H 17000 -w contigs.txt -E /Users/nputnam/JAMtestSequenceQuick/projects/Limulus_testpolyphemus/JAM-2013.04.29/kmers/edges.txt.gz -i <( cat /Users/nputnam/JAMtestSequenceQuick/projects/Limulus_testpolyphemus/JAM-2013.04.29/kmers/snpmers-filt.txt /Users/nputnam/JAMtestSequenceQuick/projects/Limulus_testpolyphemus/JAM-2013.04.29/kmers/MmTable.11slice5.txt ) -S 11:5
        
        cmd = "GenomeMmContigs -o 23 -i <( cat %s %s ) -S 11:5 -E %s/edges.txt.gz -H 170000 -w %s/contigs.txt" % (os.path.join(kmers_dir,"snpmers-filt.txt"),os.path.join(kmers_dir,"MmTable.11slice5.txt"),kmers_dir,kmers_dir) 
        print cmd
        wc = subprocess.call(["bash","-c",cmd]) 

#        cmd = "cat %s/snpmers.txt | awk '$5==3 || $5==21 || $5==12' > %s/snpmers-filt.txt"  % (kmers_dir,kmers_dir)
#        print cmd
#        wc = subprocess.call(["bash","-c",cmd]) 

        fh = open("%s/contigs.txt"%kmers_dir,"rb")
        d = fh.readlines()
        fh.close()
        cksum = hashlib.sha1("".join(d)).hexdigest()
        print "contigs.txt",self.cksum.get("contigs.txt"),cksum
        self.assertEqual(self.cksum.get("contigs.txt"),cksum)



if __name__ == '__main__':
    #unittest.main()
    suite = unittest.TestLoader().loadTestsFromTestCase(TestJamContigs)
#    unittest.TextTestRunner(verbosity=2).run(suite)
    r=unittest.TextTestRunner(verbosity=2).run(suite)
    if not r.wasSuccessful():
        exit(1)
    else:
        exit(0)
