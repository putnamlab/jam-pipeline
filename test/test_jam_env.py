import unittest
import os
import random

jam_executables=['kmerpipe.py','fastq2pfa.pl']
other_executables=['gzip','gzcat']

def which(program):

#    print program
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
#            print path
            path = path.strip('"')
            exe_file = os.path.join(path, program)
#            print exe_file
            if is_exe(exe_file):
                return exe_file

    return None

class TestJamEnv(unittest.TestCase):

    def setUp(self):
        self.jam_executables = jam_executables
        self.other_executables = other_executables

    def test_paths(self):
        for fn in self.jam_executables + self.other_executables :
            print fn
            fpath = which(fn)
#            print fpath
            self.assertTrue( os.path.isfile(fpath)  )
            self.assertTrue( os.access(fpath, os.X_OK) )

#        element = random.choice(self.seq)
#        self.assertTrue(element in self.seq)

if __name__ == '__main__':
    #unittest.main()
    suite = unittest.TestLoader().loadTestsFromTestCase(TestJamEnv)
    unittest.TextTestRunner(verbosity=2).run(suite)
