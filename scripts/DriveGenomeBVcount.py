#!/usr/bin/env python


from optparse import OptionParser
from fnmatch import fnmatch
import os, os.path
import md5
import re
import subprocess
from nhp_pbs_util import wait_for_jobs_to_finish
#import Popen
#from path import path


desc="""A driver script for GenomeBVcount, for use with JAM project directories on davinci.rice.edu.  By default, this submits jobs via qsub with options: -V -d projects_root -q serial -N batch.

"""

from datetime import date

print date.isoformat(date.today())
print date.today().strftime("%Y.%m.%d")
projects_root=os.environ.get('JAM_ROOT') # "/scratch/nhp1/projects/"

parser = OptionParser(usage='Usage: %prog <options> project_directory', description=desc)
#parser = optparse.OptionParser(usage='Usage: %prog <options>')
#parser.add_option("-d","--dir", dest="dir" )
#parser.add_option("-l","--label", dest="label" )
#parser.add_option("-D","--depth", dest="depth" )
parser.add_option("-k", dest="kmersize", default=23, type="int" , help="kmer size.  Default = 23")
parser.add_option("-f", dest="hashfact", default=17, type="int" , help="Hashing factor.  Default = 17")
parser.add_option("-H", dest="hashsize", default=10000000, type="int" , help="Hash size.  Default = 10000000")
parser.add_option("-C", dest="cmdfile", default=False , help="Instead of running, write commands to file FILE")
parser.add_option("-a", dest="analysis_dir", default="JAM-"+date.today().strftime("%Y.%m.%d"), type="string" , help="Output directory.  Default = JAM-%s"%(date.today().strftime("%Y.%m.%d")))
parser.add_option("-w","--wait", dest="wait", action="store_true" , help="Wait for the submitted pbs jobs to be done.")
parser.add_option("-s","--serial", dest="serial", action="store_true" , help="Don't submit to PBS.  Run in serial.")

#parser.add_option("-F", dest="delim", default="/" , help="File path delimeter.  Default = '/'")
#parser.add_option("-s", "--samplePat", dest="samplePat", default="offspring_", help="specify a pattern to strip off from the start of the sample names, e.g. 'offspring_'." )
parser.add_option("-n","--nothing", dest="nothing", action="store_true" , help="Dry run:  don't actually submit the trimming jobs.")
#parser.add_option("-b","--batch", dest="batch", type=int , help="Max number of reads per file.  Used for read renumbering")
#parser.add_option("-q","--soft", dest="soft", type=int , default=30, help="Softmask quality score threshold.  Default=30")
#parser.add_option("-Q","--hard", dest="hard", type=int , default=20, help="Hardmask quality score threshold.  Default=20")
#parser.add_option("-G","--genomeSize", dest="genomeSize" )
#parser.add_option("-r","--rate", dest="rate", type="float" )
parser.add_option("--debug", action="store_true")

(options, args) = parser.parse_args()
print options
analysis_dir = options.analysis_dir
m=md5.new()

print options, args

#print args[0], args[1]

sample_label_map = {"dad":"p1", "mom":"p2"}

#print "mkdir", args[1]
#print "mkdir", args[1]/sequence
#print "mkdir", args[1]/sequence/raw


read_lists={}
keep={}
digits=6
denom = (16.0**digits-1.0)

#if not options.nothing:
#    os.mkdir(args[1])

max_n_reads=0
files_to_count={}

import re

def collect_files(pattern,dir,files):

    for file in files:
        m=re.search(pattern,file)
        if m:
#            print file
            c = file.strip().split("_")
#            print c[0]
            files_to_count[c[1]] = files_to_count.get(c[1],[]) + [file]

project_dir = args[0]
#print project_dir+"/sequence/raw"
tmpdir= args[0]+"/sequence/FastaMasked/tmp"
try:
    if not options.nothing:
        os.mkdir(tmpdir)
except:
    print "%s already exists" % (tmpdir )


if options.debug:
    print project_dir+"/sequence/FastaMasked"
os.path.walk(project_dir+"/sequence/FastaMasked",collect_files, '\.fam.gz$')
print files_to_count

def ordering(x):
    if x[0]=="p":
        return((1,int(x[1:])))
    else:
        return((-1,int(x)))

def hash_to_args(h):
    samples = h.keys()
    samples.sort(key=ordering)
    r=" / ".join(map(str,[ " ".join([project_dir+"/sequence/FastaMasked/"+x for x in sorted(h[i]) ]) for i in samples ]))
    return r

filelist= hash_to_args(files_to_count)

outdir= project_dir + "/"+ analysis_dir
try:
    if not options.nothing:
        os.mkdir(outdir)
except:
    print "%s already exists" % ( args[0]+"/"+ analysis_dir )

outdir= project_dir+"/"+ analysis_dir+"/kmers"
try:
    if not options.nothing:
        os.mkdir(outdir)
except:
    print "%s already exists" % ( args[0]+"/"+ analysis_dir+"/kmers" )


cmdfile=False
if options.cmdfile:
    cmdfile = open(options.cmdfile,"wt")

#hashfact=17
#analysis_dir="JAM"
job_ids=[]
for i in range(options.hashfact):
    cmd = "GenomeBVcount -H %d -S %d:%d -d + -o %d %s > %s/GenomeBVcount.%d-%d.out 2> %s/GenomeBVcount.%d-%d.err"  % (options.hashsize, options.hashfact,i,options.kmersize,filelist,outdir,options.hashfact,i,outdir,options.hashfact,i ) #,analysis_dir,options.hashfact,i,analysis_dir,options.hashfact,i)
# ; cat ../../../%s/kmers/GenomeBVcount.%d-%d.out |  perl -ane 'print hex($F[1]); print \"\\n\"'  | perl ~/scripts/histogram2.pl - 1 1 > ../../../%s/kmers/GenomeBVcount.%d-%d.histogram.txt ; " 
#    cmd = "/home/havlak/bin/src/newGenomeMerHist/GenomeBVcount -H %d -S %d:%d -d + -o %d %s > ../../../%s/kmers/GenomeBVcount.%d-%d.out 2> ../../../%s/kmers/GenomeBVcount.%d-%d.err ; cat ../../../%s/kmers/GenomeBVcount.%d-%d.out |  perl -ane 'print hex($F[1]); print \"\\n\"'  | perl ~/scripts/histogram2.pl - 1 1 > ../../../%s/kmers/GenomeBVcount.%d-%d.histogram.txt ; " % (options.hashsize, options.hashfact,i,options.kmersize,filelist,analysis_dir,options.hashfact,i,analysis_dir,options.hashfact,i,analysis_dir,options.hashfact,i,analysis_dir,options.hashfact,i)
#    cmd = "cat ../../../%s/kmers/GenomeBVcount.%d-%d.out |  perl -ane 'print hex($F[1]); print \"\\n\"'  | perl ~/scripts/histogram2.pl - 1 1 > ../../../%s/kmers/GenomeBVcount.%d-%d.histogram.txt ; " % (analysis_dir,options.hashfact,i,analysis_dir,options.hashfact,i)
    if options.debug:
        print cmd # + " | qsub -d %s -V -N kmer%d -q serial" % (projects_root+tmpdir,i) 
    if not options.nothing:
        if options.serial:
            if options.cmdfile:
                cmdfile.write(cmd+"\n")
            else:
                subprocess.call(["bash","-c",cmd])
            
        else:
            wc=subprocess.Popen(["qsub","-V","-d",projects_root+tmpdir,"-q","serial","-N","kmer."+str(i)],stdin=subprocess.PIPE,stdout=subprocess.PIPE)
            out= wc.communicate(cmd)
            print "out",out
            job_name = out[0]
            tokens = job_name.split(".")
            print "tokens[0]",tokens[0]
            job_ids.append(tokens[0])

if cmdfile:
    cmdfile.close()

if not options.serial:
    print job_ids

    if options.wait:
        wait_for_jobs_to_finish(job_ids,10)

exit(1)


