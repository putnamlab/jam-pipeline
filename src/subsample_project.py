#!/usr/bin/env python
from optparse import OptionParser
from fnmatch import fnmatch
import os, os.path
import md5
import re
import subprocess
import re
import gzip
#import Popen
#from path import path

desc="""
This script makes a copy of an existing JAM project directory by randomly selecting a subset of the templates, or by keeping only the templates in a list.
"""

parser = OptionParser(usage='Usage: %prog <options> project_directory new_project_directory', description=desc)
#parser.add_option("-d","--dir", dest="dir" )
#parser.add_option("-d", action="store_true" )
parser.add_option("-l","--label", dest="label" )
parser.add_option("-L","--list", dest="listfile" )
#parser.add_option("-D","--depth", dest="depth" )
parser.add_option("-F", dest="delim", default="/" )
parser.add_option("-s", "--samplePat", dest="samplePat", default="offspring_" )
parser.add_option("-S", "--seed", dest="seed", default="", help="random seed" )
parser.add_option("-n","--nothing", dest="nothing", action="store_true" )
parser.add_option("-G","--genomeSize", dest="genomeSize" )
parser.add_option("-r","--rate", dest="rate", type="float" )
parser.add_option("-t","--trimmed", dest="trimmed", action="store_true", help="Subsample the data under sequence/trim; default is to subsample the raw data.")
parser.add_option("-d","--debug", action="store_true")

(options, args) = parser.parse_args()

template_list={}
f=open(options.listfile)
while True:
    l = f.readline()
    if not l:
        break
    c=l.strip().split()
    template_list[c[0]]=1
f.close()

m=md5.new()

print options, args

print args[0], args[1]

sample_label_map = {"dad":"p1", "mom":"p2"}

read_lists={}
keep={}
digits=6
denom = (16.0**digits-1.0)

if not options.nothing:
#    new_dir=args[1]+"/sequence"
    try:
        print "mkdir",args[1]
        os.mkdir(args[1])
        os.mkdir(args[1]+"/sequence")
    except:
        print "couldn't make",args[1]

max_n_reads=0

def read_readnames(pattern,dir,files):
    global max_n_reads
#    print "mkdir",dir
    try:
        print "mkdir",args[1]+dir[len(args[0]):]
        os.mkdir( args[1]+dir[len(args[0]):])
    except:
        print "failed to make ",args[1]+dir[len(args[0]):]
        pass
    nreads=[0]

    for file in files:
#        print file, pattern
        m=re.search(pattern,file)
        if m:
            seqfile = dir+"/"+file
            print seqfile

            newfile = args[1]+"/"+dir[len(args[0]):]+"/"+file
            print "newfile:",newfile
            dirtokens=dir.split(options.delim)
            sample=dirtokens[-1]

            if options.debug:
                sample_label0 = re.sub(options.samplePat,"",sample)
                sample_label = sample_label_map.get(sample_label0,sample_label0)
                print dir, sample, sample_label, seqfile, newfile, max_n_reads
            if options.nothing:
                continue

            nreads=0
            f=False
            g=False
            if seqfile[-3:]==".gz":
                f=gzip.open(seqfile,"r")
                g=gzip.open(newfile,"w")
            else:
                f = open(seqfile)
                g = open(newfile,"w")
            phase=0

            fasta=False
            if re.search(".fasta$",seqfile):
                fasta=True

            while True:
                l = f.readline()
                if not l:
                    break
                if fasta and l[0]==">":
                    namem = re.match(">(\S*).[fr] ",l)
                    name = namem.group(1)
                    keep=False
                    if options.listfile:
                        if template_list.has_key(name):
                            keep=True
                    else:
                        x = int(md5.new(options.seed+name).hexdigest()[:digits],16)/denom
                        if x<options.rate:
                            keep=True
                    if keep:
                        g.write(l)
                        l=f.readline()
                        while l and l[0]!=">":
                            g.write(l)
                            l=f.readline()
                     
                elif phase==0:
#                    name = l[
                    name = l.split("/")[0]
                    namec = name.split()
                    keep=False

                    #print namec[0][1:]
                    
                    if options.listfile:
                        if template_list.has_key(namec[0][1:]):
                            keep=True
                    else:
                        x = int(md5.new(options.seed+name).hexdigest()[:digits],16)/denom
                        if x<options.rate:
#                        x = int(md5.new(options.seed+name).hexdigest()[:digits],16)/denom
#                        if x<options.rate:
                            keep=True
                    if keep:

                        g.write(l)
                        l = f.readline()
                        g.write(l)
                        l = f.readline()
                        g.write(l)
                        l = f.readline()
                        g.write(l)
#                        l = f.readline()
                        nreads+=1
                        phase=3
#                        if not l:
#                            break


                phase = (phase+1 ) %4
            f.close()
            g.close()
            if nreads > max_n_reads:
                max_n_reads = nreads
            if options.debug:
                for sample in read_lists.keys():
                    print sample,len(read_lists[sample].keys())


src_seq_dir = args[0]+"/sequence/raw"
if options.trimmed:
    src_seq_dir = args[0]+"/sequence/trim"
    
#os.path.walk(src_seq_dir,read_readnames, "(.fq$)|(.fastq$)|(.fasta$)|(.fastq.gz$)",followlinks=True)
files = os.walk(src_seq_dir,True,None,True) #,read_readnames, "(.fq$)|(.fastq$)|(.fasta$)|(.fastq.gz$)",followlinks=True)
for f in files:
    print f
    read_readnames("(.fq$)|(.fastq$)|(.fasta$)|(.fastq.gz$)",f[0],f[2])
    
for sample in read_lists.keys():
    print sample,len(read_lists[sample].keys())

print max_n_reads, "max reads per batch"




