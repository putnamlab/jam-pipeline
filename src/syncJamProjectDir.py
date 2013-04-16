#!/usr/bin/env python
#import sys
#import getpass
import os
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-l","--login",dest="login",default="nhp1")
parser.add_option("-H","--host",dest="host",default="davinci")
parser.add_option("-r","--clusterRoot", dest="clusterRoot", default = "/scratch/nhp1/projects")
#parser.add_option("-R","--labRoot", dest="labRoot", default = "nputnam@10.68.3.163:/Volumes/Sequence/projects")
parser.add_option("-R","--labRoot", dest="labRoot", default = "/Volumes/Sequence/projects")
parser.add_option("-P","--project", dest="project" )
parser.add_option("-S","--skipSequence", dest="skipSequence", action="store_true" )
parser.add_option("-B","--skipBambus", dest="skipBambus", action="store_true" )
parser.add_option("-b","--backup", dest="backup", action="store_true" )
parser.add_option("-d","--debug", dest="debug", action="store_true" )

(options, args) = parser.parse_args()

ls = os.listdir( options.labRoot + "/" + options.project )
if options.debug:
    print ls

rsyncOptions = "-La --progress " + " ".join([ "--exclude '%s/tmp'" % (f) for f in ls ])


if options.skipSequence:
    rsyncOptions += " --exclude %s" % (   "sequence" ) 

for d in [ f for f in ls if f[:3] == "JAM" ]:
    if options.skipSequence:
        rsyncOptions += " --exclude %s/%s" % (  d, "scaffolding/3fBambus" ) 

if options.debug:
    print options,args



if options.backup:

    import subprocess
    if options.debug:
        print " ".join( ["ssh","%s@%s" % (options.login, options.host),"ls","-1DR", options.clusterRoot ]  )
    nw_proc = subprocess.Popen( ["ssh","%s@%s" % (options.login, options.host),"ls","-1DR", options.clusterRoot ] , stdout=subprocess.PIPE)
    listing = nw_proc.stdout.readlines()
    print listing
    exclusions = []
    for line in listing:
        if line[-2:]==":\n":
            dir = line[:-2]
#            print dir
#            print "x"
        elif line[-5:]=="tmp\n":
            pp = dir.find(options.project) + len(options.project)+1
            exclude = "/".join([dir[pp:],"tmp"])
#            print "##",dir[pp:],line, exclude
            exclusions.append(exclude)

    for ex in exclusions:
        rsyncOptions += " --exclude %s" % (  ex )
    
    rsyncCommand = "rsync %s %s@%s:%s/%s %s " % ( rsyncOptions, options.login, options.host, options.clusterRoot, options.project , options.labRoot)

else:
    rsyncCommand = "rsync %s %s/%s %s@%s:%s" % ( rsyncOptions, options.labRoot, options.project , options.login, options.host,options.clusterRoot )

print "###run this command:"
print rsyncCommand

#if not options.login:
#    options.login = getpass.getuser()








