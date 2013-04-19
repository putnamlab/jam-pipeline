#!/usr/bin/env python

# Script to 

from optparse import OptionParser
import os
import sqlite3
import re
import subprocess

def getType(token):
    ivalue = None
    rvalue = None
    try:
        ivalue = int(token)
        rvalue = int(token)
        if (ivalue == rvalue):
            return "integer"
        else:
            return "real"
    except:
        return "text"

# Get the table field names from the first line of the file
# ... and scan all the lines of the file (admittedly inefficient)
#     to be sure whether it's numeric (make it "real" in the descrip)
#     or other (make it "text").
def describefields(file):
    # print os.getcwd()
    f = open(file)
    first = f.readline()
    first = re.sub('#', '', first)  # get rid of comment character on table descrip line
    first = first.strip()           # being paranoid about trailing newlines
    fieldids = first.split('\t')
    ftype = ["integer"]*len(fieldids) # Assume all numeric until proven otherwise

    for line in f:
        fieldvals = line.split('\t')
        for i in range(len(fieldvals)):
            temp = getType(fieldvals[i])
            if temp > ftype[i]:
                # take advantage of string comparison "text" > "real" > "integer"
                ftype[i] = temp
    f.close()
    descrips = [""] * len(fieldids)
    for i in range(len(descrips)):
        descrips[i] = fieldids[i] + " " + ftype[i]
        if descrips[i] == 'ID integer':
            descrips[i] = 'ID INTEGER PRIMARY KEY'
    return "(" + ", ".join(descrips) + ")"

def localOrNonlocalF(pathroot, path):
    if os.path.exists(path):
        return path
    elif os.path.exists(os.path.join(pathroot, path)):
        return os.path.join(pathroot, path)
    else:
        return None

def tablename(path):
    # - for now, just strip off leading path and trailing extension
    # - TODO: should double-check filename does't have bogus chars
    (head, tail) = os.path.split(path)
    (table, ext)  = os.path.splitext(tail)
    return table

def loadtable(conn, file):
    descrip = describefields(file)
    # print file, descrip
    table = tablename(file)
    # Should check here that the format of the table matches even if not reloading
    conn.cursor().execute('CREATE TABLE IF NOT EXISTS ' + table + ' ' + descrip)
    fd = open(file)
    for line in fd:
        # Skip comment lines and chromosomes with underscores
        if re.search('(^#)|(chr(\S+)_)', line):
            continue
        # Split on tab instead of whitespace, required for correctness
        # on fields that are empty or that contain spaces
        vals = (line.rstrip('\n')).split('\t')
        conn.cursor().execute("INSERT INTO " + table + " VALUES (" \
                  + "?,"*(len(vals)-1) + "?)", vals)
    fd.close()
    conn.commit()

def main():

    
    jam_root = os.environ.get('JAM_ROOT')
    if not jam_root:
        jam_root = "."

    usage="""kmerpipe.py [options]

for example:   kmerpipe.py --dbroot . --drop --projects ../projects/ --gspec Ltest

Creates sqlite databaset file sequencing.db

"""
    parser = OptionParser(usage=usage)
    parser.add_option("--dbroot", action="store", dest="dbroot",
                      default=jam_root+"/database", help="location of the database files.  Defaults to $JAM_ROOT/database")
    # print os.path.expanduser('~/Desktop/Dropbox/SeqPipe/Database')
    parser.add_option("--projects", action="store", dest="projects",
                      default=jam_root+"/projects", help="location of the projects dir.  Defaults to $JAM_ROOT/projects")
    parser.add_option("--clusterProjects", action="store", dest="clusterProjects",
                      default=False, help="location of the projects dir on the remote compute cluster (e.g. /scratch/nhp1/projects).  Defaults to the local projects root")
    parser.add_option("--incoming", action="store", dest="incoming",
                      default=jam_root+"/incoming", help="location of the incoming sequence dir.  Defaults to $JAM_ROOT/incoming")
    # Files to get table values for lanes and sequence files
    parser.add_option("--lanes",    action="store", dest="lanesF",
                      default='lanes.tsv', help="name of the lanes table.  defaults to lanes.tsv")
    parser.add_option("--seqfiles", action="store", dest="seqfilesF",
                      default='seqfiles.tsv', help="name of the table of sequence files.  defaults to seqfiles.tsv")
    parser.add_option("--species", action="store", dest="speciesF",
                      default='species.tsv', help="name of species table file.  defaults to species.tsv")
    parser.add_option("--gspec", action="store", dest="gspec",
                      default='Aluca', help="Species to process, in 5-letter mnemonic form, eg Lpoly for Limulus polyphemus")

    # Copy files up to the (davinci) cluster
    parser.add_option("--upload", action="store_true", dest="upload", default=False, help="Attempt to upload files to the computer cluster")

    parser.add_option("--serial", action="store_true", dest="serial", default=False, help="Generate serial commands rather than qsub commands")
    parser.add_option("--debug", action="store_true", dest="debug", default=False, help="Print debugging output")

    # Use option "--reload" to drop old tables, then create & load anew
    parser.add_option("--drop",   action="store_true", dest="drop", default=False, help="drop old tables, then create & load anew")
    # Force linking, clobbering old files or links
    parser.add_option("--force",  action="store_true", dest="force", default=False, help="force creation of links, even if they already exist.")

    (options, args) = parser.parse_args()
    # print "options:", options
    # print "args: ", args

    if not options.clusterProjects:
        options.clusterProjects = options.projects

    # Connect to database
    conn = sqlite3.connect(os.path.join(options.dbroot, 'sequencing.db'))
    c = conn.cursor()

    print options.dbroot
    # Read each file f into a table 
    for f in [localOrNonlocalF(options.dbroot, options.lanesF),
              localOrNonlocalF(options.dbroot, options.seqfilesF),
              localOrNonlocalF(options.dbroot, options.speciesF)]:
        if options.debug:
            print f
        if options.drop:
	    c.execute('DROP TABLE IF EXISTS ' + (tablename(f)))
        loadtable(conn, f)

    # Take in set of patterns for batches
    # Right now default is Gspec_ck like Aluca and group by LaneID
    c.execute("SELECT ReadsPerFile, Path, File, Gspec_ck, Organism, LaneID, Library, batch, ZeroQual, direction " +
              "FROM lanes, seqfiles " +
              "WHERE lanes.ID = LaneID " +
              "AND Gspec_ck = '%s'" % options.gspec)
    cmdfile = open(options.gspec + ".cmds", "w")
    for row in c:
        (readsPerFile, path, file,
         gspec, organism, laneid, library, batch, zq, direction) = row
        organism = organism.replace(" ", "_")
        newName = ""
        readRoot = ""
        batchsize = 999999999
        if int(readsPerFile):
            readRoot = "%s%d_%s" % (gspec, laneid, library)
            newName = "%s_%d_%s" % (readRoot, batch, direction)
            batchsize = int(readsPerFile)
        else:
            readRoot = "%s%d_%s" % (gspec, laneid, library)
            newName = "%s_%s" % (readRoot, direction)
        print newName
        newDir    = os.path.join(options.projects,organism,"sequence","raw") # "%s/sequence/raw" % organism
        maskedDir = os.path.join(options.projects,organism,"sequence","FastaMasked")

        newDirRemote    = os.path.join(options.clusterProjects,organism,"sequence","raw") # "%s/sequence/raw" % organism
        maskedDirRemote = os.path.join(options.clusterProjects,organism,"sequence","FastaMasked")

        if not os.path.isdir(newDir):
            subprocess.call("mkdir -p %s" % (newDir),shell=True)
        if not os.path.isdir(maskedDir):
            subprocess.call("mkdir -p %s" % (maskedDir),shell=True)

            #subprocess.call("ssh havlak@davinci.rice.edu mkdir -p /scratch/havlak/%s" % newDir,
            #                shell=True)
            #subprocess.call("ssh havlak@davinci.rice.edu mkdir -p /scratch/havlak/SeqPipe/%s/FastaMasked" % gspec,
            #                shell=True)
            #subprocess.call("ssh havlak@davinci.rice.edu mkdir -p /scratch/havlak/SeqPipe/%s/Kmers" % gspec,
            #                shell=True)
        subprocess.call("ln -s %s %s %s.fastq.gz"
                        % ("-f" if options.force else "",
                            os.path.join(options.incoming, path, file),  os.path.join( newDir, newName)),
                        shell=True)
        if options.serial:
            cmdfile.write('gzcat  %s.fastq.gz | fastq2pfa.pl -m 20 -soft 30 -bs %d -bnum %d -p %s -suf %s -zq %d | gzip > %s.fam.gz\n' %
                          ( os.path.join( newDir, newName), batchsize,batch, readRoot, direction, zq,  os.path.join( maskedDir, newName)))
        else:
            cmdfile.write(('echo "gunzip -c %s/%s.fastq.gz ' % (newDirRemote, newName))
                      + ('| fastq2pfa.pl -m 20 -soft 30 -bs %d ' % batchsize)
                      + ('-bnum %d -p %s -suf %s -zq %d ' % (batch, readRoot, direction, zq))
                      + ('| gzip > %s/%s.fam.gz" ' % ( maskedDirRemote, newName))
                      + ('| qsub -l cput=4:00:00 -d $PWD -N %s -m e -q serial -V\n' % newName) )
    cmdfile.close()
    if options.upload:
        cmd = "rsync -LPurt %s/SeqPipe havlak@davinci.rice.edu:/scratch/havlak/" % options.seqroot
        subprocess.call(cmd, shell=True)
    conn.commit()
    conn.close()

main()


