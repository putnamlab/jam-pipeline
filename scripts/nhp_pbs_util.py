import time
import subprocess
from xml.dom.minidom import parse, parseString
import re

wc=subprocess.Popen(["whoami"],stdout=subprocess.PIPE)
out=wc.communicate()
print "out",out
username=out[0].strip()
print username

def getText(nodelist):
    rc = []
    for node in nodelist:
        if node.nodeType == node.TEXT_NODE:
            rc.append(node.data)
    return ''.join(rc)


def job2hash(job):
    h={}
    for c in job.childNodes:
#        print c,c.nodeName,getText(c.childNodes)
        h[c.nodeName]=getText(c.childNodes)
#    print h
    return(h)

def n_jobs_for_user(username=username):
        print "count jobs for user",username
        n_still_running=0
        still_running=[]
        wc=subprocess.Popen(["qstat","-x"],stdout=subprocess.PIPE)
        out = wc.communicate()
        #print out                                                                                                                                                       
        if len(out[0])==0:
            print "no output"
            return(0)
        #print out
        dom= parseString(out[0])
        #print dom                                                                                                                                                       
    #    print dom.toprettyxml()                                                                                                                                         
        jobs=dom.getElementsByTagName("Job")
        
        for job in jobs:                                                                                                                                                
#            print dir(job)
            jh=job2hash(job)
#            print "jobid",jh['Job_Id'],jh['Job_Owner']

#            print jh.keys()
            #print jh['Job_Id'],jh['Job_Owner'],jh['job_state']
            (login,host) = jh['Job_Owner'].split('@')
            if login == username:
#                print "XXXX"
                n_still_running+=1
#            owner = getText(job.getElemensByTagName("Job_Owner"))
#            print owner
    #        print job,getText(job.childNodes)                                                                                                                           
#        for jobid in jobs:
#            n_still_running+=1
        return(n_still_running)


def n_still_running(job_ids):
        nn_still_running=0
        still_running=[]
        wc=subprocess.Popen(["qstat","-x"],stdout=subprocess.PIPE)
        out = wc.communicate()
        #print out                                                                                                                                                       
        dom= parseString(out[0])
        #print dom                                                                                                                                                       
    #    print dom.toprettyxml()                                                                                                                                         
        jobs=dom.getElementsByTagName("Job_Id")
    #    for job in jobs:                                                                                                                                                
    #        print job,getText(job.childNodes)                                                                                                                           
        for jobid in job_ids:
            for job in jobs:
                if re.match(jobid,getText(job.childNodes)):
                    nn_still_running+=1
                    still_running.append(jobid)
    #        if wc != "":                                                                                                                                                
    #            n_still_running+=1                                                                                                                                      
    
        return(nn_still_running, still_running)


def wait_for_jobs_to_finish(job_ids,interval=10):
    nn_still_running=1
    still_running=[]
    while nn_still_running>0:
        print "still running:", nn_still_running,still_running

        (nn_still_running, still_running) = n_still_running(job_ids)

        if nn_still_running==0:
            break

        time.sleep(interval)


def submit_job(cmd,qsub_options={"-d":".", "-N":"JAM","-q":"serial" } ): # wd=".",label="JAM"):
    openargs = ["qsub","-V"]
    for f in qsub_options.keys():
        openargs.append(f)
        if qsub_options[f]:
            openargs.append(qsub_options[f])
    print openargs
    wc=subprocess.Popen(openargs,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
    out= wc.communicate(cmd)
    print "out",out
    job_name = out[0]
    tokens = job_name.split(".")
    print "tokens[0]",tokens[0]
    return(tokens[0])
#    job_ids.append(tokens[0])
    

def trickle_jobs(cmdlist,maxjobs,nothing=False,pbs_options={"-d":".", "-N":"JAM","-q":"serial" },interval=10):
    
    job_ids = []    

#    time.sleep(interval)
    
    n_running = n_jobs_for_user(username)
    print "total jobs for user",username,n_running

    if nothing:
        return
            

    while len(cmdlist)>0:
        if not nothing:
            n_running = n_jobs_for_user(username)
#        ( n_running, running_jobs) = n_still_running(job_ids)
        print "n_running", n_running
        for i in range(maxjobs-n_running):
            if len(cmdlist)==0:
                break
            cmd = cmdlist.pop()
            if not nothing:
                print cmd
                job_id = submit_job( cmd , pbs_options)
                time.sleep(1)
                print job_id
                job_ids.append(job_id )

#                job_ids.append( submit_job( cmd , pbs_options ) )

        time.sleep(interval)

    if not nothing:
        wait_for_jobs_to_finish(job_ids,interval)


def trickle_jobs_options(cmdlist,maxjobs,nothing=False,interval=10):

    default_options = {"-d":".", "-N":"JAM","-q":"serial" }
    job_ids = []
    
    print "total jobs for user",username,n_jobs_for_user(username)

#    time.sleep(interval)

    if nothing:
        print "dry run:"
        while len(cmdlist)>0:
            (cmd,options) = cmdlist.pop()
            print "## ",cmd, str(options)
    
    while len(cmdlist)>0:
        n_running = n_jobs_for_user(username)
#        ( n_running, running_jobs) = n_still_running(job_ids)
        print "n_running", n_running
        for i in range(maxjobs-n_running):
            if len(cmdlist)==0:
                break
            (cmd,options) = cmdlist.pop()
            if not nothing:
                print cmd
                job_id = submit_job( cmd , options)
                time.sleep(1)
                print job_id
                job_ids.append(job_id )
            

#                job_ids.append( submit_job( cmd , pbs_options ) )

        time.sleep(interval)

    if not nothing:
        wait_for_jobs_to_finish(job_ids,interval)


