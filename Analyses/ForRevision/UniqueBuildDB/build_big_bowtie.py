import os
import subprocess
import time
from LabQueue.qp import qp, fakeqp
from LabUtils.addloglevels import sethandlers
import configparser
def build_index( ind_path,input_index,bowtie_exe,offrate,threads,ind_name):

    buildCommand='%s --offrate %d --threads %d %s %s/%s' % ( bowtie_exe, offrate, \
                                                             threads, input_index, \
                                                             ind_path, ind_name)
    print ("building")
    print ("%s" % buildCommand )
    status, output = subprocess.getstatusoutput(buildCommand)
    if status != 0:
        print(output)
        print("status:", status)
        raise Exception("failed")
    return

def run(configFile):
    print ("Initializing")
    config = configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation())
    config.read(configFile)
    build_representatives = config['build_bowtie']
    basedir = build_representatives['qp_base_dir']
    if not os.path.exists(basedir):
        os.makedirs(basedir)
    #sethandlers()
    os.chdir(basedir)
    print ("Starting")

    print (time.ctime())
    with fakeqp(jobname='build', q=['himem7.q']) as q:
        q.startpermanentrun()
        ind_path = build_representatives['ind_path']
        input_index = build_representatives['input_index']
        bowtie_exe = build_representatives['bowtie_exe']
        offrate = eval(build_representatives['offrate'])
        threads = eval(build_representatives['threads'])
        ind_name = build_representatives['ind_name']
        if os.path.isdir(ind_path):
            print ("Indexing directory exists. Writing to it")
        else:
            os.mkdir(ind_path)
        waiton = [ q.method( build_index, ( ind_path,input_index,bowtie_exe,offrate,threads,ind_name) ) ]
        print ("job sent")
        q.wait(waiton)
    print(time.ctime())

if __name__=="__main__":
    run()