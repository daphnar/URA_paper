import os
import pandas as pd
from Bio import SeqIO
import time
import glob
import configparser
from LabQueue.qp import qp, fakeqp
from LabUtils.addloglevels import sethandlers


# def renameCore(SelectedSGBs,output_cores_dir,genomes_dir):
#     all_contig_lens=[]
#     all_contig_ids=[]
#     for sgb in SelectedSGBs:
#         matchingFasta=glob.glob(os.path.join(genomes_dir,'SGB_%s_*.fa'%sgb))
#         coreContigs=[]
#         if len(matchingFasta)!=1:
#             print("no core found for %s under %s"%(sgb,genomes_dir))
#             continue
#         matchingFasta=matchingFasta[0]
#         coreFasta = os.path.join(output_cores_dir, os.path.basename(matchingFasta))
#         if not os.path.exists(coreFasta):
#             count = 0
#             for record in SeqIO.parse(matchingFasta, "fasta"):
#                 genome = os.path.basename(matchingFasta).split('SGB_%s_'%sgb)[1][:-3]
#                 record.id = 'SGB_%s_c_%s_%s' % (sgb, count, genome)
#                 #record.id = 'SGB_%s_core_c_%s' % (sgb, count)
#                 all_contig_ids.append(record.id)
#                 count+=1
#                 coreContigs.append(record)
#                 all_contig_lens.append(len(record))
#             SeqIO.write(coreContigs,coreFasta ,'fasta')
#     pd.Series(index=all_contig_ids,data=all_contig_lens).to_csv(\
#         os.path.join(output_cores_dir, 'all_contigs.txt'),sep='\t',header=False)
#     print ("Done writing core chunk fasta")

def buildByCore(SelectedSGBs,output_fasta,genomes_dir):
    if os.path.exists( output_fasta ):
        os.system( "rm -f %s" % output_fasta )
    first = True
    assert len(SelectedSGBs)==len(set(SelectedSGBs))
    for sgb in SelectedSGBs:
        matchingFasta=glob.glob(os.path.join(genomes_dir,'SGB_%s_*.fa'%sgb))
        if len(matchingFasta)!=1:
            print("no core found for %s under %s"%(sgb,genomes_dir))
            continue
        matchingFasta=matchingFasta[0]
        if first:
            os.system('cat %s > %s'%(matchingFasta,output_fasta))
            first = False
        else:
            os.system('cat %s >> %s' % (matchingFasta, output_fasta))
    print ("Done writing one big fasta for bowtie index")

def run(SelectedSGBs,configFile):
    config = configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation())
    config.read(configFile)
    build_representatives = config['build_representatives']
    print ("Initializing")
    basedir = build_representatives['qp_base_dir']
    if not os.path.exists(basedir):
        os.makedirs(basedir)
    if not os.path.exists(build_representatives['output_cores_dir']):
        os.makedirs(build_representatives['output_cores_dir'])
    sethandlers()
    os.chdir(basedir)
    print ("Starting")
    print (time.ctime())
    with fakeqp(jobname='build', q=['himem7.q']) as q:
        q.startpermanentrun()
        waiton=[]
        chunk_size = eval(build_representatives['chunksize'])
        # for chunk in range(0,len(SelectedSGBs), chunk_size):
        #     waiton.append(q.method(renameCore, (SelectedSGBs[chunk:chunk+chunk_size],
        #                                  build_representatives['output_cores_dir'],
        #                                  build_representatives['genomes_dir'])))
        # q.wait(waiton)
        waiton = [q.method(buildByCore, (SelectedSGBs,
                                         build_representatives['output_fasta'],
                                         build_representatives['output_cores_dir']))]
        q.wait(waiton)
    print (time.ctime())

if __name__=='__main__':
    run()