import numpy as np
# # from Unicorn.build_pan_genome.buildUniqueDB.createDistantConfigFile import createDistantConfigFile
# from URA_paper.Analyses.ForRevision.UniqueBuildDB import build_big_bowtie, \
#     buildRepresentatives, make_fict_reads, map_fict_reads,analyse_fict_maps,unite_fict_maps
import time
import configparser
import os
import pandas as pd
import glob
from LabQueue.qp import qp, fakeqp
from LabUtils.addloglevels import sethandlers

from URA_paper.Analyses.ForRevision.UniqueBuildDB import make_fict_reads, map_fict_reads, analyse_fict_maps, \
    unite_fict_maps
from URA_paper.Analyses.ForRevision.UniqueBuildDB.createDistantConfigFile import createDistantConfigFile

def getAllSGBs(representatives,genomes_dir,all_large_or_new_sgbs):
    if not os.path.exists(representatives):
        sgbs={}
        keepSGBs=pd.Series.from_csv(all_large_or_new_sgbs).values
        for fastafile in glob.glob(os.path.join(genomes_dir,'*.fa')):
            sgb=os.path.basename(fastafile).split('_')[1]
            if int(sgb) in keepSGBs:
                rep=os.path.basename(fastafile).split(sgb+'_')[1][:-3]
                sgbs[rep]=sgb
        pd.Series(sgbs).to_csv(representatives,sep='\t')
    representatives_df = pd.read_table(representatives, header=None)
    representatives_df.columns = ['nameOfFile', 'SGB']
    return representatives_df['SGB']

def sgb_run(sgb,configFile):
    make_fict_reads.run_sample(sgb,configFile)
    map_fict_reads.run_sample(sgb,configFile)
    analyse_fict_maps.run_sample(sgb,configFile)
    make_fict_reads.clean(sgb,configFile)
    map_fict_reads.clean(sgb,configFile)
    unite_fict_maps.run_sample(sgb,configFile)

def runChuckOfSGBs(SelectedSGBs,configFile):
    config = configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation())
    config.read(configFile)
    run_pipeline = config['run_pipeline']
    score_output = run_pipeline['score_output']
    for sgb in SelectedSGBs:
        if len(glob.glob(os.path.join(score_output, 'dict_scores_SGB_%s_*' % sgb))) == 0 \
                or len(glob.glob(os.path.join(score_output, 'dict_bins_SGB_%s_*' % sgb))) == 0:
            sgb_run(sgb, configFile)
        else:
            print ("Done on %s" % sgb)

def runOnSGBs(configFile):
    config = configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation())
    config.read(configFile)
    run_pipeline = config['run_pipeline']
    SelectedSGBs=getAllSGBs(run_pipeline['representatives'],
               run_pipeline['genomes_dir'],
               run_pipeline['all_large_or_new_sgbs'])
    if not os.path.exists(run_pipeline['stage1output']):
        print ("Making representatives fasta", time.ctime())
        buildRepresentatives.run(SelectedSGBs,configFile)
        print ("Bulding Bowtie index", time.ctime())
        build_big_bowtie.run(configFile)
        with open(run_pipeline['stage1output'],'w') as donefile:
            donefile.write('Done\n')
    basedir = run_pipeline['qp_base_dir']
    score_output = run_pipeline['score_output']
    sethandlers()
    os.chdir(basedir)
    print ("Starting")
    with fakeqp(jobname='build', q=['himem7.q']) as q:
        q.startpermanentrun()
        waiton = []
        chucksize=50
        count=0
        for chunkSGBsIDs in range(0,len(SelectedSGBs),chucksize):
            chunkSGBs=SelectedSGBs.loc[count*chucksize:chucksize*(count+1)-1]
            count+=1
        q.wait(waiton)
        print ("Done running on %s SGBs"%len(waiton))
    print ("Done", time.ctime())
    return

if __name__=="__main__":
    configFile=createDistantConfigFile(forceRewrite=False)
    runOnSGBs(configFile)
