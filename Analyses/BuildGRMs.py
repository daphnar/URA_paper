###################################################################################################
# File: BuildGRMs.py
# Version: 0.0
# Date: 16.07.2018
# Noam Bar, noam.bar@weizmann.ac.il
#
# 
# Python version: 2.7
###################################################################################################

import pandas as pd
import numpy as np
# import scipy.linalg as la
import os.path
# import subprocess
# import sys
# from scipy.spatial.distance import pdist, squareform
from PNPChip.ForPaper import kernel_utils
# from PNPChip.ForPaper.mirkat import rl_skat
# from PNPChip.albi.albi_lib import fiesta_lib
import Utils
from Analyses.AnalysisHelperFunctions import change_index, make_dir_if_not_exists, remove_rare_elements

from queue.qp import qp, fakeqp
from addloglevels import sethandlers
import argparse
from datetime import datetime

#default paths
GCTA_EXE = '/net/mraid08/export/genie/Bin/gcta/gcta64'
# LMM_DIR = '/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/LMMs/'
GRM_DIR = '/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/LMMs/GRMs/' # TODO: change to your dir
# TODO: add "Segata" to this list
LEGAL_GRMs = ['IGC', 'MPA_species', 'KO', 'KEGGBacGenes', 'MPA_genera', 'IGC-COG', 'IGC-COG-05', 'MPA_species-class', 
              'MPA_species-phyla', 'MPA_species-order', 'MPA_species-kingdom', 'Diet', 'Drugs', 'DietBin', 'Genetics',
              'Age+Sex', 'Anthropometrics', 'Macronutrients', 'Phenome', 'LifeStyle', 'BLOOD',
              'CBC', 'LIVER_FUNCTIONS', 'ELECTROLYTES', 'RENAL_FUNCTIONS', 'LIPIDS', 'ENDOCRINE',
              'TCR0.1'] #TODO: add KEGGBacGenes with some other levels, think of other ways to divide IGC

IGC_PATH = '/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/microbiome_files/EMGenes.csv'
NOG_PATH = '/net/mraid08/export/jafar/Microbiome/Analyses/AllSeqProjects/DFOut/NogGeneSpidDF.dat'
NOG_PATH_DIC_BY_CORR = {'05':'/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/remove_corr_IGC/COG_IGC_rare0.1.dat'}

KEGGBacGenes_PATH = '/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/KEGG_DB_12_12_2017/KEGG_genes_Bacteria_protein_old_DF/BacGenes.csv'
KO_DF_PATH = '/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/KEGG_DB_12_12_2017/KEGG_genes_Bacteria_protein_old_DF/KO_mat_by_KEGG_10_normed.df' # TODO: maybe use other KO matrices
MPA_DF_PATH = '/net/mraid08/export/jafar/Microbiome/Analyses/AllSeqProjects/DFOut/MPASpid.dat'
MPA_PHYLO_PATH = '/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/MetaPhlan2/mpa_phylo_18072018.dat'
MPA_COMPONENT_ABUNDANCE_TH = 0.05
COG_DICT_PATH = '/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/COG2014/cog_dict.dat'
COG_LETTER_DICT_PATH = '/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/COG2014/letter_cog_dict.dat'
NOGFUNC_RANDOM_UNIQ_PATH = '/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/COG2014/NOGfunccat_uniq.dat'
SAMPLES_DICT = {'MAR17':'/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/microbiome_files/FD_SPID/MAR17.dat',
                'MAY18':'/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/microbiome_files/FD_SPID/MAY18.dat',
                'ACS':'/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/microbiome_files/FD_SPID/ACS.dat'}
RARE_ELEMENTS_TH = 0.05

DIET_DF_PATH = '/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/microbiome_files/metabolon/pnp_diet.dat'
DIETBIN_DF_PATH = '/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/microbiome_files/metabolon/pnp_diet_binary.dat'
DRUGS_DF_PATH = '/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/microbiome_files/metabolon/pnp_drugs.dat'
BLOOD_DF_PATH = '/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/microbiome_files/metabolon/pnp_blood_features.dat'
PNP_FEATURES_PATH = '/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/microbiome_files/metabolon/pnp_comb_features.dat'
TCR_01_DF_PATH = '/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/Metabolon/tmp_files/Shani_TCR_0.1.dat'

macronut = ['Daily_Carbohydrate_g', 'Daily_Energy_kcal', 'Daily_Protein_g', 'Daily_TotalLipid_g', 'Daily_Water_g']
anthropometrics = ['BMI', 'Waist', 'Hips', 'Height', 'Weight', 'WHR']
age_gender = ['Age', 'Gender']
phenome = ['HeartRate', 'systolic', 'diastolic', 'glycemic_index']

lifestyle = ['Daily_exercisetime', 'Daily_sleeptime', 'Daily_midday_sleep', 'Currently_smokes', 'Ever_smoked']
ml_cat = {'Age+Sex':age_gender, 'Anthropometrics':anthropometrics, 'Macronutrients':macronut, \
          'Phenome':phenome, 'LifeStyle':lifestyle}

BLOOD_DIC = {'CBC':['WBC', 'Hemoglobin', 'HCT', 'MCH', 'MCHC', 'MCV', 'MPV', 'RBC', 'RDW', 'Basophils %', 'Eosinophils %', 'Lymphocytes %', 'Monocytes %', 'Neutrophils %', 'Platelets'],
            'LIVER_FUNCTIONS':['ALT', 'AST', 'Albumin'],
            'ELECTROLYTES':['Potassium', 'Sodium', 'Calcium', 'Chloride', 'Phosphorus'],
            'RENAL_FUNCTIONS':['Creatinine'],
            'LIPIDS':['Cholesterol, total', 'HDL Cholesterol'],
            'ENDOCRINE':['HbA1C%', 'TSH']}

def get_FD_SPID_list(command_args):
    print ('get_FD_SPID_list')
    samples_list = []
    for samps in command_args.samples:
        samples_list += Utils.Load(SAMPLES_DICT[samps])
    return list(set(samples_list))


def get_microbiome_df(microbiome_level, command_args, grm, mc=None):

    ##########EDIT
    print ('get_microbiome_df')
    if not command_args.use_quantile_normalization:
        if microbiome_level.startswith('MPA'):
            mpa = Utils.Load(MPA_DF_PATH).fillna(0)
            tax = microbiome_level.split('_')[1][0] # MPA_species -> 's'
            return mpa.loc[tax]
        elif microbiome_level == 'IGC':
            if mc == 'COG':
                if len(grm.split('-')) == 2:
                    return Utils.Load(NOG_PATH).fillna(0).T
                return Utils.Load(NOG_PATH_DIC_BY_CORR[grm.split('-')[2]]).fillna(0).T
            return pd.read_csv(IGC_PATH, index_col=0)#.fillna(0)
        elif microbiome_level == 'KEGGBacGenes':
            bac = pd.read_csv(KEGGBacGenes_PATH, index_col=0).fillna(0)
            return bac.apply(lambda x: x/x.sum())
        elif microbiome_level == 'KO':
            return Utils.Load(KO_DF_PATH).fillna(0).T
        else:
            print ("illegal choice of microbiome")
            exit
    else:
        # TODO read matrices after quantile normalization
        pass
        
def build_grm(command_args, grm):
    command_args.output_dir += '/' + grm + '/'
    make_dir_if_not_exists(command_args.output_dir)
    print ('build_grm')
    samples = get_FD_SPID_list(command_args)
    # load relevant data
    db = grm.split('-')
    mc = None
    if len(db) > 2: mc = db[1]
    db = db[0]
    if db in ['IGC', 'MPA_species', 'MPA_genera', 'KO', 'KEGGBacGenes']:
        data = get_microbiome_df(db, command_args, grm, mc).loc[:, samples].dropna(how='all', axis=1)
        data = remove_rare_elements(data, null=False, rare_def=RARE_ELEMENTS_TH)
    elif db == 'Diet':
        data = Utils.Load(DIET_DF_PATH).loc[samples].fillna(0)
    elif db == 'DietBin':
        data = Utils.Load(DIETBIN_DF_PATH).loc[samples].fillna(0)    
    elif db == 'Drugs':
        data = Utils.Load(DRUGS_DF_PATH).loc[samples].fillna(0)
    elif db in ml_cat:
        data = Utils.Load(PNP_FEATURES_PATH).replace({'Gender':{'Female':0, 'Male':1}}).loc[samples, ml_cat[db]] # fillna(0) ?
        data = data.fillna(data.median())
    elif db == 'BLOOD':
        data = Utils.Load(BLOOD_DF_PATH).loc[samples].drop(['CRP (WIDE RANGE)', 'CRP hs'], axis=1) # fillna(0) ?
        data = data.fillna(data.median())
    elif db in BLOOD_DIC:
        data = Utils.Load(BLOOD_DF_PATH).loc[samples, BLOOD_DIC[db]] # fillna(0) ?
        data = data.fillna(data.median())
    elif db == 'TCR0.1':
        data = Utils.Load(TCR_01_DF_PATH).loc[samples].dropna()
    else:
        exit
    # check for multi-component
    if mc is not None:
        divide_into_multicomponent(data, db, mc, command_args)
    else:
        compute_grm_and_save(data, command_args, db)
    return

def compute_grm_and_save(df, command_args, name):
    print ('compute_grm_and_save')
#     norm = False
    norm = True
    if name.startswith('IGC') or name.startswith('KEGG') or name.startswith('KO'): # should MPA be here as well?
        norm = True
        df = (df.T > command_args.presence_absence_th).astype(np.float) # TODO: pick a relevant threshold for presence absence
    K_mb = _compute_as_snp(df, norm)  
    df_K = pd.DataFrame(K_mb, index=df.index, columns=df.index)
    Utils.Write(command_args.output_dir + name + '.dat', df_K)
    kernel_utils.write_grm(df_K.values, command_args.output_dir + name + '.grm.gz')
    df_K['IID'] = df_K.index
    df_K['IID'].to_csv(command_args.output_dir + name + '.grm.id', sep='\t', index_label='IID')
    del df_K['IID']
    return
    
    
def divide_into_multicomponent(grm, db, mc, command_args):
    print ('divide_into_multicomponent')
    # for some predetermined rules, divide the dataframe into component
    if db == 'IGC':
        if mc == 'COG':
            mc_IGC_COG(grm, command_args)
        else:
            pass
    elif db.startswith('MPA'):
        mc_MPA(grm, db, mc, command_args)
    elif db.startswith('Diet'):
        mc_Diet(grm, db, mc, command_args)
    else:
        pass
    
    
def mc_IGC_COG(grm, command_args):
    print ('mc_IGC_COG')
    cog_dict = Utils.Load(COG_DICT_PATH)
    letter_cog_dict = Utils.Load(COG_LETTER_DICT_PATH)
    NOGfunccat_uniq = Utils.Load(NOGFUNC_RANDOM_UNIQ_PATH)
    nog_letters = grm.T.rename(columns=NOGfunccat_uniq).copy()
    # letters to main pathways
    nog_pathways = nog_letters.rename(columns=letter_cog_dict).copy()
    # main pathways
    for cog in cog_dict:
        temp_cog = nog_pathways[cog].T
#         Utils.Write(command_args.output_dir + cog + '.dat', temp_cog)
        compute_grm_and_save(temp_cog, command_args, cog)
    # all but main pathways
    nog_pathways_rest = nog_pathways.loc[:, ~nog_pathways.columns.isin(cog_dict.keys())].T
#     Utils.Write(command_args.output_dir + '/UNKNOWN.dat', nog_pathways_rest)
    compute_grm_and_save(nog_pathways_rest, command_args, 'UNKNOWN')
    return

def mc_MPA(grm, db, mc, command_args):
    print ('mc_MPA')
    mpa_phylo = Utils.Load(MPA_PHYLO_PATH)
    working_tax_level = db.split('_')[1][0]
    component_tax_level = mc[0]

    temp_mpa_phylo = mpa_phylo[mpa_phylo[working_tax_level].isin(grm.index)].drop('t', axis=1).drop_duplicates()
    abund_tax = temp_mpa_phylo[component_tax_level].value_counts()[temp_mpa_phylo[component_tax_level].value_counts() > temp_mpa_phylo.shape[0] * MPA_COMPONENT_ABUNDANCE_TH].index
    
    for tax in abund_tax:
        temp_subtax = mpa_phylo[mpa_phylo[component_tax_level] == tax][working_tax_level].dropna()
        rel_tax = set(temp_subtax).intersection(set(grm.index))
        temp_mpa = grm.loc[rel_tax].dropna()
#         Utils.Write(command_args.output_dir + '/' + tax + '.dat', temp_mpa)
        compute_grm_and_save(temp_mpa, command_args, tax)
    temp_subtax = mpa_phylo[~mpa_phylo[component_tax_level].isin(abund_tax)][working_tax_level].dropna()
    rel_tax = set(temp_subtax).intersection(set(grm.index))
    temp_mpa = grm.loc[rel_tax].dropna()
#     Utils.Write(command_args.output_dir + '/rest.dat', temp_mpa)
    compute_grm_and_save(temp_mpa, command_args, 'rest')
    return

def mc_Diet(grm, db, mc, command_args):
    pass
    

def _compute_as_snp(df, norm=True):
    # in case of presence absence remove genes which are all 1
    df = df.loc[:, ((df==1).sum() != df.shape[0])].copy()
    # in case of presence absence remove genes which are all 0
    df = df.loc[:, ((df!=0).sum() > df.shape[0]*0.01)]
    # normalizing columns
    if norm:
        df_norm = df.apply(lambda x: (x - x.mean()) / x.std())
    else:
        df_norm = df.apply(lambda x: (x - x.mean()))
    # multiplying matrix by its transpose and dividing by n (X*X^T / n)
    df_norm = df_norm.dot(df_norm.T) / df_norm.shape[1]
    return df_norm.values
    
    
def upload_job_buildGRM(q, command_args):
    print ('upload_job_buildGRM')
    
    waiton = []
    for grm in command_args.grm:
        waiton.append(q.method(build_grm, (command_args, grm))) 
    res = q.waitforresults(waiton)
    return res

def _convert_comma_separated_to_list(s):
    return s.split(',')

def main():
    print ('main')
    parser = argparse.ArgumentParser()
    parser.add_argument('grm', help='What are the random effects. (IGC, IGC-COG, MPA_species, MPA_species-phyla, KEGGBacGenes, etc.) separated by comma', type=str)
    parser.add_argument('samples', help='What samples to use. (ACS, MAR17, MAY18, etc.) separated by comma', type=str)
    parser.add_argument('-output_dir', help='Path to output directory', type=str, default=GRM_DIR)
    parser.add_argument('-use_quantile_normalization', help='Whether to use quantile normalization over microbiome data', type=bool, default=False)
    parser.add_argument('-presence_absence_th', help='What should be the presence absence threshold for microbiome abundances', type=float, default=1e-6)
#     parser.add_argument('-IGC_corr', help='Correlation threshold to use for IGC matrix', type=float, default=None)
    
    command_args = parser.parse_args()
    
    make_dir_if_not_exists(command_args.output_dir)
    
    command_args.grm = _convert_comma_separated_to_list(command_args.grm)
    for grm in command_args.grm:
        if grm not in LEGAL_GRMs:
            print ('grm currently supported are: ' + ','.join(LEGAL_GRMs))
            exit
            
    command_args.samples = _convert_comma_separated_to_list(command_args.samples)
    for samps in command_args.samples:
        if samps not in SAMPLES_DICT:
            print ('samples currently supported are: ' + ','.join(SAMPLES_DICT.keys()))
            exit
            
    if command_args.use_quantile_normalization:
        command_args.output_dir += '/QN/'
        make_dir_if_not_exists(command_args.output_dir)
        
    command_args.output_dir += '/' + '-'.join(command_args.samples) + '/'
    
    make_dir_if_not_exists(command_args.output_dir)
    
    with open(command_args.output_dir + '/args' + str(datetime.now()), 'w') as handle:
        for arg in vars(command_args):
            handle.write(str(arg) + '\t' + str(getattr(command_args, arg)) + '\n')
    
    
#     with qp(jobname = 'GRM-build', q=['himem7.q'], mem_def = '30G', trds_def = 1, tryrerun=False,
    with fakeqp(jobname = 'GRM-build', q=['himem7.q'], mem_def = '30G', trds_def = 1, tryrerun=False,
        max_u = 25, delay_batch=15) as q:
        os.chdir("/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/temp_q_dir/") # TODO: change to your dir
        q.startpermanentrun()
        upload_job_buildGRM(q, command_args)

    return


if __name__ == '__main__':
    sethandlers()
    main()





df = #TODO: load your microbiome data

df = (df.T > presence_absence_th).astype(np.float) # TODO: pick a relevant threshold for presence absence
# in case of presence absence remove genes which are all 1
df = df.loc[:, ((df == 1).sum() != df.shape[0])].copy()
# in case of presence absence remove genes which are all 0
df = df.loc[:, ((df != 0).sum() > df.shape[0] * 0.01)]

df_norm = df.apply(lambda x: (x - x.mean()) / x.std())

K_mb = df_norm.dot(df_norm.T) / df_norm.shape[1].values

df_K = pd.DataFrame(K_mb, index=df.index, columns=df.index)
Utils.Write(output_dir + name + '.dat', df_K) # TODO: save the grm matrix
kernel_utils.write_grm(df_K.values, output_dir + name + '.grm.gz')
df_K['IID'] = df_K.index
df_K['IID'].to_csv(output_dir + name + '.grm.id', sep='\t', index_label='IID')