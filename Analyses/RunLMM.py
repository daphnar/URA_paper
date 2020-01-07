###################################################################################################
# File: RunLMM.py
# Version: 0.0
# Date: 16.07.2018
# Noam Bar, noam.bar@weizmann.ac.il
#
#
# Python version: 2.7
###################################################################################################

import pandas as pd
import numpy as np
import scipy.linalg as la
import os.path
import subprocess
import sys
from scipy.spatial.distance import pdist, squareform
from PNPChip.ForPaper import kernel_utils
from PNPChip.ForPaper.mirkat import rl_skat
from PNPChip.albi.albi_lib import fiesta_lib
import Utils
from Analyses.AnalysisHelperFunctions import change_index, make_dir_if_not_exists, remove_rare_elements
from Analyses.BuildGRMs import get_FD_SPID_list, SAMPLES_DICT, GRM_DIR, LEGAL_GRMs
import scipy.stats as stats
import re

from queue.qp import qp, fakeqp
from addloglevels import sethandlers
import argparse
from datetime import datetime

#default paths
GCTA_EXE = '/net/mraid08/export/genie/Bin/gcta/gcta64'
LMM_DIR = '/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/LMMs/Results/'
LEGAL_PHENOTYPES = ['BMI', 'Height', 'Age', 'Metabolomics_raw', 'Metabolomics_normed', 'Metabolomics_grouped']
LEGAL_FIXED_EFFECTS = ['Age', 'Gender', 'Nextera', 'Drugs']
# PNP_FEATURES_PATH = '/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/microbiome_files/pnp1_features_data_drugs_20072018.dat'
PNP_FEATURES_PATH = '/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/microbiome_files/metabolon/pnp_comb_features.dat'
ACS_FEATURES_PATH = '/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/microbiome_files/cardio_features_19062018.dat'
DRUGS_LIST = ['D.Contraception', 'D.CVD', 'D.lipid', 'D.NSAID', 'D.pain', 'D.Thyroid', 'D.Bone/Joint', 'D.GI', 'D.Psychiatric', 'D.All']
# PNP_DRUGS = Utils.Write('/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/microbiome_files/pnp_drugs.dat')

STOOL_INFO_DIC = {'MAR17':'/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/Cardio/Metabolomics/dataframes/mar17_stool_info.dat',
                  'MAY18':'/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/Cardio/Metabolomics/dataframes/may18_stool_info.dat',
                  'ACS':'/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/Cardio/Metabolomics/dataframes/cardio_stool_info.dat'}


UNJOINED_NORMED_DIC = {'MAR17':'/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/Cardio/Metabolomics/dataframes/inrun_normed1.dat',
                       'MAY18':'/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/Cardio/Metabolomics/dataframes/inrun_normed2.dat',
                       'ACS':'/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/Cardio/Metabolomics/dataframes/inrun_normed2.dat'}
UNJOINED_UNNORMED_DIC = {'MAR17':'/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/Cardio/Metabolomics/dataframes/os_1_cid.dat',
                       'MAY18':'/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/Cardio/Metabolomics/dataframes/os_2_cid.dat',
                       'ACS':'/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/Cardio/Metabolomics/dataframes/os_2_cid.dat'}
JOINED_NORMED_DIC = {'MAR17':'/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/Cardio/Metabolomics/dataframes/joined1.dat',
                     'MAY18':'/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/Cardio/Metabolomics/dataframes/joined2.dat',
                     'ACS':'/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/Cardio/Metabolomics/dataframes/joined2.dat'}
JOINED_UNNORMED_DIC = {'MAR17':'/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/Cardio/Metabolomics/dataframes/joined_not_norm_within_run1.dat',
                       'MAY18':'/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/Cardio/Metabolomics/dataframes/joined_not_norm_within_run2.dat',
                       'ACS':'/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/Cardio/Metabolomics/dataframes/joined_not_norm_within_run2.dat'}
GROUPED_DIC = {'MAR17':'/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/Metabolon/SHAP/dataframes/mar17_metabolomics_grouped085_unnormed_fillna_min_dayfromfirstsample_regressed_rzs_sample_id.dat'}

UNJOINED_METABOLOMICS_DIC = {'raw':UNJOINED_UNNORMED_DIC, 'normed':UNJOINED_NORMED_DIC, 'grouped':GROUPED_DIC}
JOINED_METABOLOMICS_DIC = {'raw':JOINED_UNNORMED_DIC, 'normed':JOINED_NORMED_DIC}
# UNJOINED_METABOLOMICS_DIC = {'MAR17':'/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/Cardio/Metabolomics/dataframes/os_1_cid.dat',
#                              'MAY18':'/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/Cardio/Metabolomics/dataframes/os_2_cid.dat',
#                              'ACS':'/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/Cardio/Metabolomics/dataframes/os_2_cid.dat'}
# JOINED_METABOLOMICS_DIC = {'MAR17':'/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/Cardio/Metabolomics/dataframes/joined1.dat',
#                            'MAY18':'/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/Cardio/Metabolomics/dataframes/joined2.dat',
#                            'ACS':'/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/Cardio/Metabolomics/dataframes/joined2.dat'}
METABOLOMICS_DIC = {True:JOINED_METABOLOMICS_DIC, False:UNJOINED_METABOLOMICS_DIC}
FEATURES_TO_NUM_DIC = {'Gender':{'Male':0, 'Female':1}, 'DM':{'Yes':1, 'No':0}, 'Smoking':{'Yes':1, 'No':0, 'Past':0}}
# .replace({'Gender':{'Female':0, 'Male':1}}).replace({True:1, False:0})


def get_phenotype_df(command_args):
    if not command_args.phenotype.startswith('Metabolomics'):
        return get_features(command_args, [command_args.phenotype])
    elif command_args.phenotype.startswith('Metabolomics'):
        joined = False
        if 'ACS' in command_args.samples and ('MAR17' in command_args.samples or 'MAY18' in command_args.samples):
            joined = True
        metabolomics = pd.DataFrame()
        for sample in command_args.samples:
            metabolomics = pd.concat((metabolomics, load_metabolomics(sample, joined, command_args.phenotype.split('_')[1])))
        return metabolomics
    else:
        print ('phenotype currently supported are: ' + ','.join(LEGAL_PHENOTYPES))
        exit

def load_metabolomics(sample, joined, phenotype):
    metabolomics = Utils.Load(METABOLOMICS_DIC[joined][phenotype][sample])
    stool_info = Utils.Load(STOOL_INFO_DIC[sample])
    metabolomics = metabolomics.loc[stool_info.SAMPLE_ID]
    if not joined and phenotype != 'grouped': metabolomics = metabolomics.apply(np.log10).dropna(how='all')
    return change_index(metabolomics, stool_info, from_idx='SAMPLE_ID', to_idx='FD_SPID')

def get_features(command_args, features2take):
    if 'Drugs' in features2take:
        if 'ACS' in command_args.samples:
            print ('Drugs currently not supported for ACS samples')
            exit
        features2take.remove('Drugs')
        features2take += DRUGS_LIST
        pnp_features = Utils.Load(PNP_FEATURES_PATH).replace(FEATURES_TO_NUM_DIC)[features2take]
        return pnp_features.replace({True:1, False:0})
#     pnp_features = Utils.Load(PNP_FEATURES_PATH).replace(FEATURES_TO_NUM_DIC)[features2take]
#     acs_features = Utils.Load(ACS_FEATURES_PATH).replace(FEATURES_TO_NUM_DIC)[features2take]
    pnp_features, acs_features = pd.DataFrame(), pd.DataFrame()
    if 'MAR17' in command_args.samples or 'MAY18' in command_args.samples:
        pnp_features = Utils.Load(PNP_FEATURES_PATH).replace(FEATURES_TO_NUM_DIC)[features2take]
    if 'ACS' in command_args.samples:
        acs_features = Utils.Load(ACS_FEATURES_PATH).replace(FEATURES_TO_NUM_DIC)[features2take]
    features = pd.concat((pnp_features, acs_features)).replace({True:1, False:0})
    for feature in features2take:
        if feature not in features.columns:
            print (feature + ' is missing from the dataframe, check if it exists for the samples you chose.')
            exit
    return features


def _parse_GCTA_results(stdout, phenotype_name):
    intercept, sig2g, sig2e, h2_reml, h2_reml_std, h2_reml_l, h2_reml_l_std = None, None, None, None, None, None, None
    fixed_effects = []
    for line in stdout.split('\n'):
        line_split = line.split()
        if (len(line_split) == 0): continue
        if (line_split[0] == 'mean'): intercept = float(line_split[1])
        if (line_split[0][:2] == 'X_'): fixed_effects.append(float(line_split[1]))
        if (line_split[0] == 'V(G)'): sig2g = float(line_split[1])
        if (line_split[0] == 'V(e)'): sig2e = float(line_split[1])
        if ('error' in line.lower()):
            print ('GCTA', phenotype_name, line)
            return intercept, sig2g, sig2e, h2_reml, h2_reml_std, h2_reml_l, h2_reml_l_std, fixed_effects
        if line.startswith('V(G)/Vp') and 'V(G)/Vp_L' not in line: h2_reml, h2_reml_std = float(line_split[1]), float(line_split[2])
        if line.startswith('V(G)/Vp_L'): h2_reml_l, h2_reml_l_std = float(line_split[1]), float(line_split[2])
    return intercept, sig2g, sig2e, h2_reml, h2_reml_std, h2_reml_l, h2_reml_l_std, fixed_effects

def _parse_GCTA_results_MC(stdout, phenotype_name):
    estimates = []
    intercept, sum_of_estimates, vp = None, None, None
    fixed_effects = []
    for line in stdout.split('\n'):
        line_split = line.split()
        if (len(line_split) == 0): continue
        if (line_split[0] == 'mean'): intercept = float(line_split[1])
        if (line_split[0][:2] == 'X_'): fixed_effects.append(float(line_split[1]))
        if ('error' in line.lower()):
            print ('GCTA', phenotype_name, line)
            return intercept, fixed_effects, estimates, sum_of_estimates
#         if re.match('V\(G[0-9]\)/Vp.*', line): estimates.append(float(line.split()[1]))
        if re.match('V\(G[0-9]\).*', line) and '/Vp' not in line: estimates.append(float(line.split()[1]))
        if line.startswith('Vp'): vp = float(line.split()[1])
        if line.startswith('Sum of'): sum_of_estimates = float(line.split()[-2])
    estimates = [e/vp for e in estimates]
    return intercept, fixed_effects, estimates, sum_of_estimates

def estimate_h2_single_kernel(command_args, df_pheno, df_covariates):
    pheno_file = os.path.join(command_args.output_dir, 'temp.phe')
    covariates_file = os.path.join(command_args.output_dir, 'temp.cov')
    grm_file = GRM_DIR + '-'.join(command_args.samples) + '/' + command_args.random[0] + '/' + command_args.random[0]

    df_covariates['FD_SPID'] = df_covariates.index.values
    df_covariates = df_covariates[['FD_SPID'] + [c for c in df_covariates.columns if c != 'FD_SPID']]
    df_covariates.to_csv(covariates_file, sep='\t', header=False, na_rep='NA')
    del df_covariates['FD_SPID']

    with open(command_args.output_dir + 'LMM_results.txt', 'w') as handle:
        handle.write('\t'.join(['Phenotype', 'microbiome-association index', '95% CI', 'P value',
                                'Sample size', 'V(G)', 'V(e)', 'mean'] + list(df_covariates.columns)) + '\n')

    cov_mat = Utils.Load(GRM_DIR + '-'.join(command_args.samples) + '/' + command_args.random[0] + '/' + command_args.random[0] + '.dat')
    df_K = pd.DataFrame(cov_mat.values, index=cov_mat.index, columns=cov_mat.index)
    df_K['IID'] = df_K.index
    del df_K['IID']

    for phenotype_name in df_pheno.columns:

        #remove outlier phenotypes
        is_outlier = np.abs(df_pheno[phenotype_name] - df_pheno[phenotype_name].mean()) > 5*df_pheno[phenotype_name].std()
        if (np.sum(is_outlier) > 0): df_pheno.loc[is_outlier, phenotype_name] = np.nan
        # check if there are enough samples
        if df_pheno[phenotype_name].notnull().sum() < 20:
            print ("Not enough samples", str(phenotype_name))
            continue

        df_pheno_temp = df_pheno[[phenotype_name]].copy()
        df_pheno_temp['FD_SPID'] = df_pheno_temp.index.values
        df_pheno_temp[['FD_SPID', phenotype_name]].to_csv(pheno_file, sep='\t', header=False, na_rep='NA', index=True)
        del df_pheno_temp['FD_SPID']

        #invoke GCTA
        cmdLine = [GCTA_EXE, '--grm-gz', grm_file, '--reml', '--pheno', pheno_file,
                   '--qcovar', covariates_file, '--reml-est-fix', '--thread-num', '1']
        proc = subprocess.Popen(cmdLine, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        stdout, stderr = proc.communicate()

        if (stderr is not None):
            print ('GCTA reml error:')
            print stderr
            raise Exception()

        intercept, sig2g, sig2e, h2_reml, h2_reml_std, h2_reml_l, h2_reml_l_std, fixed_effects = _parse_GCTA_results(stdout, phenotype_name)

        if h2_reml is None:
            print (stdout)
            continue
        #compute a p-value with RL-SKAT
        df_merged = df_K.merge(df_covariates.dropna(), left_index=True, right_index=True).copy()
        df_merged = df_merged.merge(df_pheno_temp[[phenotype_name]], left_index=True, right_index=True)
        df_merged = df_merged.loc[df_merged[phenotype_name].notnull()]
        covariates = df_merged[[c for c in df_covariates.columns if c != 'FD_SPID']].copy()
        covariates['intercept'] = 1.0
        covariates = covariates.values
        pheno_vec = df_merged[phenotype_name].values
        df_K_subset = df_K.loc[df_K.index.isin(df_merged.index), df_K.columns.isin(df_merged.index)]
        df_K_subset = df_K_subset.loc[df_merged.index, df_merged.index]
        assert (df_K_subset.index == df_K_subset.columns).all()
        assert (df_K_subset.shape[0] == df_merged.shape[0])
        assert (df_K_subset.index.isin(df_merged.index).all())
        assert (df_K_subset.index == df_merged.index).all()
        kernel = df_K_subset.values
        pvalues = rl_skat.RL_SKAT_Full_Kernel(kernel, covariates, add_intercept=False).test(np.row_stack(pheno_vec))
        #compute CIs with FIESTA
        s,U = la.eigh(kernel)
        ind = np.argsort(s)[::-1]; s=s[ind]; U=U[:,ind]
        CI = fiesta_lib.calculate_cis_general(np.array([h2_reml]), s, U, covariates,
                                              iterations=command_args.fiesta_iterations, alpha=0.05,
                                              tau=0.4, use_convergence_criterion=False,
                                              use_progress_bar=False) # iterations=1000
        line =  '\t'.join([str(c) for c in [
                         phenotype_name,
                         '%0.3g'%(h2_reml), '%0.3g - %0.3g'%(CI[0][0], CI[0][1]),
                         '%0.3g'%(pvalues[0]), pheno_vec.shape[0], '%0.3g'%(sig2g), '%0.3g'%(sig2e),
                         '%0.3f'%(intercept)] + ['%0.3g'%(f) for f in fixed_effects]]) + '\n'
        with open(command_args.output_dir + 'LMM_results.txt', 'a') as handle:
            handle.write(line)
        sys.stdout.flush()


def estimate_h2_multi_kernel(command_args, df_pheno, covariates, idx):
    covariates_file = os.path.join(command_args.output_dir, 'temp.cov')
    with open(command_args.output_dir + '/multi_grm.txt', 'r') as handle:
        grms = handle.readlines()
        grms = [g.strip().split('/')[-1] for g in grms]

    with open(command_args.output_dir + 'results' + str(idx), 'w') as handle:
        columns = ['Phenotype', 'Sample size', 'mean', 'Jackknife'] + \
        [item for sublist in [[f + ' : ' + s for s in ['microbiome-association index', '95% CI', 'P value', 'SE']] for f in grms + ['Sum of components']] for item in sublist] + \
        list(covariates[1:])
        handle.write('\t'.join(columns) + '\n')

    for phenotype_name in df_pheno.columns:
        pheno_file = os.path.join(command_args.output_dir, str(phenotype_name) + '.phe')

        #remove outlier phenotypes
        is_outlier = np.abs(df_pheno[phenotype_name] - df_pheno[phenotype_name].mean()) > 5*df_pheno[phenotype_name].std()
        if (np.sum(is_outlier) > 0): df_pheno.loc[is_outlier, phenotype_name] = np.nan
        # check if there are enough samples
        if df_pheno[phenotype_name].notnull().sum() < 100:
            print ("Not enough samples", str(phenotype_name))
            continue

        df_pheno_temp = df_pheno[[phenotype_name]].copy()
        df_pheno_temp = df_pheno_temp.dropna()
        df_pheno_temp['FD_SPID'] = df_pheno_temp.index.values
        df_pheno_temp[['FD_SPID', phenotype_name]].to_csv(pheno_file, sep='\t', header=False, na_rep='NA', index=True)
        del df_pheno_temp['FD_SPID']

        # perform jackknife sampling
        jackknife_estimates_df = pd.DataFrame(columns=grms)
        sum_of_estimates_list = []
        indices = np.arange(df_pheno_temp.shape[0])
        keep_file = command_args.output_dir + str(phenotype_name) + '_keep.list'
#         for i in np.random.choice(range(df_pheno_temp.shape[0]), command_args.jackknife_iterations, replace=False): #TODO: think if we want to limit the number of jackknife rounds
        for i in range(df_pheno_temp.shape[0]):
            # take one sample out TODO: maybe should be at random
            sampls2keep = df_pheno_temp.reset_index().iloc[indices != i, 0]
            pd.concat((sampls2keep, sampls2keep), axis=1).to_csv(keep_file, sep='\t', header=False, na_rep='NA', index=False)
            cmdLine = [GCTA_EXE, '--mgrm-gz', command_args.output_dir + '/multi_grm.txt', '--reml',
                       '--pheno', pheno_file, '--reml-est-fix', '--thread-num', '1', '--keep', keep_file,
                       '--qcovar', covariates_file]
            proc = subprocess.Popen(cmdLine, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            stdout, stderr = proc.communicate()
            if (stderr is not None):
                print ('GCTA reml error:')
                print stderr
                raise Exception()
            sum_of_estimates = None
            intercept, fixed_effects, estimates, sum_of_estimates = _parse_GCTA_results_MC(stdout, phenotype_name)
            if sum_of_estimates is None:
                print (stdout)
                continue
            jackknife_estimates_df.loc[jackknife_estimates_df.shape[0]] = estimates
            sum_of_estimates_list.append(sum_of_estimates)

        # run GCTA over all samples
        cmdLine = [GCTA_EXE, '--mgrm-gz', command_args.output_dir + '/multi_grm.txt', '--reml',
           '--pheno', pheno_file, '--reml-est-fix', '--thread-num', '1', '--qcovar', covariates_file]
        proc = subprocess.Popen(cmdLine, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        stdout, stderr = proc.communicate()
        intercept, fixed_effects, estimates, sum_of_estimates = _parse_GCTA_results_MC(stdout, phenotype_name)
        if sum_of_estimates is None:
            print (stdout)
            os.remove(keep_file)
            os.remove(pheno_file)
            continue
        jackknife_estimates_df.loc['B'] = estimates

        for col in jackknife_estimates_df.columns:
            jackknife_estimates_df.loc['SE', col] = compute_SE_using_jackknife(jackknife_estimates_df[col].dropna().values)
            jackknife_estimates_df.loc['P value', col] = compute_pval_using_jackknife(jackknife_estimates_df.loc['B', col], jackknife_estimates_df.loc['SE', col])
            jackknife_estimates_df.loc['CI', col] = compute_CI_using_jackknife(jackknife_estimates_df.loc['B', col], jackknife_estimates_df.loc['SE', col])
        sum_of_estimates_SE = compute_SE_using_jackknife(sum_of_estimates_list)
        sum_of_estimates_pval = compute_pval_using_jackknife(sum_of_estimates, sum_of_estimates_SE)
        sum_of_estimates_CI = compute_CI_using_jackknife(sum_of_estimates, sum_of_estimates_SE)

        line =  '\t'.join([str(c) for c in [
                         phenotype_name, df_pheno_temp.shape[0], '%0.3f'%(intercept), len(sum_of_estimates_list)] +
                         ['%0.3g'%(item) for sublist in [[jackknife_estimates_df.loc[s, f] for s in ['B', 'CI', 'P value', 'SE']] for f in jackknife_estimates_df.columns] for item in sublist] +
                          ['%0.3g'%(item) for item in [sum_of_estimates, sum_of_estimates_CI, sum_of_estimates_pval, sum_of_estimates_SE]] +
                          ['%0.3g'%(f) for f in fixed_effects]]) + '\n'
        with open(command_args.output_dir + 'results' + str(idx), 'a') as handle:
            handle.write(line)
        sys.stdout.flush()

        try:
            os.remove(keep_file)
            os.remove(pheno_file)
        except:
            print ('Couldnt remove temp files for some reason: ' + str(phenotype_name))
    return

def compute_SE_using_jackknife(arr, k=1):
    return np.sqrt(np.var(arr, ddof=0) * (len(arr)-1))

def compute_pval_using_jackknife(b_i, se_i):
    statistic = b_i / se_i
    pvalue = stats.chi2(1).sf(statistic**2) / 2.0
    return pvalue

def compute_CI_using_jackknife(b_i, se_i):
    return ' - '.join([str(b_i - 1.96 * se_i), str(b_i + 1.96 * se_i)])

def create_multi_grm_file(command_args):
    with open(command_args.output_dir + '/multi_grm.txt', 'w') as handle:
        for f in command_args.random:
            if len(f.split('-')) > 1:
                mc_grms = os.listdir(GRM_DIR + '-'.join(command_args.samples) + '/' + f + '/')
                mc_grms = [mc for mc in mc_grms if mc.endswith('.grm.gz')]
                for mc in mc_grms:
                    handle.write(GRM_DIR + '-'.join(command_args.samples) + '/' + f + '/' + mc.split('.grm.gz')[0] + '\n')
            else:
                handle.write(GRM_DIR + '-'.join(command_args.samples) + '/' + f + '/' + f + '\n')
    return

def merge_temp_res_files(command_args):
    temp_res_files = os.listdir(command_args.output_dir)
    temp_res_files = [f for f in temp_res_files if f.startswith('results')]
    temp_res_files.sort(key=lambda x: int(x.split('results')[1]))
    temp_res_files = [command_args.output_dir + f for f in temp_res_files]
    final_results = pd.DataFrame()
    for f in temp_res_files:
        temp_res = pd.read_csv(f, sep='\t')
        final_results = pd.concat((final_results, temp_res))
        os.remove(f)
    final_results.to_csv(command_args.output_dir + 'LMM_results.txt', sep='\t')
    return


def upload_lmm_per_n_phenotypes(q, command_args):

    # read samples FD_SPID
    samples = get_FD_SPID_list(command_args)
    # get phenotype dataframe
    df_pheno = get_phenotype_df(command_args).loc[samples].dropna(how='all')

    # build the fixed effect matrix
    df_covariates = get_features(command_args, command_args.fixed).loc[samples].dropna(how='all')
    # split into single and multi-component
    if len(command_args.random) == 1 and len(command_args.random[0].split('-')) == 1:
        # run as single component
        q.waitforresults([q.method(estimate_h2_single_kernel, (command_args, df_pheno, df_covariates))])
        return
    else:
        create_multi_grm_file(command_args)
#         pheno_file = os.path.join(command_args.output_dir, 'temp.phe')
        covariates_file = os.path.join(command_args.output_dir, 'temp.cov')
        df_covariates['FD_SPID'] = df_covariates.index.values
        df_covariates = df_covariates[['FD_SPID'] + [c for c in df_covariates.columns if c != 'FD_SPID']]
        df_covariates.to_csv(covariates_file, sep='\t', header=False, na_rep='NA')

        waiton = []
        # upload n phenotypes per job
        for idx in range(0, df_pheno.shape[1], command_args.n_phenotypes_per_job):
            waiton.append(q.method(estimate_h2_multi_kernel, (command_args, df_pheno.iloc[:, idx:idx+command_args.n_phenotypes_per_job], df_covariates.columns, idx)))
        res = q.waitforresults(waiton)
        # merge the temp results files
        merge_temp_res_files(command_args)
    return

def _convert_comma_separated_to_list(s):
    return s.split(',')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('phenotype', help='Which phenotype/s to take as y. (Metabolomics_raw, Metabolomics_normed, BMI, Cohort, etc.)', type=str)
    parser.add_argument('samples', help='What samples to use. (ACS, MAR17, MAY18, etc.)', type=str)
    parser.add_argument('-output_dir', help='Path to output directory', type=str, default=LMM_DIR)
    parser.add_argument('-random', help='What are the random effects. (IGC, IGC-COG, MPA_species, MPA_species-phyla, KEGGBacGenes, Diet, etc.) separated by comma', type=str, default='IGC')
    parser.add_argument('-fixed', help='What are the fixed effects. (Age, Gender, BMI, Nextera, etc.) separated by comma', type=str, default='Age,Gender,Nextera')
    parser.add_argument('-n_phenotypes_per_job', help='Number of phenotypes per job', type=int, default=50)
    parser.add_argument('-use_quantile_normalization', help='Whether to use quantile normalization over microbiome data', type=bool, default=False)
    parser.add_argument('-jackknife_iterations', help='Number of Jackknife iterations', type=int, default=100)
    parser.add_argument('-output_name', help='Specify an output directory name', type=str, default=None)


#     parser.add_argument('-covariates', help='Which covariates to consider, separated by comma', type=str, default='Age,Gender')
#     parser.add_argument('-prevalence', help='In case of a case-control trait, what is the prevalence in the general population', type=float, default=None)
    parser.add_argument('-fiesta_iterations', help='Number of iterations to be made by FIESTA', type=int, default=100)

    command_args = parser.parse_args()
    make_dir_if_not_exists(command_args.output_dir)

    if command_args.output_name is not None:
        command_args.output_dir += '/' + command_args.output_name + '/'
        make_dir_if_not_exists(command_args.output_dir)

    if command_args.phenotype not in LEGAL_PHENOTYPES:
        print ('phenotype currently supported are: ' + ','.join(LEGAL_PHENOTYPES))
        return

    command_args.random = _convert_comma_separated_to_list(command_args.random)
    for grm in command_args.random:
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
    command_args.output_dir += '/' + command_args.phenotype + '/'
    make_dir_if_not_exists(command_args.output_dir)
    command_args.output_dir += '/' + '+'.join(command_args.random) + '/'
    make_dir_if_not_exists(command_args.output_dir)


    with open(command_args.output_dir + '/args' + str(datetime.now()), 'w') as handle:
        for arg in vars(command_args):
            handle.write(str(arg) + '\t' + str(getattr(command_args, arg)) + '\n')

    command_args.fixed = _convert_comma_separated_to_list(command_args.fixed)

    with qp(jobname = 'LMMs', q=['himem7.q'], mem_def = '1G', trds_def = 1, tryrerun=False,
#     with fakeqp(jobname = 'LMMs', q=['himem7.q'], mem_def = '1G', trds_def = 1, tryrerun=False,
        max_u = 220, delay_batch=15) as q:
        os.chdir("/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/temp_q_dir/")
        q.startpermanentrun()
        upload_lmm_per_n_phenotypes(q, command_args)

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



pheno_file = os.path.join(command_args.output_dir, 'temp.phe')
covariates_file = os.path.join(command_args.output_dir, 'temp.cov')
# saving the covariates
df_covariates.to_csv(covariates_file, sep='\t', header=False, na_rep='NA')
# saving the phenotype
df_pheno_temp[['FD_SPID', phenotype_name]].to_csv(pheno_file, sep='\t', header=False, na_rep='NA', index=True)

# running gcta
cmdLine = [GCTA_EXE, '--grm-gz', grm_file, '--reml', '--pheno', pheno_file,
           '--qcovar', covariates_file, '--reml-est-fix', '--thread-num', '1']
proc = subprocess.Popen(cmdLine, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
stdout, stderr = proc.communicate()
# the results will be in stdout

# load the grm from before (df_K)
cov_mat = Utils.Load(
    GRM_DIR + '-'.join(command_args.samples) + '/' + command_args.random[0] + '/' + command_args.random[0] + '.dat')
df_K = pd.DataFrame(cov_mat.values, index=cov_mat.index, columns=cov_mat.index)
df_K['IID'] = df_K.index
del df_K['IID']

#compute a p-value with RL-SKAT
df_merged = df_K.merge(df_covariates.dropna(), left_index=True, right_index=True).copy()
df_merged = df_merged.merge(df_pheno_temp[[phenotype_name]], left_index=True, right_index=True)
df_merged = df_merged.loc[df_merged[phenotype_name].notnull()]
covariates = df_merged[[c for c in df_covariates.columns if c != 'FD_SPID']].copy()
covariates['intercept'] = 1.0
covariates = covariates.values
pheno_vec = df_merged[phenotype_name].values
df_K_subset = df_K.loc[df_K.index.isin(df_merged.index), df_K.columns.isin(df_merged.index)]
df_K_subset = df_K_subset.loc[df_merged.index, df_merged.index]
assert (df_K_subset.index == df_K_subset.columns).all()
assert (df_K_subset.shape[0] == df_merged.shape[0])
assert (df_K_subset.index.isin(df_merged.index).all())
assert (df_K_subset.index == df_merged.index).all()
kernel = df_K_subset.values
pvalues = rl_skat.RL_SKAT_Full_Kernel(kernel, covariates, add_intercept=False).test(np.row_stack(pheno_vec))
#compute CIs with FIESTA
s,U = la.eigh(kernel)
ind = np.argsort(s)[::-1]; s=s[ind]; U=U[:,ind]
# take the h2_reml from the gcta results
CI = fiesta_lib.calculate_cis_general(np.array([h2_reml]), s, U, covariates,
                                      iterations=1000, alpha=0.05,
                                      tau=0.4, use_convergence_criterion=False,
                                      use_progress_bar=False) # iterations=1000