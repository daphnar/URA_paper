import pandas as pd
import numpy as np
import os
import glob
from scipy.stats import pearsonr
import statsmodels.formula.api as sm
import matplotlib.pyplot as plt
coverageFilesPath='/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Analyses/MetaPhlan100Uniq/bacteria_tmp/data'
translateFile='/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/databases/MetaPhlan/'+\
              'MetaPhlanNamesWithoutMissingGenomes.csv'#MetaPhlanNames.csv'
outputPath='/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Analyses/'+ \
           'MetaPhlan100Uniq/bacteria_tmp/abundances_by_sample'
MPApath='/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/databases/MetaPhlan/MPASpid.dat'
plotPath='/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Analyses/'+\
    'MetaPhlan100Uniq/bacteria_tmp/compare_to_metaphlan'
translate_df=pd.read_csv(translateFile).set_index('group')
MPA=pd.read_pickle(MPApath).reset_index()
drop_multiStrainGroups=False
plot=False
metaphlanThreshold=5e-5
def transformSpecies(species):
    res='s__'+species
    return res.lower()
MPA['Tax']=MPA['Tax'].apply(transformSpecies)
MPA=MPA[MPA['TaxLevel']=='s'].set_index('Tax')
corrs = {}
bacteriaCoverage=[]
methods=['leftmean_abnd','meannon0cov2_abnd','meannon0cov3_abnd','mean_abnd']
for method in methods:
    corrs[method]=[]
    if not os.path.exists(os.path.join(plotPath,method)):
        os.makedirs(os.path.join(plotPath,method))
for coveragefile in glob.glob(os.path.join(coverageFilesPath,'s__*')):
    groupsDic={}
    coverage_df=pd.read_csv(coveragefile,header=None)
    coverage_df['bac']="_".join(os.path.basename(coveragefile).replace('.csv','').split('_')[:-2])
    coverage_df['coverage_type']="_".join(os.path.basename(coveragefile).replace('.csv','').split('_')[-2:])
    coverage_df=coverage_df.rename(index=str,columns={0:'sample',1:'coverage'})
    bacteriaCoverage.append(coverage_df)
bacteriaCoverage_df=pd.concat(bacteriaCoverage)
FDs=set(bacteriaCoverage_df['sample'].values)
for sample in FDs:

    sample_coverage_df=bacteriaCoverage_df[bacteriaCoverage_df['sample']==sample]
    sample_coverage_df['species']=sample_coverage_df['bac'].apply(lambda x: x.split('_GCF')[0])

    for coverage_method in methods:
        currentMethod=sample_coverage_df['coverage_type']==coverage_method
        sample_coverage_df.loc[currentMethod,'relative_abundance']=sample_coverage_df.loc[currentMethod,'coverage']/sample_coverage_df.loc[currentMethod,'coverage'].sum()
    if drop_multiStrainGroups:
        sample_coverage_df.drop_duplicates(subset=['species','coverage_type'],keep=False,inplace=True)
    groupCoverage_df=sample_coverage_df.groupby(['species','coverage_type']).sum()['relative_abundance']#.reset_index()
    groupCoverage_df.reset_index().to_csv(os.path.join(outputPath, sample + '.csv'))

    for coverage_method in methods:
        species=groupCoverage_df.index.levels[0]
        method_groupCoverage_df=groupCoverage_df.reset_index()
        method_groupCoverage_df=method_groupCoverage_df[method_groupCoverage_df['coverage_type']==coverage_method]
        method_groupCoverage_df=method_groupCoverage_df.set_index('species')
        MPAshared=list(set(MPA.loc[species,sample].dropna().index))
        regression_df=pd.DataFrame(index=MPAshared,data={'MPA':MPA.loc[MPAshared,sample].values,
                              'Unicorn':method_groupCoverage_df.loc[MPAshared,'relative_abundance'].values})#\
                    #translate_df.reset_index().set_index('species').loc[shared,'group']})
        regression_df[regression_df<metaphlanThreshold]=metaphlanThreshold
        regression_df=np.log10(regression_df)
        cor, pval = pearsonr(regression_df['MPA'].values, regression_df['Unicorn'].values)
        corrs[coverage_method].append(cor)
        if plot:
            stat=sm.ols(formula='MPA ~ -1+Unicorn', data=regression_df).fit()
            beta=stat.params['Unicorn']
            plt.figure()
            plt.scatter(regression_df['Unicorn'].values,regression_df['MPA'].values,c='b')
            plt.plot(regression_df['Unicorn'].values,regression_df['Unicorn'].values*beta,c='r')
            plt.xlabel('Unicorn relative abundance')
            plt.ylabel('MetaPhlan relative abundance')
            plt.title(sample+' b=%.2f'%beta)
            plt.savefig(os.path.join(plotPath,coverage_method,sample+'.png'))
            plt.close('all')
for method in methods:
    print method, ' ', np.nanmean(corrs[method])
