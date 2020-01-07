import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import os
import matplotlib.gridspec as gridspec
FIGURES_DIR = '/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/figures'
OUTPUT_DIR = os.path.join(FIGURES_DIR,"Abundances-IL_US")
params = {
    'axes.labelsize': 12,
    'font.size': 12,
    'legend.fontsize': 12,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'figure.dpi': 300,
    'axes.linewidth': 0.5}

fontsize = 12
NUM_BACTERIA='all'
plt.rcParams.update(params)
illustraitor_il=(34./255,181./255,115./255)
illustraitor_il_validation=(41./255,171./255,226./255)
illustraitor_us=(102./255,45./255,145./255)
colors_rgb = [illustraitor_il,illustraitor_il_validation,illustraitor_us]
us_sgb=pd.read_csv(os.path.join(FIGURES_DIR,'sgb__us.csv'),index_col=0).apply(np.log10)
il_sgb=pd.read_csv(os.path.join(FIGURES_DIR,'sgb__il.csv'),index_col=0).apply(np.log10)
il_validation_sgb=pd.read_csv(os.path.join(FIGURES_DIR,'sgb__il_validation.csv'),index_col=0).apply(np.log10)

def keepname(name,return_sgb):
    name=name.split('|')
    sgb = name[-1]
    species = name[-4]
    genus = name[-5]
    family = name[-6]
    order = name[-7]
    sclass = name[-8]
    phylum = name[-9]
    if species!='s__unknown' and species[-3:]!='_sp' and 'unclassified' not in species:
        ncbi=species[3:].split('_')
        if ncbi[1]== 'sp':
            pass
        ncbi = ncbi[0][0]+'.'+" ".join(ncbi[1:])+' (s)'
        return_val =  ncbi
    elif genus!='g__unknown' and 'unclassified' not in genus:
        return_val =  genus[3:].replace('_', ' ') + ' (g)'
    elif family!='f__unknown' and 'unclassified' not in family:
        return_val =  family[3:].replace('_', ' ') + ' (f)'
    elif order!='o__unknown' and 'unclassified' not in order:
        return_val =  order[3:].replace('_', ' ') + ' (o)'
    elif sclass!='c__unknown' and 'unclassified' not in sclass:
        return_val =  order[3:].replace('_', ' ') + ' (c)'
    elif phylum!='p__unknown' and 'unclassified' not in phylum:
        return_val =  order[3:].replace('_', ' ') + ' (p)'
        if return_val == 'Firmicutes unclassified (p)':
            return_val = 'unknown (p)'
    else:
        return_val =  'unknown (p)'
    if return_sgb:
        return_val = 'SGB ' + sgb[6:] + ': '+ return_val
    return return_val


phenotypes = ['hba1c','bmi','age']
labels = ['Train-IL','Test1-IL','Test2-US']

def_data= pd.read_csv(os.path.join(FIGURES_DIR, 'Figures - il_vs_us_sgb_pvals.csv'))
def_data=def_data.set_index('species')
for phenotype in phenotypes:
    pheno_data=def_data[def_data['pheno']==phenotype]
    pheno_data=pheno_data.sort_values('spear_pval')
    if NUM_BACTERIA!='all':
        bacteria = pheno_data.index.values[:NUM_BACTERIA]
    else:
        bacteria = pheno_data.index.values
    US_bacteria=us_sgb[bacteria]
    IL_bacteria=il_sgb[bacteria]
    num_US_bacteria = 1.0*(US_bacteria>-4).sum()/US_bacteria.shape[0]
    num_IL_bacteria = 1.0*(IL_bacteria > -4).sum()/IL_bacteria.shape[0]
    plt.xlabel('IL presence')
    plt.ylabel('US presence')
    plt.title(phenotype)
    plt.plot([0,1],[0,1],'k')
    plt.xlim([0,1])
    plt.ylim([0,1])
    plt.scatter(num_IL_bacteria,num_US_bacteria)
    for i, txt in enumerate([keepname(sgb,True) for sgb in bacteria]):
        plt.gca().annotate(txt, (num_IL_bacteria[i], num_US_bacteria[i]), fontsize=8)
    plt.savefig(os.path.join(OUTPUT_DIR, 'figure5-%s-presence_%s.png' %(phenotype,NUM_BACTERIA)), format='png',
                bbox_inches='tight')
    plt.close()
    plt.scatter(pheno_data.loc[bacteria,'spear'],pheno_data.loc[bacteria,'spear_us'])
    plt.xlabel('IL Spearmann')
    plt.ylabel('US Spearmann')
    # if NUM_BACTERIA!='all':
    if True:
        for i, txt in enumerate([keepname(sgb,True) for sgb in bacteria[:10]]):
            plt.gca().annotate(txt, (pheno_data.loc[bacteria,'spear'].iloc[i], pheno_data.loc[bacteria,'spear_us'].iloc[i]),
                               fontsize=8)
    plt.title(phenotype)
    plt.savefig(os.path.join(OUTPUT_DIR, 'figure5-%s-spearmann_%s.png' %(phenotype,NUM_BACTERIA)), format='png',
                bbox_inches='tight')
    plt.close()
    # for sgb in bacteria:
    #     IL=il_sgb[sgb]
    #     IL_validation=il_validation_sgb[sgb]
    #     US=us_sgb[sgb]
    #     sns.distplot(IL,norm_hist=True, kde = True,color=colors_rgb[0],label=labels[0],kde_kws = {'linewidth': 2})
    #     sns.distplot(IL_validation,norm_hist=True, kde = True,color=colors_rgb[1],label=labels[1],kde_kws = {'linewidth': 2})
    #     sns.distplot(US,norm_hist=True, kde = True,color=colors_rgb[2],label=labels[2],kde_kws = {'linewidth': 2})
    #     plt.xlabel('Log abundance')
    #     plt.ylabel('Frequency', fontsize=fontsize - 2)
    #     plt.title("%s - %s IL p-value %s, US p-value %s"%(phenotype,keepname(sgb,True),
    #                 pheno_data.loc[sgb,'spear_pval'],pheno_data.loc[sgb,'spear_pval_us']))
    #     plt.legend(labels, loc='upper right', bbox_to_anchor=(1, 1), frameon=False, fontsize=fontsize - 2)
    #     plt.savefig(os.path.join(OUTPUT_DIR, 'figure5-%s-%s.png'%(phenotype,keepname(sgb,False))), format='png', bbox_inches='tight')
    #     plt.close()
