from SHAP import summary_plot
import os
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib as mpl
import nature_guidline_utils
import pandas as pd
import numpy as np
# import plotp_fixed as plotp
import matplotlib.cm as cm
import seaborn as sns
from scipy.stats import gaussian_kde,spearmanr
from sklearn.metrics import r2_score
#FIGURES_DIR = '/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/figures'
FIGURES_DIR = '/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Cohort_Paper/revision_Analyses/figures'

num_bacs = 3

plt.figure(figsize=(nature_guidline_utils.two_columns(),
                    nature_guidline_utils.full_page()*1/2), dpi=300)
params = {
    'axes.labelsize': 10,
    'font.size': 10,
    'legend.fontsize': 10,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'figure.dpi': 300,
    'axes.linewidth': 0.5,
}
plt.rcParams.update(params)

outer_grid = gridspec.GridSpec(2, 1, height_ratios=[0.5,0.5])
outer_grid.update(hspace=0.7)
abc_grid = gridspec.GridSpecFromSubplotSpec(1,4, outer_grid[0,:],width_ratios=[ 0.33,0.33,0.33,0.02],wspace=0.5)
def_grid = gridspec.GridSpecFromSubplotSpec(1,4, outer_grid[1,:],width_ratios=[ 0.33,0.33,0.33,0.02],wspace=0.5)

outer_grid.update(wspace=1)
color= "coolwarm"
colors=sns.color_palette('RdBu')
colors= [colors[0],colors[-1]]

def colorbar(ax):
    m = cm.ScalarMappable(cmap=color)
    m.set_array([0, 1])
    cb = plt.colorbar(m, ticks=[0, 1],cax = ax)#, aspect=1000)
    cb.set_ticklabels(['P=1',r'P<$10^{-20}}$'])
    cb.set_label('P value (Spearman correlation)', size=10, labelpad=0)
    cb.ax.tick_params(labelsize=10, length=0)
    # cb.set_alpha(1)
    cb.outline.set_visible(False)
    # bbox = cb.ax.get_window_extent().transformed(plt.gcf().dpi_scale_trans.inverted())
    # cb.ax.set_aspect((bbox.height - 0.9) * 20)

def keepname(name,return_sgb):
    name = name.split('|')
    sgb = name[-1][6:]
    corrected_names={'6936':'V. atypica',
                     '4826': 'Blautia sp',
                     '17244': 'B. adolescentis',
                     '10068': 'E. coli',
                     '10064': 'E. marmotae',
                     '10115': 'K. pneumoniae',
                     '4705': 'Clostridium sp',
                     '15369': 'Faecalibacterium sp',
                     '4964' : 'Unknown Eubacteriaceae',
                     '14311' : 'Unknown Clostridiaceae'}
                     # '8069' : 'S. parasanguinis'}
    if sgb in corrected_names.keys():
        name =corrected_names[sgb]
    else:
        name = "?"
    if return_sgb:
        return 'SGB %s: %s'%(sgb,name)
    else:
        return name

def keepnameOLD(name, return_sgb):
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
def_data= pd.read_csv(os.path.join(FIGURES_DIR, 'Figures - il_vs_us_sgb_pvals.csv'))
d_ax = plt.subplot(abc_grid[0,0])
plt.sca(d_ax)
d_ax.set_zorder(1)
plt.text(-0.15, 1.1, 'a', ha='center', va='center', transform=d_ax.transAxes, fontsize=16)
plt.title('Age',fontsize=10)
def_data= pd.read_csv(os.path.join(FIGURES_DIR, 'Figures - il_vs_us_sgb_pvals.csv'))
def_data=def_data.set_index('species')
pheno_data=def_data[def_data['pheno']=='age']
pheno_data=pheno_data.sort_values('spear_pval')
bacteria = pheno_data.index.values
plt.scatter(pheno_data.loc[bacteria,'spear'],pheno_data.loc[bacteria,'spear_us'],# c=colors[1],
            c=np.log10(pheno_data['spear_pval']), cmap=plt.cm.get_cmap('RdBu'),
            norm=mpl.colors.Normalize(vmin=-20, vmax=0,clip=True),
            edgecolor='',alpha=0.5)
r=spearmanr(pheno_data.loc[bacteria,'spear'],pheno_data.loc[bacteria,'spear_us'])
d_ax.annotate('R=%.2f\n'
              r'P<$10^{%d}}$'%(r[0],np.log10(r[1])-1),(-.18,0.12),fontsize=8)

plt.xlabel('Train-IL\nSpearmann correlation')
plt.ylabel('Test2-US\nSpearmann correlation')
for i, txt in enumerate([keepname(sgb, False) for sgb in bacteria[:num_bacs+2]]):
    if txt=='?':
        continue
    # elif txt=='V. atypica':
    #     t = plt.gca().annotate(txt,
    #              (pheno_data.loc[bacteria, 'spear'].iloc[i],# - 0.02,
    #     pheno_data.loc[bacteria, 'spear_us'].iloc[i]),# - 0.035),
    #               fontsize=8)
        # t.set_bbox(dict(facecolor='red', alpha=0.5, edgecolor='red'))
    elif txt == 'Blautia sp':
        plt.gca().annotate(txt,
                           (pheno_data.loc[bacteria, 'spear'].iloc[i] + 0.015,
                            pheno_data.loc[bacteria, 'spear_us'].iloc[i] - 0.04),  # - 0.01),
                           fontsize=8)
    else:
        plt.gca().annotate(txt,
                       (pheno_data.loc[bacteria, 'spear'].iloc[i] + 0.015,
                        pheno_data.loc[bacteria, 'spear_us'].iloc[i]),# - 0.01),
                       fontsize=8)
    plt.scatter(pheno_data.loc[bacteria[i], 'spear'], pheno_data.loc[bacteria[i], 'spear_us'], c=colors[0],
                edgecolors='k')

d_ax.tick_params(axis='both', which='major', pad=3)
# plt.xlim(-0.20,0.20)
# plt.ylim(-0.20,0.20)
plt.xticks([-0.20,0,0.20])
d_ax.set_xticklabels(["-0.2","0","0.2"])
plt.yticks([-0.2,0,0.2])
d_ax.set_yticklabels(["-0.2","0","0.2"])
e_ax = plt.subplot(abc_grid[0,1])
plt.sca(e_ax)
plt.text(-0.15, 1.1, 'b', ha='center', va='center', transform=e_ax.transAxes, fontsize=16)
plt.title('HbA1C%',fontsize=10)
def_data= pd.read_csv(os.path.join(FIGURES_DIR, 'Figures - il_vs_us_sgb_pvals.csv'))
def_data=def_data.set_index('species')
pheno_data=def_data[def_data['pheno']=='hba1c']
pheno_data=pheno_data.sort_values('spear_pval')
bacteria = pheno_data.index.values
plt.scatter(pheno_data.loc[bacteria,'spear'],pheno_data.loc[bacteria,'spear_us'],
            c=np.log10(pheno_data['spear_pval']), cmap=plt.cm.get_cmap('RdBu'),
            norm=mpl.colors.Normalize(vmin=-20, vmax=0,clip=True),
            edgecolor='',alpha=0.5)
keep_rev=pheno_data.loc[bacteria,'spear_us'].dropna().index
#r=spearmanr(pheno_data.loc[bacteria,'spear'].dropna(),pheno_data.loc[bacteria,'spear_us'])
r=spearmanr(pheno_data.loc[keep_rev,'spear'].dropna(),pheno_data.loc[keep_rev,'spear_us'])
e_ax.annotate('R=%.2f\n'
              r'P<$10^{%d}$'%(r[0],int(np.log10(r[1])-1)),(-.18,0.17),fontsize=8)
plt.xlabel('Train-IL\nSpearmann correlation')
for i, txt in enumerate([keepname(sgb, False) for sgb in bacteria[:num_bacs]]):
    plt.scatter(pheno_data.loc[bacteria[i], 'spear'], pheno_data.loc[bacteria[i], 'spear_us'], c=colors[0],
                edgecolors='k')
    # if txt=='K. pneumoniae' or
    #
    if txt=='E. marmotae':
        plt.gca().annotate(txt,
                           (pheno_data.loc[bacteria, 'spear'].iloc[i],# - 0.07,
                            pheno_data.loc[bacteria, 'spear_us'].iloc[i]-0.045),
                           fontsize=8)
    else:
        plt.gca().annotate(txt,
        (pheno_data.loc[bacteria, 'spear'].iloc[i] - 0.02,
        pheno_data.loc[bacteria, 'spear_us'].iloc[i] + 0.018),fontsize=8)
e_ax.tick_params(axis='both', which='major', pad=3)
# plt.xlim(-0.2,0.2)
# plt.ylim(-0.2,0.2)
plt.xticks([-0.2,0,0.25])
e_ax.set_xticklabels(["-0.2","0","0.25"])
plt.yticks([-0.2,0,0.25])
e_ax.set_yticklabels(["-0.2","0","0.25"])

f_ax = plt.subplot(abc_grid[0,2])
plt.sca(f_ax)
plt.text(-0.15, 1.1, 'c', ha='center', va='center', transform=f_ax.transAxes, fontsize=16)
plt.title('BMI',fontsize=10)
def_data= pd.read_csv(os.path.join(FIGURES_DIR, 'Figures - il_vs_us_sgb_pvals.csv'))
def_data=def_data.set_index('species')
pheno_data=def_data[def_data['pheno']=='bmi']
pheno_data=pheno_data.sort_values('spear_pval')
bacteria = pheno_data.index.values
plt.scatter(pheno_data.loc[bacteria,'spear'],pheno_data.loc[bacteria,'spear_us'],
            c=np.log10(pheno_data['spear_pval']), cmap=plt.cm.get_cmap('RdBu'),
            norm=mpl.colors.Normalize(vmin=-20, vmax=0,clip=True),
            edgecolor='',alpha=0.5)
r=spearmanr(pheno_data.loc[bacteria,'spear'],pheno_data.loc[bacteria,'spear_us'])
f_ax.annotate('R=%.2f\n'
              r'P<$10^{%d}}$'%(r[0],np.log10(r[1])-1),(-.23,0.12),fontsize=8)
plt.xlabel('Train-IL\nSpearmann correlation')
for i, txt in enumerate([keepname(sgb, False) for sgb in bacteria[:num_bacs]]):
    plt.scatter(pheno_data.loc[bacteria[i], 'spear'], pheno_data.loc[bacteria[i], 'spear_us'], c=colors[0],
                edgecolors='k')
    # plt.gca().annotate(txt, (pheno_data.loc[bacteria, 'spear'].iloc[i]+0.01, pheno_data.loc[bacteria, 'spear_us'].iloc[i]),
    #                    fontsize=8)
    if txt=='?':
        continue
    if txt=='Faecalibacterium sp':
        plt.gca().annotate(txt, (pheno_data.loc[bacteria, 'spear'].iloc[i] - 0.07,
        pheno_data.loc[bacteria, 'spear_us'].iloc[i]-0.04), fontsize=8)
    elif txt=='Clostridium sp':
        plt.gca().annotate(txt, (pheno_data.loc[bacteria, 'spear'].iloc[i] - 0.11,
        pheno_data.loc[bacteria, 'spear_us'].iloc[i] + 0.02), fontsize=8)
    else:
        plt.gca().annotate(txt,(pheno_data.loc[bacteria, 'spear'].iloc[i] - 0.07,
        pheno_data.loc[bacteria, 'spear_us'].iloc[i] - 0.04), fontsize=8)
f_ax.tick_params(axis='both', which='major', pad=3)
# plt.xlim(-0.25,0.2)
# plt.ylim(-0.25,0.2)
plt.xticks([-0.25,0,0.2])
f_ax.set_xticklabels(["-0.25","0","0.2"])
plt.yticks([-0.25,0,0.2])
f_ax.set_yticklabels(["-0.25","0","0.2"])
colorbar_ax = plt.subplot(abc_grid[0,3])
colorbar(colorbar_ax)
d_ax.spines['right'].set_visible(False)
d_ax.spines['top'].set_visible(False)
d_ax.tick_params(top=False, right=False)
for item in ([d_ax.xaxis.label, d_ax.yaxis.label]):
    item.set_fontsize(10)
e_ax.spines['right'].set_visible(False)
e_ax.spines['top'].set_visible(False)
e_ax.tick_params(top=False, right=False)
for item in ([e_ax.xaxis.label, e_ax.yaxis.label]):
    item.set_fontsize(10)
f_ax.spines['right'].set_visible(False)
f_ax.spines['top'].set_visible(False)
f_ax.tick_params(top=False, right=False)
for item in ([f_ax.xaxis.label, f_ax.yaxis.label]):
    item.set_fontsize(10)



def_data= pd.read_csv(os.path.join(FIGURES_DIR, 'Figures - il_vs_us_spearman_saturation.csv'))


d_ax = plt.subplot(def_grid[0,0])
plt.sca(d_ax)
plt.title('Age',fontsize=10)
plt.text(-0.15, 1.1, 'd', ha='center', va='center', transform=d_ax.transAxes, fontsize=16)
data = def_data[def_data['target']=='age']
d_ax.errorbar(data['size'],data['spearman_r'],yerr=data['std'],color=colors[1])
d_ax.spines['right'].set_visible(False)
d_ax.spines['top'].set_visible(False)
d_ax.tick_params(top=False, right=False)
plt.xticks([0,2000,4000])
plt.xlim([0,4000])
plt.ylim([0,0.6])
plt.yticks([0,0.3,0.6])
d_ax.set_yticklabels(['0','0.3','0.6'])
plt.xlabel('Cohort sizes (IL & US)')
plt.ylabel('Spearman correlation\nbetween cohorts')

e_ax = plt.subplot(def_grid[0,1])
plt.sca(e_ax)
plt.title('HbA1C%',fontsize=10)
plt.text(-0.15, 1.1, 'e', ha='center', va='center', transform=e_ax.transAxes, fontsize=16)
data = def_data[def_data['target']=='hba1c']
e_ax.errorbar(data['size'],data['spearman_r'],yerr=data['std'],color=colors[1])
e_ax.spines['right'].set_visible(False)
e_ax.spines['top'].set_visible(False)
e_ax.tick_params(top=False, right=False)
plt.xticks([0,1000,2000])
plt.xlim([0,2300])
plt.yticks([0,0.25,0.5])
e_ax.set_yticklabels(['0','0.25','0.5'])
plt.xlabel('Cohort sizes (IL & US)')

f_ax = plt.subplot(def_grid[0,2])
plt.sca(f_ax)
plt.title('BMI',fontsize=10)
plt.text(-0.15, 1.1, 'f', ha='center', va='center', transform=f_ax.transAxes, fontsize=16)
data = def_data[def_data['target']=='bmi']
f_ax.errorbar(data['size'],data['spearman_r'],yerr=data['std'],color=colors[1])
f_ax.spines['right'].set_visible(False)
f_ax.spines['top'].set_visible(False)
f_ax.tick_params(top=False, right=False)
plt.xticks([0,2000,4000])
plt.xlim([0,4000])
plt.ylim([0,0.7])
plt.xlabel('Cohort sizes (IL & US)')
plt.yticks([0.0,0.35,0.7])
f_ax.set_yticklabels(['0','0.35','0.7'])
plt.savefig(os.path.join(FIGURES_DIR, 'figure5.pdf'), bbox_inches='tight', format='pdf')
plt.savefig(os.path.join(FIGURES_DIR, 'figure5.png'), bbox_inches='tight', format='png')
