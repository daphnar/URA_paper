import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import os
import matplotlib.gridspec as gridspec
import nature_guidline_utils

sns.set_style("ticks", {'axes.edgecolor': 'black'})
pd.set_option('display.width', 1000)
np.set_printoptions(precision=4, linewidth=200)
FIGURES_DIR = '/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/figures'
colors=sns.color_palette('RdBu')
colors= [colors[0],colors[-1]]
colors_rgb=sns.color_palette('colorblind', 5)

params = {
    'axes.labelsize': 12,
    'font.size': 12,
    'legend.fontsize': 12,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'figure.dpi': 300,
    'axes.linewidth': 0.5}

fontsize = 12
plt.rcParams.update(params)
# figsize=(183/10/2.54,247/10/2.54)
figsize=(183/10/2.54,185/10/2.54)

plt.figure(figsize=(nature_guidline_utils.two_columns(),
                    nature_guidline_utils.full_page()*(1./4)), dpi=300)

outer_grid = gridspec.GridSpec(1,2)
outer_grid.update(wspace=0.3)

a_ax = plt.subplot(outer_grid[0,0])
b_ax = plt.subplot(outer_grid[0,1])


plt.sca(b_ax)
def_data= pd.read_csv(os.path.join(FIGURES_DIR, 'Figures - e_coli_pathways.csv'))
peach =  (205./255,112./255,106./255)
violet =  (106./255,111./255,205./255)

plt.title('HbA1C%',fontsize=10)
plt.scatter(range(def_data.shape[0]),(-1)*np.log10(def_data['pval']),c='gray',
            edgecolor='',alpha=0.5,marker='.')
annotate1=plt.scatter(def_data[def_data['is_type4']==1].index.values,(-1)*np.log10(def_data[def_data['is_type4']==1]['pval']),c=colors_rgb[0],#[1],
            edgecolor='',marker='.')
annotate2=plt.scatter(def_data[def_data['is_colicin']==1].index.values,(-1)*np.log10(def_data[def_data['is_colicin']==1]['pval']),c=colors_rgb[1],#[4],
            edgecolor='',marker='.')
annotate3=plt.scatter(def_data[def_data['is_ABC']==1].index.values,(-1)*np.log10(def_data[def_data['is_ABC']==1]['pval']),c=colors_rgb[2],#[6],
            edgecolor='',marker='.')
plt.xlabel('SGB 10068: E.coli window num.')
plt.ylabel('P-value (-log10)')
b_ax.spines['right'].set_visible(False)
b_ax.spines['top'].set_visible(False)
b_ax.tick_params(top=False, right=False)
for item in ([b_ax.xaxis.label, b_ax.yaxis.label]):
    item.set_fontsize(10)

plt.xlim([0,14000])
plt.ylim([0,35])
plt.xticks([0,3500,7000,10500,14000])
plt.legend([annotate1,annotate2,annotate3],['T4SS','Colicin','ABC transporter'],fontsize=10,loc=2,frameon=False)

plt.sca(a_ax)
def_data = pd.read_csv(os.path.join(FIGURES_DIR, 'Figures - kpan pathways.csv'))
plt.title('HbA1C%',fontsize=10)
plt.scatter(range(def_data.shape[0]),(-1)*np.log10(def_data['pval']),c='gray',
            edgecolor='',alpha=0.5,marker='.')
annotate1=plt.scatter(def_data[def_data['is_type_4']==1].index.values,(-1)*np.log10(def_data[def_data['is_type_4']==1]['pval']),c=colors_rgb[0],#[1],
            edgecolor='',marker='.')
# annotate2=plt.scatter(def_data[def_data['is_fb']==1].index.values,(-1)*np.log10(def_data[def_data['is_fb']==1]['pval']),c=colors_rgb[3],#[3],
#             edgecolor='',marker='.')
annotate3=plt.scatter(def_data[def_data['is_aero']==1].index.values,(-1)*np.log10(def_data[def_data['is_aero']==1]['pval']),c=colors_rgb[4],#[0],
            edgecolor='',marker='.')
plt.xlabel('SGB 10115: K.pneumoniae window num.')
plt.ylabel('P-value (-log10)')
a_ax.spines['right'].set_visible(False)
a_ax.spines['top'].set_visible(False)
a_ax.tick_params(top=False, right=False)
for item in ([a_ax.xaxis.label, a_ax.yaxis.label]):
    item.set_fontsize(10)
plt.legend([annotate1,annotate3],['T4SS','Aerobactin assembly'],#annotate2,'Fimbriae assembly'
           bbox_to_anchor=(1.1, 1),fontsize=10,loc=1,frameon=False)

plt.xlim([0,8000])
plt.ylim([0,10])
plt.xticks([0,2000,4000,6000,8000])
plt.savefig(os.path.join(FIGURES_DIR, 'figure6.pdf'), bbox_inches='tight', format='pdf')
plt.savefig(os.path.join(FIGURES_DIR, 'figure6.png'), bbox_inches='tight', format='png')