import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import os
import matplotlib as mpl
from URA_paper.Figures import nature_guidline_utils
from mne.stats import fdr_correction
from dark import cmap_map
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
from scipy.stats import gaussian_kde
from sklearn.metrics import r2_score
sns.set_style("ticks", {'axes.edgecolor': 'black'})
pd.set_option('display.width', 1000)
np.set_printoptions(precision=4, linewidth=200)
FIGURES_DIR = '/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Cohort_Paper/revision_Analyses/figures'
params = {
    'axes.labelsize': 10,
    'font.size': 10,
    'legend.fontsize': 10,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'figure.dpi': 300,
    'axes.linewidth': 0.5,
}
fontsize = 10
plt.rcParams.update(params)
fig = plt.figure(figsize=(nature_guidline_utils.two_columns(),
                          nature_guidline_utils.full_page()*0.25), dpi=300)  # m2inch(165)

import matplotlib.gridspec as gridspec

# col_list = ["warm grey", "gunmetal", "dusky blue",
#             "cool blue", "deep teal",
#             "viridian", "turtle green", "twilight blue"][::-1]
#
# col_list = ["gunmetal", "dusky blue", "viridian"][::-1]
#
# colors_rgb = sns.xkcd_palette(col_list)
blue=(43./255,74./255,130./255)
green=(47./255,142./255,52./255)
pink=(191./255,64./255,140./255)
lightBlue=(42./255,107./255,126./255)
peach =  (205./255,112./255,106./255)
violet =  (106./255,111./255,205./255)
red = (54./255,23./255,69./255)
nothing = (0,0,0)
colors_rgb=sns.color_palette('Paired', 7)
cmap_paired=LinearSegmentedColormap.from_list(
        'Paired', colors_rgb, N=7)
new_cmap = cmap_map(lambda x: x*0.9, cmap_paired)
colors_rgb_maker = plt.get_cmap(new_cmap,7)
colors_rgb = [colors_rgb_maker(x) for x in np.linspace(0, 1, 7)]
two_colors = [violet,green]#colors_rgb[-2:]
three_colors = [colors_rgb[x] for x in [1,3,5]]
colors_rgb_age = [nothing]+colors_rgb[-3:]+[nothing]
colors_rgb_bmi_hba1c = colors_rgb[:4]+[colors_rgb[-1]]
# colors_rgb_age = [nothing,green,lightBlue,peach,nothing]
# colors_rgb_bmi_hba1c = [blue,red,pink,violet,peach]
# colors_rgb = [blue,red,pink,violet,green,lightBlue,peach]


outer_grid = gridspec.GridSpec(1,1)
ax__ab = outer_grid[0,0]
ab_grid = gridspec.GridSpecFromSubplotSpec(1,2, ax__ab, wspace=0.2, width_ratios=[17./21,4./21])
axa__quantitive_phenotypes = plt.subplot(ab_grid[0,0])
axa__quantitive_phenotypes.set_zorder(1)
axb__binary_phenotypes = plt.subplot(ab_grid[0,1])
axb__binary_phenotypes.set_zorder(1)





rename = {'age':'Age',
          'bmi':'BMI',
          'hba1c':'HbA1C%',
          'bowel_movement_frequency':'Bowel movement',
          'bt__fasting_glucose':'Glucose',
          'bt__fasting_triglycerides':'Triglycerides',
          'bt__hdl_cholesterol': 'HDL cholesterol',
          'height':'Height',
          'bt__cholesterol': 'Total cholesterol',
          'bt__protein':"Protein",
          'bt__albumin':"Albumin",
          'bt__sgot': "SGOT",
          'bt__sgpt': "SGPT",
          'bt__alkaline_phosphatase':"Alkaline phosphatase",
          'bt__tsh':'TSH ',
          'bt__inr':'INR  ',
          'T2D':'Type 2 diabetes',
          'currently_smokes':'Currently smokes',
          'ever_smoked':'Ever smoked',
          'type_2_diabetes': 'Type 2 diabetes',
          'gender': 'Gender'
          }
plt.sca(axa__quantitive_phenotypes)
plt.text(-.175, 1.1, 'a', ha='center', va='center', transform=axa__quantitive_phenotypes.transAxes, fontsize=16)
correlations_df = pd.read_csv(os.path.join(FIGURES_DIR,'Figures - prediction_other_phenotypes.csv')).set_index('target')
# phenotypes = correlations_df[['pearson','stdev']].dropna().sort_values(by='pearson',ascending=False).mul(100)
phenotypes = correlations_df[['pearson','pearson_linear','stdev','stdev_linear']].dropna().sort_values(by='pearson',ascending=False).mul(100)
phenotypes.columns = ['GBDT','Ridge','std_xgboost','std_ridge']
phenotypes.index=phenotypes.index.map(lambda x: rename[x] if x in rename else '')#phenotypes.index.str.replace('bt__',"")
phenotypes=phenotypes.loc[phenotypes.index[phenotypes.index !='']]
ind = np.array(range(len(phenotypes)))


# axa__quantitive_phenotypes.bar(ind,phenotypes['pearson'],yerr=phenotypes['stdev'], ecolor='black',#edgecolor='none',
#        zorder=1,color=colors_rgb[0],align='center')
estimation = phenotypes[['GBDT','Ridge']]
yerr = phenotypes[['std_xgboost','std_ridge']].T.values
estimation.plot.bar(yerr=yerr,ax = axa__quantitive_phenotypes,
                              zorder=1,#edgecolor='none',
                              color=two_colors,legend=True,
                              width=0.8)
axa__quantitive_phenotypes.legend(frameon=False)#loc=10,

axa__quantitive_phenotypes.set_xlim(ind.min()-0.5,ind.max()+0.5)
axa__quantitive_phenotypes.spines['bottom'].set_position('zero')
axa__quantitive_phenotypes.spines['left'].set_bounds(0, 31)
plt.yticks([0,10,20,30])
plt.xticks(ind,phenotypes.index,rotation=45,ha='right')
axa__quantitive_phenotypes.tick_params(top='off',right='off',pad=0,labelsize=fontsize)
axa__quantitive_phenotypes.yaxis.set_ticks_position('left')
axa__quantitive_phenotypes.xaxis.set_ticks_position('bottom')
sns.despine(ax=axa__quantitive_phenotypes)
plt.ylabel('$R^{2}$ (%)')
plt.xlabel('')
# plt.ylim(0,31)


plt.sca(axb__binary_phenotypes)
correlations_df = pd.read_csv(os.path.join(FIGURES_DIR,'Figures - classification.csv')).set_index('target')
phenotypes = correlations_df[['auc','auc_linear','stdev','stdev_linear']].dropna().sort_values(by='auc',ascending=False)
phenotypes.index=phenotypes.index.map(lambda x: rename[x] if x in rename else '')#phenotypes.index.str.replace('bt__',"")
plt.text(-.6, 1.1, 'b', ha='center', va='center', transform=axb__binary_phenotypes.transAxes, fontsize=16)
ind = range(len(phenotypes))
estimation = phenotypes[['auc','auc_linear']]
yerr = phenotypes[['stdev','stdev_linear']].T.values
estimation.plot.bar(yerr=yerr,ax = axb__binary_phenotypes,
                              zorder=1,#edgecolor='none',
                              color=two_colors,legend=False,
                              width=0.8)

# axb__binary_phenotypes.bar(ind,phenotypes['auc'], yerr=phenotypes['stdev'], ecolor='black',
#        zorder=1,color=two_colors[0],align='center')
plt.xticks(ind,phenotypes.index,rotation=45,ha='right')
axb__binary_phenotypes.spines['right'].set_visible(False)
axb__binary_phenotypes.spines['top'].set_visible(False)
axb__binary_phenotypes.tick_params(top='off',right='off',pad=2,labelsize=fontsize)
axb__binary_phenotypes.yaxis.set_ticks_position('left')
axb__binary_phenotypes.xaxis.set_ticks_position('bottom')
sns.despine(ax=axb__binary_phenotypes)

plt.ylabel('AUC')
plt.xlabel('')
plt.ylim(0.5,1)

plt.savefig(os.path.join(FIGURES_DIR, 'figureS1.pdf'), bbox_inches='tight', format='pdf')
plt.savefig(os.path.join(FIGURES_DIR, 'figureS1.png'), bbox_inches='tight', format='png')