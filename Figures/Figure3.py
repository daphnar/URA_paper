import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import os
import matplotlib.gridspec as gridspec

# import nature_guidline_utils
# from mne.stats import fdr_correction
#####TODO!!!!! One list of phenotypes!!!

rename = {'age':'Age',
          'bmi':'BMI',
          'hba1c':'HbA1C%',
          'bt__fasting_glucose':'Glucose',
          'bt__fasting_triglycerides':'Triglycerides',
          'bt__hdl_cholesterol': 'HDL cholesterol',
          'height':'Height',
          'bt__cholesterol': 'Total cholesterol',
          'bt__albumin':"Albumin",
          'bt__sgot': "SGOT",
          'bt__alkaline_phosphatase':"Alkaline phosphatase",
          'bt__tsh':'TSH',
          'bt__inr':'INR',
          'type_2_diabetes':'Type 2 diabetes',
          'currently_smokes':'Currently smokes',
          'ever_smoked':'Ever smoked',
          'bt__bilirubin':'Bilirubin'
          }


# 'bowel_movement_frequency':'Bowel movement',
#
# 'bt__sgpt': "SGPT",
# 'bt__protein':"Protein",

sns.set_style("ticks", {'axes.edgecolor': 'black'})
pd.set_option('display.width', 1000)
np.set_printoptions(precision=4, linewidth=200)
FIGURES_DIR = '/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/figures'
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
figsize=(183/10/2.54,247/10/2.54)
fig = plt.figure(figsize=figsize, dpi=300)


blue=(43./255,74./255,130./255)
green=(47./255,142./255,52./255)
pink=(191./255,64./255,140./255)
lightBlue=(42./255,107./255,126./255)
peach =  (205./255,112./255,106./255)
violet =  (106./255,111./255,205./255)

outer_grid = gridspec.GridSpec(1,3,width_ratios=[0.25,0.25,0.5])
# outer_grid.update(wspace=0.5)
ax__c = outer_grid[:,2]
# c_grid = gridspec.GridSpecFromSubplotSpec(1,2, ax__c,width_ratios=[0.3,0.7])

ax0 = plt.subplot(outer_grid[:,0])
ax0.set_axis_off()
ax1 = plt.subplot(outer_grid[:,1])
plt.text(-0.25, 1, 'a', ha='center', va='center', transform=ax1.transAxes, fontsize=16)

illustraitor_il=(34./255,181./255,115./255)
illustraitor_il_validation=(41./255,171./255,226./255)
illustraitor_us=(102./255,45./255,145./255)
ax2 = plt.subplot(outer_grid[:,2])
labels = ['Train+test1 (IL)','Test2 (US)']
colors_rgb = [(np.array(illustraitor_il)*9+np.array(illustraitor_il_validation))/10,
              illustraitor_us]

plt.text(-0.11, 1, 'b', ha='center', va='center', transform=ax2.transAxes, fontsize=16)

BEX_PATH_il = '/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/figures/LMM_results_corrected_binary-il.csv'
BEX_PATH_us = '/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/figures/LMM_results_corrected_binary-us.csv'

plt.sca(ax2)

h2CI_il = pd.read_csv(BEX_PATH_il).set_index('Phenotype')
h2CI_il.index=h2CI_il.index.map(lambda x: rename[x] if x in rename else '')#phenotypes.index.str.replace('bt__',"")
h2CI_il=h2CI_il.loc[h2CI_il.index[h2CI_il.index !='']]
h2CI_us = pd.read_csv(BEX_PATH_us).set_index('Phenotype')
h2CI_us.index=h2CI_us.index.map(lambda x: rename[x] if x in rename else '')#phenotypes.index.str.replace('bt__',"")
h2CI_us=h2CI_us.loc[h2CI_us.index[h2CI_us.index !='']]

phenotypes = h2CI_il.index.tolist()
phenotypes=phenotypes[::-1]
if 'Type 2 diabetes' not in phenotypes:
    phenotypes.append('Type 2 diabetes')
if 'Age' not in phenotypes:
    phenotypes.append('Age')
# phenotypes=phenotypes.tolist() ##############TODO!!!!!!!!!! remove this line!!!
# phenotypes.remove('currently_smokes') ##############TODO!!!!!!!!!! remove this line!!!
# # phenotypes.remove('impaired_glucose_tolerance_or_impaired_fasting_glucose') ##############TODO!!!!!!!!!! remove this line!!!
h2CI_il=h2CI_il.loc[phenotypes]
h2CI_il[["CI_high","CI_low"]]=h2CI_il[["CI_high","CI_low"]].fillna(0)
h2CI_us=h2CI_us.loc[phenotypes]
h2CI_us[["CI_high","CI_low"]]=h2CI_us[["CI_high","CI_low"]].fillna(0)
# # h2CI["pvalue"] = fdr_correction(h2CI["pvalue"].values)[1]

N = len(phenotypes)
width = 0.2
height = 0.05
ind = np.arange(N)





# h2CI.rename(index={'Hips cimcumference': 'Hip circumference'}, inplace=True)

ax2.barh(ind+width,(h2CI_il["CI_high"] - h2CI_il["CI_low"]).mul(100), height=0.05,edgecolor='none', left=h2CI_il['CI_low'] * 100,
                                                  zorder=1,
                                                  color=colors_rgb[0],tick_label=phenotypes)

ax2.barh(ind,(h2CI_us["CI_high"] - h2CI_us["CI_low"]).mul(100), height=0.05,edgecolor='none', left=h2CI_us['CI_low'] * 100,
                                                 zorder=1,
                                                 color=colors_rgb[-1])

ax2.legend(labels,loc=4,ncol=1, fontsize=fontsize - 2,frameon=False)#bbox_to_anchor=(-1.2, -0.095)


plt.scatter(h2CI_il["H2"].values * 100, ind+width+height/2, marker='o', color='black', zorder=2, linewidth=1, s=10)
plt.scatter(h2CI_us["H2"].values * 100, ind+height/2, marker='o', color='black', zorder=2, linewidth=1, s=10)


# ax2.barh(ind+width*2,(h2CI["CI_high"] - h2CI["CI_low"]).mul(100), height=height,edgecolor='none', left=h2CI['CI_low'] * 100,
#                                                   zorder=1,
#                                                   color=colors_rgb[0])
#
# ax2.barh(ind+width,(h2CI["CI_high"] - h2CI["CI_low"]).mul(100), height=0.05,edgecolor='none', left=h2CI['CI_low'] * 100,
#                                                   zorder=1,
#                                                   color=colors_rgb[-1],tick_label=phenotypes)
#
# ax2.barh(ind,(h2CI["CI_high"] - h2CI["CI_low"]).mul(100), height=0.05,edgecolor='none', left=h2CI['CI_low'] * 100,
#                                                  zorder=1,
#                                                  color=colors_rgb[-1])
#
# ax2.legend(labels,loc=4, bbox_to_anchor=(1.1, 0),ncol=1, fontsize=fontsize - 2,frameon=False)
#
#
# plt.scatter(h2CI["H2"].values * 100, ind+width*2+height/2, marker='o', color='black', zorder=2, linewidth=1, s=10)
# plt.scatter(h2CI["H2"].values * 100, ind+width+height/2, marker='o', color='black', zorder=2, linewidth=1, s=10)
# plt.scatter(h2CI["H2"].values * 100, ind+height/2, marker='o', color='black', zorder=2, linewidth=1, s=10)

plt.xlim(0, 65)
plt.xticks(np.arange(0, 65, 5))
ax = plt.gca()
box = ax.get_position()
ax.tick_params(top=False, right=False, pad=2, labelsize=fontsize - 2)
plt.xlabel('Microbiome-association index\n(% explained variance)', fontsize=fontsize - 2) #,labelpad=10
plt.ylabel('')
ax.set_yticklabels([])
# ax.tick_params(direction="outward", top='off', right='off', pad=2, labelsize=fontsize - 2)

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.xticks([0,10,20,30,40,50,60])



#######################ALPHA DIVERSITY

alpha_PATH_il = '/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/figures/Alpha-explained-variance-il_correct_binary.csv'
alpha_PATH_us = '/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/figures/Alpha-explained-variance-us_correct_binary.csv'
h2CI_il = pd.read_csv(alpha_PATH_il,index_col=0)
h2CI_il.index=h2CI_il.index.map(lambda x: rename[x] if x in rename else '')#phenotypes.index.str.replace('bt__',"")
h2CI_il=h2CI_il.loc[h2CI_il.index[h2CI_il.index !='']]
h2CI_us = pd.read_csv(alpha_PATH_us,index_col=0)
h2CI_us.index=h2CI_us.index.map(lambda x: rename[x] if x in rename else '')#phenotypes.index.str.replace('bt__',"")
h2CI_us=h2CI_us.loc[h2CI_us.index[h2CI_us.index !='']]

h2CI_il=h2CI_il.loc[phenotypes]
h2CI_us=h2CI_us.loc[phenotypes]
N = len(phenotypes)
width = 0.2
height = 0.05
ind = np.arange(N)




yticklabels = pd.Series(phenotypes)
yticklabels=yticklabels.apply(lambda x: x.replace('bt__',""))
yticklabels=yticklabels.apply(lambda x: x.replace('_'," "))
yticklabels=yticklabels.apply(lambda x: x.replace('frequency',""))
yticklabels=yticklabels.apply(lambda x: x.replace('freq',""))

ax1.barh(ind+width,(h2CI_il["CI_high"] - h2CI_il["CI_low"]).mul(100), height=0.05,edgecolor='none', left=h2CI_il['CI_low'] * 100,
                                                  zorder=1,
                                                  color=colors_rgb[0],tick_label=yticklabels)

ax1.barh(ind,(h2CI_us["CI_high"] - h2CI_us["CI_low"]).fillna(0).mul(100), height=0.05,edgecolor='none', left=h2CI_us['CI_low'].fillna(0).values * 100,
                                                 zorder=1,
                                                 color=colors_rgb[-1])

# ax1.legend(labels,loc=4, bbox_to_anchor=(1.1, 0.01),ncol=1, fontsize=fontsize - 2,frameon=False)

plt.sca(ax1)
plt.scatter(h2CI_il["H2"].values * 100, ind+width+height/2, marker='o', color='black', zorder=2, linewidth=1, s=10)
plt.scatter(h2CI_us["H2"].values * 100, ind+height/2, marker='o', color='black', zorder=2, linewidth=1, s=10)

plt.xlim(0, 15)
plt.xticks(np.arange(0, 16, 5))
box = ax1.get_position()
ax1.tick_params(top=False, right=False, pad=2, labelsize=fontsize - 2)
plt.xlabel('Alpha diversity\n(% explained variance)', fontsize=fontsize - 2) #,labelpad=10
plt.ylabel('')
# ax.tick_params(direction="outward", top='off', right='off', pad=2, labelsize=fontsize - 2)

ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
# ax.set_xticklabels(phenotypes)


plt.savefig(os.path.join(FIGURES_DIR, 'figure3_new.pdf'), format='pdf', bbox_inches='tight')
plt.savefig(os.path.join(FIGURES_DIR, 'figure3_new.png'), format='png', bbox_inches='tight')