import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import os
import matplotlib.gridspec as gridspec
from scipy.stats import pearsonr
# import nature_guidline_utils
# from mne.stats import fdr_correction

sns.set_style("ticks", {'axes.edgecolor': 'black'})
pd.set_option('display.width', 1000)
np.set_printoptions(precision=4, linewidth=200)
FIGURES_DIR = '/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/figures'
us_phenotypes=pd.read_csv(os.path.join(FIGURES_DIR,'all_phenotypes__us.csv'),index_col=0)
il_phenotypes=pd.read_csv(os.path.join(FIGURES_DIR,'all_phenotypes__il.csv'),index_col=0)
il_validation_phenotypes=pd.read_csv(os.path.join(FIGURES_DIR,'all_phenotypes__il_validation.csv'),index_col=0)

il_age_sex=pd.read_csv(os.path.join(FIGURES_DIR,'covariates__il.csv'),index_col=0)
il_validation_age_sex=pd.read_csv(os.path.join(FIGURES_DIR,'covariates__il_validation.csv'),index_col=0)
us_age_sex=pd.read_csv(os.path.join(FIGURES_DIR,'covariates__us.csv'),index_col=0)

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
figsize=(183/10/2.54,0.75*198/10/2.54) #80% of figure length

fig = plt.figure(figsize=figsize, dpi=300)


blue=(43./255,74./255,130./255)
green=(47./255,142./255,52./255)
pink=(191./255,64./255,140./255)
lightBlue=(42./255,107./255,126./255)
peach =  (205./255,112./255,106./255)
violet =  (106./255,111./255,205./255)
# colors_rgb = [green,lightBlue,violet]
illustraitor_il=(34./255,181./255,115./255)
illustraitor_il_validation=(41./255,171./255,226./255)
illustraitor_us=(102./255,45./255,145./255)
colors_rgb = [illustraitor_il,illustraitor_il_validation,illustraitor_us]

labels = ['Train-IL','Test1-IL','Test2-US']

# outer_grid = gridspec.GridSpec(1,2,width_ratios=[0.5,0.5])
# outer_grid.update(wspace=0.7)
# ax__abd = outer_grid[0,0]
# abd_outer_grid = gridspec.GridSpecFromSubplotSpec(2,1, ax__abd, hspace=0.2,height_ratios=[0.47,0.53])
abd_outer_grid = gridspec.GridSpec(3,3,width_ratios=[0.4,0.3,0.3],wspace=0.7,hspace=0.7)

# abd_grid = gridspec.GridSpecFromSubplotSpec(2,1, abd_outer_grid[1,0], hspace=0.35,height_ratios=[0.5,0.5])

ax_c = plt.subplot(abd_outer_grid[0,2])#abd_grid[0,0])
plt.sca(ax_c)
plt.text(-.25, 1.2, 'c', ha='center', va='center', transform=ax_c.transAxes, fontsize=16)
ax_c.spines['right'].set_visible(False)
ax_c.spines['top'].set_visible(False)
ax_c.tick_params(top=False,right=False,pad=2,labelsize=fontsize-2)
ax_c.yaxis.set_label_coords(-0.33,0.45)
# IL=np.random.normal(size=10000)
# US=np.random.normal(size=2000)
IL=il_phenotypes['hba1c'].dropna()
IL_validation=il_validation_phenotypes['hba1c'].dropna()
US=us_phenotypes['hba1c'].dropna()
# plt.hist(IL,bins=20,color=colors_rgb[0],density=True,histtype='step')#alpha=0.5,
# plt.hist(IL_validation,bins=20,color=colors_rgb[1],density=True,histtype='step')#alpha=0.5,
# plt.hist(US,bins=20,color=colors_rgb[2],density=True,histtype='step')#alpha=0.5,
sns.distplot(IL,hist=False, kde = True,color=colors_rgb[0],kde_kws = {'linewidth': 1})
sns.distplot(IL_validation,hist=False, kde = True,color=colors_rgb[1],kde_kws = {'linewidth': 1})
sns.distplot(US,hist=False, kde = True,color=colors_rgb[2],kde_kws = {'linewidth': 1})

plt.xlabel('HbA1C%',fontsize=fontsize - 2)
plt.ylabel('Frequency',fontsize=fontsize - 2)
#ax_c.yaxis.set_label_coords(-0.1,0.5)
plt.yticks([0,0.9])
plt.xticks([2,4,6,8,10])
ax_c.set_yticklabels([0,0.9])
ax_c.yaxis.set_label_coords(-0.33,0.45)
# plt.legend([labels[0],labels[-1]],loc='upper right',bbox_to_anchor=(1.1, 1),frameon=False,fontsize=fontsize - 2)
#plt.legend(labels,loc='upper right',bbox_to_anchor=(1.1, 1),frameon=False,fontsize=fontsize - 2)

# ax_b = plt.subplot(abd_grid[1,0])
ax_d = plt.subplot(abd_outer_grid[1,1])#abd_grid[0,0])
plt.sca(ax_d)
plt.text(-.25, 1.2, 'd', ha='center', va='center', transform=ax_d.transAxes, fontsize=16)
ax_d.spines['right'].set_visible(False)
ax_d.spines['top'].set_visible(False)
ax_d.tick_params(top=False,right=False,pad=2,labelsize=fontsize-2)
ax_d.yaxis.set_label_coords(-0.33,0.45)

IL=il_phenotypes['bmi'].dropna()
IL_validation=il_validation_phenotypes['bmi'].dropna()
US=us_phenotypes['bmi'].dropna()
plt.xticks([15,25,35,45,55])
plt.xlim([15,55])
# plt.hist(IL,bins=20,color=colors_rgb[0],density=True,histtype='step')
# plt.hist(IL_validation,bins=20,color=colors_rgb[1],density=True,histtype='step')#alpha=0.5,
# plt.hist(US,bins=20,color=colors_rgb[2],density=True,histtype='step')
sns.distplot(IL,hist=False, kde = True,color=colors_rgb[0],kde_kws = {'linewidth': 1})
sns.distplot(IL_validation,hist=False, kde = True,color=colors_rgb[1],kde_kws = {'linewidth': 1})
sns.distplot(US,hist=False, kde = True,color=colors_rgb[2],kde_kws = {'linewidth': 1})

plt.xlabel('BMI',fontsize=fontsize - 2)
plt.ylabel('Frequency',fontsize=fontsize - 2)
#ax_d.yaxis.set_label_coords(-0.1,0.5)
plt.yticks([0,0.09])
ax_d.set_yticklabels([0,0.09])
# plt.legend([labels[0],labels[-1]],loc='upper right',bbox_to_anchor=(1.1, 1),frameon=False,fontsize=fontsize - 2)
#plt.legend(labels,loc='upper right',bbox_to_anchor=(1.1, 1),frameon=False,fontsize=fontsize - 2)

ax_b = plt.subplot(abd_outer_grid[0,1])#abd_grid[0,0])
plt.sca(ax_b)
plt.text(-.25, 1.2, 'b', ha='center', va='center', transform=ax_b.transAxes, fontsize=16)
ax_b.spines['right'].set_visible(False)
ax_b.spines['top'].set_visible(False)
ax_b.tick_params(top=False,right=False,pad=2,labelsize=fontsize-2)

# IL=np.random.normal(size=10000,scale=2)
# US=np.random.normal(size=2000,scale=2)
IL=il_age_sex['age'].dropna()
IL_validation=il_validation_age_sex['age'].dropna()

# IL_male=il_age_sex['gender'].sum()
# IL_female=il_age_sex.shape[0]-IL_male
US=us_age_sex['age'].dropna()
# US_male=us_age_sex['gender'].sum()
# US_female=us_age_sex.shape[0]-US_male
# plt.hist(IL,bins=20,color=colors_rgb[0],density=True,histtype='step')
# plt.hist(IL_validation,bins=20,color=colors_rgb[1],density=True,histtype='step')
# plt.hist(US,bins=20,color=colors_rgb[2],density=True,histtype='step')
sns.distplot(IL,hist=False, kde = True,color=colors_rgb[0],label=labels[0],kde_kws = {'linewidth': 1})
sns.distplot(IL_validation,hist=False, kde = True,color=colors_rgb[1],label=labels[1],kde_kws = {'linewidth': 1})
sns.distplot(US,hist=False, kde = True,color=colors_rgb[2],label=labels[2],kde_kws = {'linewidth': 1})
plt.legend(labels,loc='upper right',bbox_to_anchor=(1.57, 1.39),
           frameon=False,fontsize=fontsize - 2,labelspacing=0.2,handletextpad=0.2)

plt.xlabel('Age',fontsize=fontsize - 2)
plt.xlim(15,100)
plt.xticks([18,40,60,80])
plt.ylabel('Frequency',fontsize=fontsize - 2)
ax_b.yaxis.set_label_coords(-0.33,0.45)
plt.yticks([0,0.03])
ax_b.set_yticklabels([0,0.03])
# plt.legend([labels[0]+'\nF: 7,338\nM: 4,637',
#             labels[-1]+'\nF: 1,240\nM: 940'],loc='upper right',bbox_to_anchor=(1.1, 1.1),frameon=False,fontsize=fontsize - 2)
# ax_b.legend()
ax_a = plt.subplot(abd_outer_grid[0,0])
plt.sca(ax_a)
plt.text(-.25, 1.2, 'a', ha='center', va='center', transform=ax_a.transAxes, fontsize=16)
ax_a.set_axis_off()

ax_e = plt.subplot(abd_outer_grid[1,2])
plt.sca(ax_e)
plt.text(-.25, 1.2, 'e', ha='center', va='center', transform=ax_e.transAxes, fontsize=16)
IL = pd.Series.from_csv(os.path.join(FIGURES_DIR,'Alpha-Shannon__il.csv'))
IL_validation = pd.Series.from_csv(os.path.join(FIGURES_DIR,'Alpha-Shannon__il_validation.csv'))
US = pd.Series.from_csv(os.path.join(FIGURES_DIR,'Alpha-Shannon__us.csv'))

sns.distplot(IL,hist=False, kde = True,color=colors_rgb[0],kde_kws = {'linewidth': 1})
sns.distplot(IL_validation,hist=False, kde = True,color=colors_rgb[1],kde_kws = {'linewidth': 1})
sns.distplot(US,hist=False, kde = True,color=colors_rgb[2],kde_kws = {'linewidth': 1})

plt.xlabel('Alpha diversity',fontsize=fontsize - 2)
plt.ylabel('Frequency',fontsize=fontsize - 2)
ax_e.spines['right'].set_visible(False)
ax_e.spines['top'].set_visible(False)
ax_e.tick_params(top=False,right=False,pad=2,labelsize=fontsize-2)
plt.yticks([0,0.5])
ax_e.set_yticklabels([0,0.5])
ax_e.yaxis.set_label_coords(-0.33,0.45)

cohort_df = pd.read_csv(os.path.join(FIGURES_DIR,'mean_species_abundance.csv'),index_col=0)
ax_f = plt.subplot(abd_outer_grid[2,1])
plt.sca(ax_f)
plt.text(-.25, 1.2, 'f', ha='center', va='center', transform=ax_f.transAxes, fontsize=16)
# plt.plot([-4,-1],[-4,-1],'k')
plt.plot([0.0001,0.1],[0.0001,0.1],'k')
plt.scatter(10**cohort_df['Train_IL'],10**cohort_df['Test_IL'],c=blue,alpha = 0.6,edgecolor='')
plt.xlabel('Train-IL species\nrelative abundance',fontsize=fontsize - 2)
plt.ylabel('Test1-IL species\nrelative abundance',fontsize=fontsize - 2)
r,p=pearsonr(cohort_df['Train_IL'],cohort_df['Test_IL'])
ax_f.annotate('R=%.2f'%r,(0.00015,0.05),fontsize=8)
# ax_f.yaxis.set_label_coords(-0.1,0.5)
plt.yscale('log')
plt.xscale('log')
plt.xlim([0.0001,0.1])
plt.ylim([0.0001,0.1])
plt.xticks([0.0001,0.001,0.01,0.1])
plt.yticks([0.0001,0.001,0.01,0.1])
plt.minorticks_off()
ax_f.spines['right'].set_visible(False)
ax_f.spines['top'].set_visible(False)
ax_f.tick_params(top=False,right=False,pad=2,labelsize=fontsize-2)
ax_f.yaxis.set_label_coords(-0.33,0.45)

ax_g = plt.subplot(abd_outer_grid[2,2])
plt.sca(ax_g)
plt.text(-.25, 1.2, 'g', ha='center', va='center', transform=ax_g.transAxes, fontsize=16)
plt.plot([0.0001,0.1],[0.0001,0.1],'k')
plt.scatter(10**cohort_df['Train_IL'],10**cohort_df['Test_US'],c=blue,alpha = 0.6,edgecolor='')
plt.xlabel('Train-IL species\nrelative abundance',fontsize=fontsize - 2)
plt.ylabel('Test2-US species\nrelative abundance',fontsize=fontsize - 2)
# ax_g.yaxis.set_label_coords(-0.1,0.5)
r,p=pearsonr(cohort_df['Train_IL'],cohort_df['Test_US'])
ax_g.annotate('R=%.2f'%r,(0.00015,0.05),fontsize=8)

plt.yscale('log')
plt.xscale('log')
plt.xlim([0.0001,0.1])
plt.ylim([0.0001,0.1])
plt.xticks([0.0001,0.001,0.01,0.1])
plt.yticks([0.0001,0.001,0.01,0.1])
plt.minorticks_off()

ax_g.spines['right'].set_visible(False)
ax_g.spines['top'].set_visible(False)
ax_g.tick_params(top=False,right=False,pad=2,labelsize=fontsize-2)
ax_g.yaxis.set_label_coords(-0.33,0.45)



plt.savefig(os.path.join(FIGURES_DIR, 'figure1_new.pdf'), format='pdf', bbox_inches='tight')
plt.savefig(os.path.join(FIGURES_DIR, 'figure1_new.png'), format='png', bbox_inches='tight')

# col_list = ["warm grey", "gunmetal", "dusky blue",
#             "cool blue", "deep teal",
#             "viridian", "turtle green", "twilight blue"][::-1]
#
# col_list = ["gunmetal", "dusky blue", "viridian"][::-1]
#
# colors_rgb = sns.xkcd_palette(col_list)


# ax_legend = plt.subplot(abd_grid[0,0])
# plt.sca(ax_legend)
# # plt.text(-.25, 1.2, 'c', ha='center', va='center', transform=ax_legend.transAxes, fontsize=16)
# ax = plt.gca()
# for i in [0,-1]:
#     plt.plot([], [], c=colors_rgb[i] , linewidth=5,
#                alpha=1, label=labels[i])#marker=markers[i],facecolors=colors_rgb[i]
# ax.set_axis_off()
# # ax.legend(bbox_to_anchor=(0.85, 1.4), ncol=3, fontsize=fontsize - 2)  # %loc='upper center',
# ax.legend(bbox_to_anchor=(0.85, 1.4), ncol=1, fontsize=fontsize - 2,frameon=False)
# plt.tick_params(labelbottom='off', labelleft='off')

# outer_grid = gridspec.GridSpec(2,1,height_ratios=[0.01,0.99])
# ax_legend = plt.subplot(c_grid[0,1])
#
# plt.sca(ax_legend)
# # plt.text(-.25, 1.2, 'c', ha='center', va='center', transform=ax_legend.transAxes, fontsize=16)
#
# ax = plt.gca()
#
#
# for i in range(len(labels)):
#     plt.plot([], [], c=colors_rgb[i] , linewidth=5,
#                alpha=1, label=labels[i])#marker=markers[i],facecolors=colors_rgb[i]
# ax.set_axis_off()
# # ax.legend(bbox_to_anchor=(0.85, 1.4), ncol=3, fontsize=fontsize - 2)  # %loc='upper center',
# ax.legend(bbox_to_anchor=(0.85, 1.4), ncol=1, fontsize=fontsize - 2,frameon=False)
# plt.tick_params(labelbottom='off', labelleft='off')

# c_grid = gridspec.GridSpecFromSubplotSpec(1,2, ax__c, width_ratios=[0.4,0.6])