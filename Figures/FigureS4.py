import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import os
import matplotlib as mpl
#from Unicorn.Figures import nature_guidline_utils

from dark import cmap_map

from URA_paper.Figures import nature_guidline_utils

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
                          nature_guidline_utils.full_page()*0.8*0.2), dpi=300)  # m2inch(165)

import matplotlib.gridspec as gridspec


green=(47./255,142./255,52./255)
violet =  (106./255,111./255,205./255)
two_colors = [violet,green]#colors_rgb[-2:]


axdrest_grid_bottom = gridspec.GridSpec(1,3,wspace=0.7,width_ratios=[1./3,1./3,1./3])

ax__age_saturation = plt.subplot(axdrest_grid_bottom[0,0])
ax__hba1c_saturation = plt.subplot(axdrest_grid_bottom[0,1])
ax__bmi_saturation = plt.subplot(axdrest_grid_bottom[0,2])


plt.sca(ax__age_saturation)
plt.text(-.53, 1.15, 'a', ha='center', va='center', transform=ax__age_saturation.transAxes, fontsize=16)
saturation_df=pd.read_csv(os.path.join(FIGURES_DIR,'Figures - age_saturation_US.csv'))
saturation_df.sort_values(by='cohort_size',inplace=True)
saturation_df.loc[:,['mean_pearson','mean_std']] = saturation_df.loc[:,['mean_pearson','mean_std']].mul(100)
ax__age_saturation.errorbar(saturation_df['cohort_size'],saturation_df['mean_pearson'],
                             yerr=saturation_df['mean_std'],color=two_colors[0])
saturation_df.loc[:,['mean_pearson_linear','mean_std_linear']] = saturation_df.loc[:,['mean_pearson_linear','mean_std_linear']].mul(100)
ax__age_saturation.errorbar(saturation_df['cohort_size'],saturation_df['mean_pearson_linear'],
                             yerr=saturation_df['mean_std_linear'],color=two_colors[1])
ax__age_saturation.tick_params(top='off',right='off',pad=2,labelsize=fontsize)
ax__age_saturation.yaxis.set_ticks_position('left')
ax__age_saturation.xaxis.set_ticks_position('bottom')
ax__age_saturation.spines['right'].set_visible(False)
ax__age_saturation.spines['top'].set_visible(False)
ax__age_saturation.set_xlim([0,20000])
ax__age_saturation.set_xticks([0,10000,20000])

plt.ylim(0,20)
ax__age_saturation.set_yticks([0,10,20])
plt.xlabel('Sample size')
plt.ylabel('age $R^{2}$ (%)')
plt.legend(['GBDT','Ridge'],bbox_to_anchor=(1.15, 0.0),frameon=False,loc=4)

plt.sca(ax__hba1c_saturation)
plt.text(-0.2, 1.15, 'b', ha='center', va='center', transform=ax__hba1c_saturation.transAxes, fontsize=16)
saturation_df=pd.read_csv(os.path.join(FIGURES_DIR,'Figures - hba1c_saturation_US.csv'))
saturation_df.sort_values(by='cohort_size',inplace=True)
saturation_df.loc[:,['mean_pearson','mean_std']] = saturation_df.loc[:,['mean_pearson','mean_std']].mul(100)
ax__hba1c_saturation.errorbar(saturation_df['cohort_size'],saturation_df['mean_pearson'],
                             yerr=saturation_df['mean_std'],color=two_colors[0])
saturation_df.loc[:,['mean_pearson_linear','mean_std_linear']] = saturation_df.loc[:,['mean_pearson_linear','mean_std_linear']].mul(100)
ax__hba1c_saturation.errorbar(saturation_df['cohort_size'],saturation_df['mean_pearson_linear'],
                             yerr=saturation_df['mean_std_linear'],color=two_colors[1])
ax__hba1c_saturation.tick_params(top='off',right='off',pad=2,labelsize=fontsize)
ax__hba1c_saturation.yaxis.set_ticks_position('left')
ax__hba1c_saturation.xaxis.set_ticks_position('bottom')
ax__hba1c_saturation.spines['right'].set_visible(False)
ax__hba1c_saturation.spines['top'].set_visible(False)
ax__hba1c_saturation.set_xlim([0,16000])
ax__hba1c_saturation.set_xticks([0,7500,15000])
plt.ylim(0,10)
ax__hba1c_saturation.set_yticks([0,5,10])#,20])
plt.xlabel('Sample size')
plt.ylabel('HbA1C% $R^{2}$ (%)')

plt.sca(ax__bmi_saturation)
plt.text(-.25, 1.15, 'c', ha='center', va='center', transform=ax__bmi_saturation.transAxes, fontsize=16)
saturation_df=pd.read_csv(os.path.join(FIGURES_DIR,'Figures - bmi_saturation_US.csv'))
saturation_df.sort_values(by='cohort_size',inplace=True)
saturation_df.loc[:,['mean_pearson','mean_std']] = saturation_df.loc[:,['mean_pearson','mean_std']].mul(100)
ax__bmi_saturation.errorbar(saturation_df['cohort_size'],saturation_df['mean_pearson'],
                             yerr=saturation_df['mean_std'],color=two_colors[0])
saturation_df.loc[:,['mean_pearson_linear','mean_std_linear']] = saturation_df.loc[:,['mean_pearson_linear','mean_std_linear']].mul(100)
ax__bmi_saturation.errorbar(saturation_df['cohort_size'],saturation_df['mean_pearson_linear'],
                             yerr=saturation_df['mean_std_linear'],color=two_colors[1])
ax__bmi_saturation.tick_params(top='off',right='off',pad=2,labelsize=fontsize)
ax__bmi_saturation.yaxis.set_ticks_position('left')
ax__bmi_saturation.xaxis.set_ticks_position('bottom')
ax__bmi_saturation.spines['right'].set_visible(False)
ax__bmi_saturation.spines['top'].set_visible(False)
ax__bmi_saturation.set_xlim([0,20000])
ax__bmi_saturation.set_xticks([0,10000,20000])
ax__bmi_saturation.set_yticks([0,5,10,15])
plt.xlabel('Sample size')
# plt.ylabel('')
plt.ylabel('BMI $R^{2}$ (%)')
plt.ylim(0,15)

plt.savefig(os.path.join(FIGURES_DIR, 'figureS_4.pdf'), bbox_inches='tight', format='pdf')
plt.savefig(os.path.join(FIGURES_DIR, 'figureS_4.png'), bbox_inches='tight', format='png')