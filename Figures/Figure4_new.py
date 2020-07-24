import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import os
import matplotlib as mpl
#from Unicorn.Figures import nature_guidline_utils
from mne.stats import fdr_correction
from dark import cmap_map
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
from scipy.stats import gaussian_kde
from sklearn.metrics import r2_score

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
                          nature_guidline_utils.full_page()), dpi=300)  # m2inch(165)

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
# colors_rgb_age = [nothing]+colors_rgb[-3:]+[nothing]
# colors_rgb_bmi_hba1c = colors_rgb[:4]+[colors_rgb[-1]]
colors_rgb_age = two_colors#colors_rgb
colors_rgb_bmi_hba1c = two_colors#colors_rgb
two_colors_9=two_colors*9
pattern_mb = ['/']
pattern_age_sex = ['\\']
pattern_both = ['x']
# illustraitor_il=(34./255,181./255,115./255)
# illustraitor_il_validation=(41./255,171./255,226./255)
# illustraitor_us=(102./255,45./255,145./255)
# colors_rgb_age = [nothing]+[illustraitor_il,illustraitor_il_validation,illustraitor_us][::-1]+[nothing]
# patterns_age = ["/","/","/"]
# colors_rgb_bmi_hba1c = [illustraitor_il,illustraitor_il,illustraitor_il,illustraitor_il_validation,illustraitor_us][::-1]
# patterns_rgb_bmi_hba1c = ["'/","\\","x","x","x"]
# patterns = [ "/" , "\\" , "|" , "-" , "+" , "x", "o", "O", ".", "*" ]

# colors_rgb_age = [nothing,green,lightBlue,peach,nothing]
# colors_rgb_bmi_hba1c = [blue,red,pink,violet,peach]
# colors_rgb = [blue,red,pink,violet,green,lightBlue,peach]


outer_grid = gridspec.GridSpec(2,1,height_ratios=[0.20,0.80])
outer_grid.update(hspace=0.45)
ax__ab = outer_grid[0,0]
ax__rest = outer_grid[1,0]
ab_grid = gridspec.GridSpecFromSubplotSpec(1,2, ax__ab, wspace=0.2, width_ratios=[17./21,4./21])
axa__quantitive_phenotypes = plt.subplot(ab_grid[0,0])
axa__quantitive_phenotypes.set_zorder(1)
axb__binary_phenotypes = plt.subplot(ab_grid[0,1])
axb__binary_phenotypes.set_zorder(1)
#rest_grid = gridspec.GridSpecFromSubplotSpec(3,2, ax__rest,hspace=0.6, wspace=0.3,height_ratios=[0.4,0.3,0.3])
rest_grid = gridspec.GridSpecFromSubplotSpec(3,1, ax__rest,hspace=0.5, height_ratios=[0.2,0.4,0.2])

axc_grid = gridspec.GridSpecFromSubplotSpec(1,3, rest_grid[0,0], wspace=0.7,width_ratios=[1./3,1./3,1./3])

axdrest_grid_top=gridspec.GridSpecFromSubplotSpec(1,2,rest_grid[1,0],wspace=0.3,width_ratios=[0.5,0.5])
axc__hba1c_bmi_bar = plt.subplot(axdrest_grid_top[0,0])
axc__hba1c_bmi_legend = plt.subplot(axdrest_grid_top[0,1])
axdrest_grid_bottom=gridspec.GridSpecFromSubplotSpec(1,3,rest_grid[2,0],wspace=0.7,width_ratios=[1./3,1./3,1./3])


ax__age_scatter = plt.subplot(axc_grid[0,0])
ax__age_saturation = plt.subplot(axdrest_grid_bottom[0,0])
ax__hba1c_scatter = plt.subplot(axc_grid[0,1])
ax__hba1c_saturation = plt.subplot(axdrest_grid_bottom[0,1])
ax__bmi_scatter = plt.subplot(axc_grid[0,2])
ax__bmi_saturation = plt.subplot(axdrest_grid_bottom[0,2])


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
          'bt__inr':'INR ',
          'T2D':'Type 2 diabetes',
          'currently_smokes':'Currently smokes',
          'ever_smoked':'Ever smoked',
          'type_2_diabetes': 'Type 2 diabetes',
          'gender': 'Gender'
          }

plt.sca(axa__quantitive_phenotypes)
plt.text(-.15, 1.1, 'a', ha='center', va='center', transform=axa__quantitive_phenotypes.transAxes, fontsize=16)
correlations_df = pd.read_csv(os.path.join(FIGURES_DIR,'Figures - prediction_other_phenotypes.csv')).set_index('target')
# phenotypes = correlations_df[['pearson','stdev']].dropna().sort_values(by='pearson',ascending=False).mul(100)
phenotypes = correlations_df[['pearson','stdev','pearson_male','stdev_male','pearson_female','stdev_female']].dropna().sort_values(by='pearson',ascending=False).mul(100)
phenotypes.columns = ['All','std_all','Male','std_male','Female','std_female']
phenotypes.index=phenotypes.index.map(lambda x: rename[x] if x in rename else '')#phenotypes.index.str.replace('bt__',"")
phenotypes=phenotypes.loc[phenotypes.index[phenotypes.index !='']]
ind = np.array(range(len(phenotypes)))


# axa__quantitive_phenotypes.bar(ind,phenotypes['pearson'],yerr=phenotypes['stdev'], ecolor='black',#edgecolor='none',
#        zorder=1,color=colors_rgb[0],align='center')
estimation = phenotypes[['All','Male','Female']]
yerr = phenotypes[['std_all','std_male','std_female']].T.values
estimation.plot.bar(yerr=yerr,ax = axa__quantitive_phenotypes,
                              zorder=1,#edgecolor='none',
                              color=three_colors,legend=True,
                              width=0.8)
axa__quantitive_phenotypes.legend(frameon=False)#loc=10,
plt.ylim(-10,36.5)
axa__quantitive_phenotypes.spines['left'].set_bounds(0, 36.5)
plt.yticks([0,10,20,30])
axa__quantitive_phenotypes.spines['bottom'].set_position('zero')
axa__quantitive_phenotypes.set_xlim(ind.min()-0.5,ind.max()+0.5)
plt.xticks(ind,phenotypes.index,rotation=45,ha='right')
# axa__quantitive_phenotypes.spines['right'].set_visible(False)
# axa__quantitive_phenotypes.spines['top'].set_visible(False)
axa__quantitive_phenotypes.tick_params(top='off',right='off',pad=0,labelsize=fontsize)
axa__quantitive_phenotypes.yaxis.set_ticks_position('left')
axa__quantitive_phenotypes.xaxis.set_ticks_position('bottom')
sns.despine(ax=axa__quantitive_phenotypes)
plt.ylabel('$R^{2}$ (%)')
plt.xlabel('')


plt.sca(axb__binary_phenotypes)
correlations_df = pd.read_csv(os.path.join(FIGURES_DIR,'Figures - classification.csv')).set_index('target')
phenotypes = correlations_df[['auc','stdev','pearson_male','stdev_male','pearson_female','stdev_female']].sort_values(by='auc',ascending=False)
phenotypes.columns = ['All','std_all','Male','std_male','Female','std_female']
phenotypes.index=phenotypes.index.map(lambda x: rename[x] if x in rename else '')#phenotypes.index.str.replace('bt__',"")
plt.text(-.6, 1.1, 'b', ha='center', va='center', transform=axb__binary_phenotypes.transAxes, fontsize=16)
ind = range(len(phenotypes))
estimation = phenotypes[['All','Male','Female']]
yerr = phenotypes[['std_all','std_male','std_female']].T.values
estimation.plot.bar(yerr=yerr,ax = axb__binary_phenotypes,
                              zorder=1,#edgecolor='none',
                              color=three_colors,legend=False,
                              width=0.8)

# axb__binary_phenotypes.bar(ind,phenotypes['auc'], yerr=phenotypes['stdev'], ecolor='black',
#        zorder=1,color=two_colors[0],align='center')
plt.xticks(ind,phenotypes.index,rotation=45,ha='right')
axb__binary_phenotypes.spines['right'].set_visible(False)
axb__binary_phenotypes.spines['top'].set_visible(False)
axb__binary_phenotypes.tick_params(top='off',right='off',pad=2,labelsize=fontsize)
axb__binary_phenotypes.yaxis.set_ticks_position('left')
axb__binary_phenotypes.xaxis.set_ticks_position('bottom')
plt.ylabel('AUC')
plt.xlabel('')
plt.ylim(0.5,1)

plt.sca(axc__hba1c_bmi_bar)

df_age_hba1c_bmi = pd.read_csv(os.path.join(FIGURES_DIR,'Figures - new_4c.csv'))
df_age_hba1c_bmi['phenotype']=df_age_hba1c_bmi['pheno'].apply(lambda x: rename[x])
df_age_hba1c_bmi[['r2','stdev']]*=100
df_age_hba1c_bmi[['mean_pearson','sixe_fixed_stdev']]=df_age_hba1c_bmi[['r2','stdev']]
labels=[]
for i in range(df_age_hba1c_bmi.shape[0]):
    initial=[]
    if df_age_hba1c_bmi.iloc[i]['cohort']=='il':
        initial='CV(Train-IL): '
    elif df_age_hba1c_bmi.iloc[i]['cohort']=='validation_il' or\
            df_age_hba1c_bmi.iloc[i]['cohort']=='IL_validation':
        initial = 'Test1-IL: '
    else:
        initial = 'Test2-US: '
    features = df_age_hba1c_bmi.iloc[i]['features']
    features=features.replace('sex','gender')
    features = features.replace('_', '+')
    features = features.replace('sgbs', 'microbiome')
    features = features.replace('microbiome+age+gender','age+gender+microbiome')
    labels.append(initial+features)
df_age_hba1c_bmi['labels']=labels
df_hba1c_bmi=df_age_hba1c_bmi.groupby(['phenotype','labels']).aggregate('first').unstack()

df_hba1c_bmi=df_hba1c_bmi.loc[['BMI','HbA1C%','Age'],]
plt.text(-.28, 1.0, 'f', ha='center', va='center', transform=axc__hba1c_bmi_bar.transAxes, fontsize=16)
# phenotypes = df_hba1c_bmi[['mean_pearson','std_pearson']].sort_values(by='mean_pearson').mul(100)
phenotypes = df_hba1c_bmi[['mean_pearson','sixe_fixed_stdev']]
phenotypes=phenotypes
ind = range(len(phenotypes))

order=[u'Test1-IL: age+gender',u'Test2-US: age+gender',u'Test1-IL: microbiome',u'Test2-US: microbiome',
 u'Test1-IL: age+gender+microbiome',u'Test2-US: age+gender+microbiome'][::-1]

# orderAge=order
# orderBMI_HbA1C=order
#
# onlyage_prediction=phenotypes['mean_pearson'].loc[:,orderAge]
# onlyage_prediction.loc[['BMI','HbA1C%'],:]=np.nan
# onlyage_err=phenotypes['sixe_fixed_stdev'].loc[:,orderAge]
# onlyage_err.loc[['BMI','HbA1C%'],:]=np.nan
#
# onlyage_prediction.plot.barh(xerr=onlyage_err,ax = axc__hba1c_bmi_bar,
#                               zorder=1,edgecolor='none',
#                              color=colors_rgb_age,
#                               legend=False,
#                               width=0.8)
#
# noage_prediction=phenotypes['mean_pearson'].loc[:,orderBMI_HbA1C]
# noage_prediction.loc['Age',:]=np.nan
# noage_err=phenotypes['sixe_fixed_stdev'].loc[:,orderBMI_HbA1C]
# noage_err.loc['Age',:]=np.nan
#
# noage_prediction.plot.barh(xerr=noage_err,ax = axc__hba1c_bmi_bar,
#                               zorder=1,edgecolor='none',
#                               color=colors_rgb_bmi_hba1c,
#                               legend=False,
#                               width=0.8)


phenotypes['mean_pearson'].loc[:,order].plot.barh(xerr=phenotypes['sixe_fixed_stdev'].loc[:,order],ax = axc__hba1c_bmi_bar,
                              zorder=1,#,edgecolor=two_colors,
                              color=colors_rgb,legend=False,
                              width=0.8)

# bars = axc__hba1c_bmi_bar.patches
# patterns =( "...","..." , "OOO","OOO"  , ".O.",".O.")
# hatches = [patterns[i%len(patterns)] for i in range(len(bars))]
# for bar, hatch in zip(bars, hatches):
#     bar.set_hatch(hatch,scale=0.5)
#
# bars = axc__hba1c_bmi_bar.patches
# hatches = pattern_age_sex*2+pattern_mb*2+pattern_both*2
# for bar, hatch in zip(bars, hatches*3):
#     bar.set_hatch(hatch)

plt.ylabel("")
# axc__hba1c_bmi_bar.set_yticklabels('HbA1C%','BMI'])
plt.xlabel('$R^{2}$ (%)')
plt.xlim(-2,40)
plt.ylim(-0.5,2.5)
axc__hba1c_bmi_bar.spines['left'].set_position('zero')

axc__hba1c_bmi_bar.tick_params(top='off',right='off',pad=2,labelsize=fontsize)
axc__hba1c_bmi_bar.yaxis.set_ticks_position('left')
axc__hba1c_bmi_bar.xaxis.set_ticks_position('bottom')
axc__hba1c_bmi_bar.spines['right'].set_visible(False)
axc__hba1c_bmi_bar.spines['top'].set_visible(False)


plt.sca(axc__hba1c_bmi_legend)
labels=df_hba1c_bmi.columns.levels[1].values
labels=order#labels[order]
for i in range(len(labels))[::-1]:
    plt.plot([], [], c=colors_rgb[i] , linewidth=5,
               alpha=1, label=labels[i])
axc__hba1c_bmi_legend.set_axis_off()
axc__hba1c_bmi_legend.legend(bbox_to_anchor=(-0.35, 0.5),loc=6,ncol=1, fontsize=fontsize,frameon=False)#loc=10,
axc__hba1c_bmi_legend.tick_params(labelbottom='off', labelleft='off')

plt.sca(ax__age_scatter)
plt.text(-.63, 1.1, 'c', ha='center', va='center', transform=ax__age_scatter.transAxes, fontsize=16)
age_scatter_df = pd.read_csv(os.path.join(FIGURES_DIR,'Figures - age_xy_scatter.csv'))
x = age_scatter_df['x\n']
y = age_scatter_df['y\n']
xy = np.vstack([x,y])
z = gaussian_kde(xy)(xy)
idx = z.argsort()
x, y, z = x[idx], y[idx], z[idx]
res = ax__age_scatter.scatter(x,y, c=z, cmap=plt.cm.get_cmap('Blues_r'),s=2, edgecolor='')
r_square=r2_score(x,y)
ax__age_scatter.annotate('$R^{2}$=%.2f'%r_square,(12,90),fontsize=8)
cbar = plt.colorbar(res,ticks=[0.0002,  0.0025],shrink=1.15,pad=0.025)
cbar.ax.tick_params(axis='both', which='major', pad=1)
cbar.outline.set_visible(False)
params = {'mathtext.default': 'regular' }
plt.rcParams.update(params)
cbar.ax.set_yticklabels(['$2x10^{-4}$', '$2.5x10^{-3}$'])
# plt.tick_params(labelleft='off')
# ax__age_scatter.scatter(age_scatter_df['x'],age_scatter_df['y'],color=colors_rgb[0],s=2)
ax__age_scatter.tick_params(top='off',right='off',pad=2,labelsize=fontsize)
ax__age_scatter.yaxis.set_ticks_position('left')
ax__age_scatter.xaxis.set_ticks_position('bottom')
ax__age_scatter.spines['right'].set_visible(False)
ax__age_scatter.spines['top'].set_visible(False)
# ax__age_scatter.set_yticks(([30,50,70]))
ax__age_scatter.set_yticks(([10,30,50,70,90]))
ax__age_scatter.set_xticks(([10,30,50,70,90]))
ax__age_scatter.plot([10,30,50,70,90],[10,30,50,70,90],'k',linewidth=0.5)
plt.ylabel('Predicted age')
plt.xlabel('Actual age')

plt.sca(ax__age_saturation)
plt.text(-.53, 1.15, 'g', ha='center', va='center', transform=ax__age_saturation.transAxes, fontsize=16)
saturation_df=pd.read_csv(os.path.join(FIGURES_DIR,'Figures - age_saturation2.csv'))
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

plt.ylim(0,35)
ax__age_saturation.set_yticks([0,10,20,30])
plt.xlabel('Sample size')
plt.ylabel('age $R^{2}$ (%)')
plt.legend(['GBDT','Ridge'],bbox_to_anchor=(1.15, 0.0),frameon=False,loc=4)

plt.sca(ax__hba1c_scatter)
plt.text(-.2, 1.1, 'd', ha='center', va='center', transform=ax__hba1c_scatter.transAxes, fontsize=16)
hba1c_scatter_df = pd.read_csv(os.path.join(FIGURES_DIR,'Figures - hba1c_xy_scatter.csv'))
# ax__hba1c_scatter.scatter(hba1c_scatter_df['x'],hba1c_scatter_df['y'],color=colors_rgb[0],s=2)
x = hba1c_scatter_df['x']
y = hba1c_scatter_df['y']
xy = np.vstack([x,y])
z = gaussian_kde(xy)(xy)
idx = z.argsort()
x, y, z = x[idx], y[idx], z[idx]
res = ax__hba1c_scatter.scatter(x,y, c=z, cmap=plt.cm.get_cmap('Blues_r'),s=2, edgecolor='')
r_square=r2_score(x,y)
ax__hba1c_scatter.annotate('$R^{2}$=%.2f'%r_square,(2.2,10),fontsize=8)
cbar = plt.colorbar(res,ticks=[0.1, 0.5, 0'c'.9],shrink=1.15,pad=0.025)
cbar.ax.tick_params(axis='both', which='major', pad=1)
cbar.outline.set_visible(False)

ax__hba1c_scatter.tick_params(top='off',right='off',pad=2,labelsize=fontsize)
ax__hba1c_scatter.yaxis.set_ticks_position('left')
ax__hba1c_scatter.xaxis.set_ticks_position('bottom')
ax__hba1c_scatter.spines['right'].set_visible(False)
ax__hba1c_scatter.spines['top'].set_visible(False)
plt.ylabel('Predicted\nHbA1C%')
plt.xlabel('Actual HbA1C%')
ax__hba1c_scatter.set_xticks([2,4,6,8,10])
ax__hba1c_scatter.set_yticks([2,4,6,8,10])
ax__hba1c_scatter.plot([2,4,6,8,10],[2,4,6,8,10],'k',linewidth=0.5)
#ax__hba1c_scatter.set_yticks([4.5,5.5,6.5,7.5])
plt.sca(ax__bmi_scatter)
plt.text(-.3, 1.1, 'e', ha='center', va='center', transform=ax__bmi_scatter.transAxes, fontsize=16)
bmi_scatter_df = pd.read_csv(os.path.join(FIGURES_DIR,'Figures - bmi_xy_scatter.csv'))
# ax__bmi_scatter.scatter(bmi_scatter_df['x'],bmi_scatter_df['y'],color=colors_rgb[0],s=2)
x = bmi_scatter_df['x']
y = bmi_scatter_df['y']
xy = np.vstack([x,y])
z = gaussian_kde(xy)(xy)
idx = z.argsort()
x, y, z = x[idx], y[idx], z[idx]
res = ax__bmi_scatter.scatter(x,y, c=z, cmap=plt.cm.get_cmap('Blues_r'),s=2, edgecolor='')
r_square=r2_score(x,y)
ax__bmi_scatter.annotate('$R^{2}$=%.2f'%r_square,(10.5,50),fontsize=8)
cbar = plt.colorbar(res,ticks=[0.002,0.01,0.018],shrink=1.15,pad=0.025)
cbar.ax.tick_params(axis='both', which='major', pad=1)
cbar.outline.set_visible(False)
ax__bmi_scatter.tick_params(top='off',right='off',pad=2,labelsize=fontsize)
ax__bmi_scatter.yaxis.set_ticks_position('left')
ax__bmi_scatter.xaxis.set_ticks_position('bottom')
ax__bmi_scatter.spines['right'].set_visible(False)
ax__bmi_scatter.spines['top'].set_visible(False)
ax__bmi_scatter.set_xticks([10,20,30,40,50])
#ax__bmi_scatter.set_yticks([20,25,30,35])
ax__bmi_scatter.set_yticks([10,20,30,40,50])
ax__bmi_scatter.plot([10,20,30,40,50],[10,20,30,40,50],'k',linewidth='0.5')
plt.ylabel('Predicted BMI')
plt.xlabel('Actual BMI')

plt.sca(ax__hba1c_saturation)
plt.text(-0.2, 1.15, 'h', ha='center', va='center', transform=ax__hba1c_saturation.transAxes, fontsize=16)
saturation_df=pd.read_csv(os.path.join(FIGURES_DIR,'Figures - hba1c_saturation2.csv'))
saturation_df.sort_values(by='cohort_size',inplace=True)
saturation_df.loc[:,['mean_pearson','mean_std']] = saturation_df.loc[:,['mean_pearson','mean_std']].mul(100)
ax__hba1c_saturation.errorbar(saturation_d'c'f['cohort_size'],saturation_df['mean_pearson'],
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
plt.ylim(0,30)
ax__hba1c_saturation.set_yticks([0,10,20,30])#,20])
plt.xlabel('Sample size')
plt.ylabel('HbA1C% $R^{2}$ (%)')

plt.sca(ax__bmi_saturation)
plt.text(-.25, 1.15, 'i', ha='center', va='center', transform=ax__bmi_saturation.transAxes, fontsize=16)
saturation_df=pd.read_csv(os.path.join(FIGURES_DIR,'Figures - bmi_saturation2.csv'))
saturation_df.sort_values(by='cohort_size',inplace=True)
saturation_df.loc[:,['mean_pearson','mean_std']] = saturation_df.loc[:,['mean_pearson','mean_std']].mul(100)
ax__bmi_saturation.errorbar(saturation_df['cohort_size'],saturation_df['mean_pearson'],
                             yerr=saturation_df['mean_std'],color=two_colors[0])
saturation_df.loc[:,['mean_pearson_linear','mean_std_linear']] = saturation_df.loc[:,['mean_pearson_linear','mean_std_linear']].mul(100)
ax__bmi_saturation.errorbar(saturation_df['c''c''cohort_size'],saturation_df['mean_pearson_linear'],
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
plt.ylim(0,17)
# axf__hba1c_saturation.set_yticks([0,5,10,15,20])
plt.savefig(os.path.join(FIGURES_DIR, 'figure4_colors.pdf'), bbox_inches='tight', format='pdf')
plt.savefig(os.path.join(FIGURES_DIR, 'figure4_colors.png'), bbox_inches='tight', format='png')