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
from scipy.stats import gaussian_kde
FIGURES_DIR = '/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/figures'

num_bacs = 3

plt.figure(figsize=(nature_guidline_utils.two_columns(),
                    nature_guidline_utils.full_page()*(3./4)), dpi=300)
params = {
    'axes.labelsize': 12,
    'font.size': 12,
    'legend.fontsize': 12,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'figure.dpi': 300,
    'axes.linewidth': 0.5,
}
plt.rcParams.update(params)

outer_grid = gridspec.GridSpec(2, 1, height_ratios=[ 1./3,2./3])
outer_grid.update(hspace=0.2)
abc_grid = gridspec.GridSpecFromSubplotSpec(1,6, outer_grid[1,:])
def_grid = gridspec.GridSpecFromSubplotSpec(1,3, outer_grid[0,:],wspace=0.4)

outer_grid.update(wspace=1)
color= "coolwarm"
colors=sns.color_palette('RdBu')
colors= [colors[0],colors[-1]]

def colorbar():
    m = cm.ScalarMappable(cmap=color)
    m.set_array([0, 1])
    cb = plt.colorbar(m, ticks=[0, 1])#, aspect=1000)
    cb.set_ticklabels(['Low', 'High'])
    cb.set_label('Feature value - bacterial abundance', size=10, labelpad=0)
    cb.ax.tick_params(labelsize=10, length=0)
    # cb.set_alpha(1)
    cb.outline.set_visible(False)
    # bbox = cb.ax.get_window_extent().transformed(plt.gcf().dpi_scale_trans.inverted())
    # cb.ax.set_aspect((bbox.height - 0.9) * 20)
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
a_and_b=False
if a_and_b:
    age_ax = plt.subplot(abc_grid[0, 1])
    plt.sca(age_ax)
    b_bmi_shap_values = pd.read_pickle(os.path.join(FIGURES_DIR, 'shap_values_age_linear.pkl'))
    b_bmi_mb = pd.read_pickle(os.path.join(FIGURES_DIR, 'df_for_shap_age_linear.pkl'))
    b_bmi_mb.columns = b_bmi_mb.columns.map(lambda x: keepname(x, False))
    summary_plot(b_bmi_shap_values, b_bmi_mb, show=False, auto_size_plot=True,
                 color_bar=False, color=color, phenotype='Age')
    for item in ([age_ax.title, age_ax.xaxis.label, age_ax.yaxis.label] +
                 age_ax.get_xticklabels() + age_ax.get_yticklabels()):
        item.set_fontsize(8)
    for item in ([age_ax.xaxis.label, age_ax.yaxis.label]):
        item.set_fontsize(10)
    age_ax.xaxis.set_label_coords(0,-0.07)


    b_hba1c_shap_values = pd.read_pickle(os.path.join(FIGURES_DIR,'shap_values_hba1c_linear.pkl'))
    b_hba1c_mb = pd.read_pickle(os.path.join(FIGURES_DIR,'df_for_shap_hba1c_linear.pkl'))
    b_hba1c_mb.columns=b_hba1c_mb.columns.map(lambda x: keepname(x,False))

    panelhba1c_ax = plt.subplot(abc_grid[0,0])
    plt.sca(panelhba1c_ax)
    plt.text(-0.15, 1, 'd', ha='center', va='center', transform=panelhba1c_ax.transAxes, fontsize=16)
    panelhba1c_ax.set_axis_off()

    hba1c_ax = plt.subplot(abc_grid[0,3])
    plt.sca(hba1c_ax)
    summary_plot(b_hba1c_shap_values, b_hba1c_mb,show=False,
                      auto_size_plot=True,color_bar=False,color=color,phenotype='HbA1C%')#plot_type='bar',
    for item in ([hba1c_ax.title, hba1c_ax.xaxis.label, hba1c_ax.yaxis.label] +
                  hba1c_ax.get_xticklabels() + hba1c_ax.get_yticklabels()):
        item.set_fontsize(8)
    for item in ([hba1c_ax.xaxis.label, hba1c_ax.yaxis.label]):
        item.set_fontsize(10)
    hba1c_ax.xaxis.set_label_coords(0,-0.07)
    hba1c_ax.set_xticks([-0.2,0,0.2,0.4])

    panelbmi_ax = plt.subplot(abc_grid[0,2])
    plt.sca(panelbmi_ax)
    plt.text(-0.15, 1, 'e', ha='center', va='center', transform=panelbmi_ax.transAxes, fontsize=16)
    panelbmi_ax.set_axis_off()
    panelbmi_ax = plt.subplot(abc_grid[0, 4])
    plt.sca(panelbmi_ax)
    plt.text(-0.15, 1, 'f', ha='center', va='center', transform=panelbmi_ax.transAxes, fontsize=16)
    panelbmi_ax.set_axis_off()

    bmi_ax = plt.subplot(abc_grid[0,5])
    plt.sca(bmi_ax)
    b_bmi_shap_values = pd.read_pickle(os.path.join(FIGURES_DIR,'shap_values_bmi_linear.pkl'))
    b_bmi_mb = pd.read_pickle(os.path.join(FIGURES_DIR,'df_for_shap_bmi_linear.pkl'))
    b_bmi_mb.columns=b_bmi_mb.columns.map(lambda x: keepname(x,False))
    summary_plot(b_bmi_shap_values, b_bmi_mb,show=False,auto_size_plot=True,
                      color_bar=False,color=color,phenotype='BMI')
    for item in ([bmi_ax.title, bmi_ax.xaxis.label, bmi_ax.yaxis.label] +
                  bmi_ax.get_xticklabels() + bmi_ax.get_yticklabels()):
        item.set_fontsize(8)
    for item in ([ bmi_ax.xaxis.label, bmi_ax.yaxis.label]):
            item.set_fontsize(10)
    colorbar()
    bmi_ax.set_xticks([-0.6,0,0.6])
    bmi_ax.xaxis.set_label_coords(0,-0.07)


def_data= pd.read_csv(os.path.join(FIGURES_DIR, 'Figures - il_vs_us_sgb_pvals.csv'))
d_ax = plt.subplot(def_grid[0,0])
plt.sca(d_ax)
d_ax.set_zorder(1)
plt.text(-0.15, 1.1, 'a', ha='center', va='center', transform=d_ax.transAxes, fontsize=16)
plt.title('Age',fontsize=10)
# age_pvals = def_data[def_data['pheno']=='age'][['pears_pval_us','pears_pval','species']].sort_values('pears_pval')
# age_pvals['species']=age_pvals['species'].apply((lambda x: keepname(x,False)))
# # plt.scatter(-np.log10(age_pvals['pears_pval']),-np.log10(age_pvals['pears_pval_us']),marker='.')
# x = -np.log10(age_pvals['pears_pval']).values
# y = -np.log10(age_pvals['pears_pval_us']).values
# xy = np.vstack([x,y])
# z = gaussian_kde(xy)(xy)
# idx = z.argsort()
# x, y, z = x[idx], y[idx], z[idx]
# res = d_ax.scatter(x,y, c=z, cmap=plt.cm.get_cmap('Blues_r'),s=2, edgecolor='')
# # cbar = plt.colorbar(res,ticks=[0.0002,  0.0025])
# # params = {'mathtext.default': 'regular' }
# # plt.rcParams.update(params)
# # cbar.ax.set_yticklabels(['$2x10^{-4}$', '$2.5x10^{-3}$'])
# plt.xlabel('P values (-log 10) IL test1')
# plt.ylabel('P values (-log 10) US test2')
# plt.xlim([0,70])
# plt.ylim([0,15])
#################################################################################################
# pheno_data = pd.read_csv(os.path.join(FIGURES_DIR, 'Figure5-%s-precent_agreement.csv'%'age'))
# pheno_data['Test-US-significant-percent']*=100
# plt.scatter(pheno_data.loc[pheno_data['Test-US-significant'].values == 0,'Train-IL-p-value'],
#             pheno_data.loc[pheno_data['Test-US-significant'].values == 0,'Test-US-significant-percent'],
#             c=colors[1],edgecolor='',alpha=0.5)#marker='.',s=20
# plt.scatter(pheno_data.loc[pheno_data['Test-US-significant'].values == 1,'Train-IL-p-value'],
#             pheno_data.loc[pheno_data['Test-US-significant'].values == 1,'Test-US-significant-percent'],
#             c=colors[0],edgecolor='',alpha=0.5)#marker='.',s=20
# # plt.text(-40,100,'Blautia (g), B.adolescentis (s)',fontsize=8)
# plt.xlabel('P-values train-IL (log10)')
# plt.ylabel('Cumulative significant species in test-US (%)')
# plt.xlim(pheno_data['Train-IL-p-value'][0]-5,0)
# plt.ylim(0,110)
#################################################################################################
def_data= pd.read_csv(os.path.join(FIGURES_DIR, 'Figures - il_vs_us_sgb_pvals.csv'))
def_data=def_data.set_index('species')
pheno_data=def_data[def_data['pheno']=='age']
pheno_data=pheno_data.sort_values('spear_pval')
bacteria = pheno_data.index.values
plt.scatter(pheno_data.loc[bacteria,'spear'],pheno_data.loc[bacteria,'spear_us'],# c=colors[1],
            c=np.log10(pheno_data['spear_pval']), cmap=plt.cm.get_cmap('RdBu'),
            norm=mpl.colors.Normalize(vmin=-20, vmax=0,clip=True),
            edgecolor='',alpha=0.5)

plt.xlabel('Train-IL Spearmann correlation')
plt.ylabel('Test2-US Spearmann correlation')
for i, txt in enumerate([keepname(sgb, False) for sgb in bacteria[:num_bacs]]):
    plt.scatter(pheno_data.loc[bacteria[i], 'spear'], pheno_data.loc[bacteria[i], 'spear_us'], c=colors[0])
    if txt=='S.parasanguinis (s)':
        plt.gca().annotate(txt,
                           (pheno_data.loc[bacteria, 'spear'].iloc[i] - 0.07,
                            pheno_data.loc[bacteria, 'spear_us'].iloc[i]+0.02),
                           fontsize=8)
    else:
        plt.gca().annotate(txt,
                       (pheno_data.loc[bacteria, 'spear'].iloc[i] + 0.01, pheno_data.loc[bacteria, 'spear_us'].iloc[i]),
                       fontsize=8)
plt.xlim(-0.20,0.20)
plt.ylim(-0.20,0.20)
plt.xticks([-0.20,0,0.20])
d_ax.set_xticklabels(["-0.2","0","0.2"])
plt.yticks([-0.2,0,0.2])
d_ax.set_yticklabels(["-0.2","0","0.2"])

e_ax = plt.subplot(def_grid[0,1])
plt.sca(e_ax)
plt.text(-0.15, 1.1, 'b', ha='center', va='center', transform=e_ax.transAxes, fontsize=16)
plt.title('HbA1C%',fontsize=10)
# hba1c_pvals = def_data[def_data['pheno']=='bmi'][['pears_pval_us','pears_pval','species']].sort_values('pears_pval')
# hba1c_pvals['species']=hba1c_pvals['species'].apply((lambda x: keepname(x,False)))
# # plt.scatter(-np.log10(hba1c_pvals['pears_pval']),-np.log10(hba1c_pvals['pears_pval_us']),marker='.')
# x = -np.log10(hba1c_pvals['pears_pval']).values
# y = -np.log10(hba1c_pvals['pears_pval_us']).values
# xy = np.vstack([x,y])
# z = gaussian_kde(xy)(xy)
# idx = z.argsort()
# x, y, z = x[idx], y[idx], z[idx]
# res = e_ax.scatter(x,y, c=z, cmap=plt.cm.get_cmap('Blues_r'),s=2, edgecolor='')
# # cbar = plt.colorbar(res,ticks=[0.0002,  0.0025])
# # params = {'mathtext.default': 'regular' }
# # plt.rcParams.update(params)
# # cbar.ax.set_yticklabels(['$2x10^{-4}$', '$2.5x10^{-3}$'])
# plt.xlabel('P values (-log 10) IL test1')
# plt.ylabel('P values (-log 10) US test2')
# plt.xlim([0,60])
# plt.ylim([0,20])
###########################################################################################################
# pheno_data = pd.read_csv(os.path.join(FIGURES_DIR, 'Figure5-%s-precent_agreement.csv'%'hba1c'))
# pheno_data['Test-US-significant-percent']*=100
# plt.scatter(pheno_data.loc[pheno_data['Test-US-significant'].values == 0,'Train-IL-p-value'],
#             pheno_data.loc[pheno_data['Test-US-significant'].values == 0,'Test-US-significant-percent'],
#             c=colors[1],edgecolor='',alpha=0.5)#marker='.',s=20
# plt.scatter(pheno_data.loc[pheno_data['Test-US-significant'].values == 1,'Train-IL-p-value'],
#             pheno_data.loc[pheno_data['Test-US-significant'].values == 1,'Test-US-significant-percent'],
#             c=colors[0],edgecolor='',alpha=0.5)#marker='.',s=20
# plt.xlabel('P-values train-IL (log10)')
# plt.xlim(pheno_data['Train-IL-p-value'][0]-5,0)
# plt.ylim(0,110)
# plt.xticks([-80,-60,-40,-20,0])
###########################################################################################################
def_data= pd.read_csv(os.path.join(FIGURES_DIR, 'Figures - il_vs_us_sgb_pvals.csv'))
def_data=def_data.set_index('species')
pheno_data=def_data[def_data['pheno']=='hba1c']
pheno_data=pheno_data.sort_values('spear_pval')
bacteria = pheno_data.index.values
plt.scatter(pheno_data.loc[bacteria,'spear'],pheno_data.loc[bacteria,'spear_us'],
            c=np.log10(pheno_data['spear_pval']), cmap=plt.cm.get_cmap('RdBu'),
            norm=mpl.colors.Normalize(vmin=-20, vmax=0,clip=True),
            edgecolor='',alpha=0.5)
plt.xlabel('Train-IL Spearmann correlation')
for i, txt in enumerate([keepname(sgb, False) for sgb in bacteria[:num_bacs]]):
    plt.scatter(pheno_data.loc[bacteria[i], 'spear'], pheno_data.loc[bacteria[i], 'spear_us'], c=colors[0])
    plt.gca().annotate(txt,
    (pheno_data.loc[bacteria, 'spear'].iloc[i] + 0.01, pheno_data.loc[bacteria, 'spear_us'].iloc[i]),
     fontsize=8)
plt.xlim(-0.2,0.2)
plt.ylim(-0.2,0.2)
plt.xticks([-0.2,0,0.25])
e_ax.set_xticklabels(["-0.2","0","0.25"])
plt.yticks([-0.2,0,0.25])
e_ax.set_yticklabels(["-0.2","0","0.25"])

f_ax = plt.subplot(def_grid[0,2])
plt.sca(f_ax)
plt.text(-0.15, 1.1, 'c', ha='center', va='center', transform=f_ax.transAxes, fontsize=16)
plt.title('BMI',fontsize=10)
#################################################################################################################3
# pheno_data = pd.read_csv(os.path.join(FIGURES_DIR, 'Figure5-%s-precent_agreement.csv'%'bmi'))
# pheno_data['Test-US-significant-percent']*=100
# plt.scatter(pheno_data.loc[pheno_data['Test-US-significant'].values == 0,'Train-IL-p-value'],
#             pheno_data.loc[pheno_data['Test-US-significant'].values == 0,'Test-US-significant-percent'],
#             c=colors[1],edgecolor='',alpha=0.5)#marker='.',s=20
# plt.scatter(pheno_data.loc[pheno_data['Test-US-significant'].values == 1,'Train-IL-p-value'],
#             pheno_data.loc[pheno_data['Test-US-significant'].values == 1,'Test-US-significant-percent'],
#             c=colors[0],edgecolor='',alpha=0.5)#marker='.',s=20
# plt.xlabel('P-values train-IL (log10)')
# plt.xlim(pheno_data['Train-IL-p-value'][0]-5,0)
# plt.ylim(0,110)
# plt.xticks([-70,-45,-25,0])
###################################################################################################################
def_data= pd.read_csv(os.path.join(FIGURES_DIR, 'Figures - il_vs_us_sgb_pvals.csv'))
def_data=def_data.set_index('species')
pheno_data=def_data[def_data['pheno']=='bmi']
pheno_data=pheno_data.sort_values('spear_pval')
bacteria = pheno_data.index.values
plt.scatter(pheno_data.loc[bacteria,'spear'],pheno_data.loc[bacteria,'spear_us'],
            c=np.log10(pheno_data['spear_pval']), cmap=plt.cm.get_cmap('RdBu'),
            norm=mpl.colors.Normalize(vmin=-20, vmax=0,clip=True),
            edgecolor='',alpha=0.5)

plt.xlabel('Train-IL Spearmann correlation')
for i, txt in enumerate([keepname(sgb, False) for sgb in bacteria[:num_bacs]]):
    plt.scatter(pheno_data.loc[bacteria[i], 'spear'], pheno_data.loc[bacteria[i], 'spear_us'], c=colors[0])
    plt.gca().annotate(txt, (pheno_data.loc[bacteria, 'spear'].iloc[i]+0.01, pheno_data.loc[bacteria, 'spear_us'].iloc[i]),
                       fontsize=8)
plt.xlim(-0.25,0.2)
plt.ylim(-0.25,0.2)
plt.xticks([-0.25,0,0.2])
f_ax.set_xticklabels(["-0.25","0","0.2"])
plt.yticks([-0.25,0,0.2])
f_ax.set_yticklabels(["-0.25","0","0.2"])
colorbar()
# bmi_pvals = def_data[def_data['pheno']=='bmi'][['pears_pval_us','pears_pval','species']].sort_values('pears_pval')
# bmi_pvals['species']=bmi_pvals['species'].apply((lambda x: keepname(x,False)))
# # plt.scatter(-np.log10(bmi_pvals['pears_pval']),-np.log10(bmi_pvals['pears_pval_us']),marker='.')
# x = -np.log10(bmi_pvals['pears_pval']).values
# y = -np.log10(bmi_pvals['pears_pval_us']).values
# xy = np.vstack([x,y])
# z = gaussian_kde(xy)(xy)
# idx = z.argsort()
# x, y, z = x[idx], y[idx], z[idx]
# res = f_ax.scatter(x,y, c=z, cmap=plt.cm.get_cmap('Blues_r'),s=2, edgecolor='')
# # cbar = plt.colorbar(res,ticks=[0.0002,  0.0025])
# # params = {'mathtext.default': 'regular' }
# # plt.rcParams.update(params)
# # cbar.ax.set_yticklabels(['$2x10^{-4}$', '$2.5x10^{-3}$'])
# plt.xlabel('P values (-log 10) IL test1')
# plt.ylabel('P values (-log 10) US test2')
# plt.xlim([0,60])
# plt.ylim([0,20])
# c_hba1c = pd.read_csv(os.path.join(FIGURES_DIR,'Figures - hba1c qq plot.csv'))
# c_hba1c['species']=c_hba1c['species'].apply(keepname)
# c_hba1c=c_hba1c.reset_index()
# plt.scatter(c_hba1c[c_hba1c.spear>0].index,-np.log10(c_hba1c[c_hba1c.spear>0].pval.values),marker='+',c=colors[0],s=30)
# plt.scatter(c_hba1c[c_hba1c.spear<0].index,-np.log10(c_hba1c[c_hba1c.spear<0].pval.values),marker='_',c=colors[1],s=30)
# plt.xlabel('Ordered species')
# plt.ylabel('P values (-log 10)')
# plt.gca().set_ylim(0,90)
# plt.gca().set_xlim(-70,max(c_hba1c[c_hba1c.spear<0].index))
# plt.text(70,80,'SGB 10068: E.coli (s)',fontsize=8)
# plt.text(70,39,'SGB 10064: E.marmotae (s)',fontsize=8)
# plt.text(70,26,'SGB 10115: K.pneumoniae (s)',fontsize=8)
# plt.text(70.15,20,'SGB 15286: Ruminococcacae (f)',fontsize=8)
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

# d_ax = plt.subplot(cd_grid[0,1])
# plt.sca(d_ax)
# plt.text(-0.15, 1.05, 'd', ha='center', va='center', transform=d_ax.transAxes, fontsize=16)
# d_bmi = pd.read_csv(os.path.join(FIGURES_DIR,'Figures - bmi qq plot.csv'))
# plt.scatter(d_bmi[d_bmi.spear>0].index,-np.log10(d_bmi[d_bmi.spear>0].pval.values),marker='+',c=colors[0],s=30)
# plt.scatter(d_bmi[d_bmi.spear<0].index,-np.log10(d_bmi[d_bmi.spear<0].pval.values),marker='_',c=colors[1],s=30)
# plt.xlabel('Ordered species')
# plt.ylabel('P values (-log 10)')
# plt.gca().set_ylim(0,75)
# plt.gca().set_xlim(-70,max(d_bmi[d_bmi.spear<0].index))
#
# # plotp.qqplot(d_bmi.pval.values, h1=d_ax,addlambda=False,xlim=[0,4],minpval=1e-90,
# #              group1=d_bmi.spear>0,group2=d_bmi.spear<0,markersize=7,color=colors)
# # d_ax.set_xlim(0,4)
# d_ax.spines['right'].set_visible(False)
# d_ax.spines['top'].set_visible(False)
# d_ax.tick_params(top=False, right=False)
# plt.text(70,70,'SGB 4964: Eubacteriaceae (f)',fontsize=8)
# plt.text(70,60,'SGB 4705: Clostridiceae (f)',fontsize=8)
# plt.text(70,51,' SGB 15369: F.sp_CAG_74(s)',fontsize=8)
# plt.text(70,46,'SGB 15249: F.sp_CAG_129 (s)',fontsize=8)
# for item in ([d_ax.xaxis.label, d_ax.yaxis.label]):
#     item.set_fontsize(10)

plt.savefig(os.path.join(FIGURES_DIR, 'figure5.pdf'), bbox_inches='tight', format='pdf')
plt.savefig(os.path.join(FIGURES_DIR, 'figure5.png'), bbox_inches='tight', format='png')
