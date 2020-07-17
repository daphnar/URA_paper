import pandas as pd
import os
# from Unicorn.Figures import nature_guidline_utils
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
from scipy.stats import ranksums,ks_2samp,pearsonr
from mne.stats import fdr_correction
import numpy as np
import matplotlib.gridspec as gridspec
from scipy.stats import gaussian_kde

sns.set_style("ticks", {'axes.edgecolor': 'black'})
# pd.set_option('display.width', 1000)
# np.set_printoptions(precision=4, linewidth=200)
illustraitor_il=(34./255,181./255,115./255)
illustraitor_il_validation=(41./255,171./255,226./255)
illustraitor_us=(102./255,45./255,145./255)
labels = ['Train+test1 (IL)','Test2 (US)']

params = {'axes.labelsize': 10,
          'axes.titlesize':8,#'text.fontsize': 8,
          'legend.fontsize': 10,
          'xtick.labelsize': 8,
          'ytick.labelsize': 8,
          'lines.linewidth' : 1,
          'lines.markersize' : 2}
blue=(43./255,74./255,130./255)
green=(47./255,142./255,52./255)
pink=(191./255,64./255,140./255)
lightBlue=(42./255,107./255,126./255)
peach =  (205./255,112./255,106./255)
violet =  (106./255,111./255,205./255)
colors_rgb = [blue,green]
colors_rgb = [(np.array(illustraitor_il)*9+np.array(illustraitor_il_validation))/10,
              illustraitor_us]

matplotlib.rcParams.update(params)
basepath='/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Cohort_Paper/Analyses/Alpha-Diversity'
# FIGURES_DIR = '/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/figures'
FIGURES_DIR = '/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Cohort_Paper/revision_Analyses/figures'


limits = {'age':[20,80],
          'bmi':[20,40],
          'hba1c':[4,8],
          'bowel_movement_frequency':[0,5],
          'bt__fasting_glucose':[70,170],
          'bt__fasting_triglycerides':[40,320],
          'bt__hdl_cholesterol':[30,90],
          'alpha':[3,7]}

rename = {'age':'Age (years)',
          'bmi':'BMI (Kg/m$^2$)',
          'hba1c':'HbA1C%',
          'bowel_movement_frequency':'Bowel movement frequency',
          'bt__fasting_glucose':'Blood glucose (mg/dl)',
          'bt__fasting_triglycerides':'Blood triglycerides (mg/dl)',
          'bt__hdl_cholesterol': 'HDL cholesterol (mg/dl)',
          'alpha': 'Alpha diversity'}

def params_for_subplots(ax,plot_asterix=False,remove_x=True):
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_title('')
    if remove_x:
        ax.set_xlabel('')
        ax.set_xticklabels([])
    ax.tick_params(top=False, right=False, pad=2)
    if plot_asterix:
        ymin, ymax = ax.get_ylim()
        x1, x2 = 0, 9
        y, h, col = ymin + (1.05 + 2 / 100.) * (ymax - ymin), 4.0 / 100 * (ymax - ymin), 'k'
        ll = plt.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1, c=col)
        plt.text((x1 + x2) * .5, y + h, "***", ha='center', va='bottom', color=col)
        ll[0].set_clip_on(False)

def plot_alpha_deciles_vs_pheontypes(ax_outer,pheno_df=None):
    if pheno_df is None:
        pheno_df= pd.read_csv(os.path.join(basepath, 'Phenotype-Alpha-Shannon__il__il_validation.csv'),index_col=0)
    pheno_df = pheno_df.sort_values('alpha')
    pheno_df['alpha-decile']=pd.qcut(pheno_df['alpha'],10,labels=[str(x) for x in range(1,11)])
    pheno_df=pheno_df[['alpha-decile', 'age', 'bmi', 'hba1c', 'bt__fasting_glucose',
              'bt__fasting_triglycerides', 'bt__hdl_cholesterol', 'alpha']]
    pheno_df['bt__fasting_triglycerides']=pheno_df['bt__fasting_triglycerides'].apply(lambda x: 10**x)
    ax_all = gridspec.GridSpecFromSubplotSpec(pheno_df.shape[1]-1, 1, ax_outer,hspace=0.55)
    ax_a = plt.subplot(ax_all[0, 0])
    plt.text(-.35, 1.1, 'a', ha='center', va='center', transform=ax_a.transAxes, fontsize=16)

    phenotype='alpha'
    ax_alpha = plt.subplot(ax_all[pheno_df.shape[1]-2, 0])
    ax_alpha = sns.boxplot(x=pheno_df['alpha-decile'].values.astype(int), y=pheno_df[phenotype].values,
                           color='white',fliersize=0,whis=[5, 95],width=0.5)
    ax_alpha.set_xlabel('Alpha diversity decile',labelpad=2)
    ax_alpha.set_ylabel('Alpha\ndiversity',labelpad=2)
    ax_alpha.set_ylim([1,7])
    ax_alpha.set_yticks([1,4, 7])
    ax_alpha.set_yticklabels([1, 4, 7])
    ax_alpha.spines['right'].set_visible(False)
    ax_alpha.spines['top'].set_visible(False)
    ax_alpha.set_title('')
    ax_alpha.tick_params(top=False, right=False, pad=2)
    pvals=[]
    stats = {}
    for i,phenotype in enumerate(pheno_df.columns):
        if phenotype == 'alpha-decile' or phenotype=='alpha':
            continue
        ax_p = plt.subplot(ax_all[i-1, 0])
        decile_df=pheno_df[['alpha-decile', phenotype]].pivot_table(values=phenotype,
               index=pheno_df[['alpha-decile', phenotype]].index,
               columns='alpha-decile', aggfunc='first')
        # print(decile_df[['1','10']].describe())
        all_stats = {}
        for j in range(10):
            all_stats[j] = [0]*10
            for k in range(j):
                res_rank = ranksums(decile_df[str(k+1)].dropna(), decile_df[str(j+1)].dropna())
                all_stats[j][k] = res_rank[1]
        pd.DataFrame(all_stats).to_csv(os.path.join(FIGURES_DIR,"fig2_stats_%s.csv"%phenotype))
        res_rank = ranksums(decile_df['1'].dropna(),decile_df['10'].dropna())
        res_ks = ks_2samp(decile_df['1'].dropna(),decile_df['10'].dropna())
        stats[phenotype] = [res_rank[1],res_ks[1]]
        ax_p = sns.boxplot(x=pheno_df['alpha-decile'], y=pheno_df[phenotype],
                           color='white',fliersize=0,whis=[5, 95],width=0.6)
        ax_p.set_ylabel(rename[phenotype].replace(' ','\n'),labelpad=2)
        ax_p.set_yticks([limits[phenotype][0],(limits[phenotype][0]+limits[phenotype][1])/2,limits[phenotype][1]])
        ax_p.set_ylim(limits[phenotype])
        params_for_subplots(ax_p,plot_asterix=True)
        plt.subplots_adjust(left=0.3)
    stats_df=pd.DataFrame(stats,index=['RankSum_Pvalue','KS_Pvalue']).T
    stats_df['RankSum_Qvalue']=fdr_correction(stats_df['RankSum_Pvalue'].values)[1]
    stats_df['KS_Qvalue']=fdr_correction(stats_df['KS_Pvalue'].values)[1]
    stats_df.to_csv(os.path.join(FIGURES_DIR,"Figure2_stats.csv"))



def calcMeanWindow(vector,s=1000,jump=100):
    current=0
    res = []
    while current<len(vector):
        if ((current + s) > len(vector)):
            break
        res.append(np.mean(vector[current:current+s]))
        current = current + jump
    res.append(np.mean(vector[-s:]))
    return res


def plot_pheotype_trend(ax_outer):
    s = 1000
    ax_all = gridspec.GridSpecFromSubplotSpec(7, 1, ax_outer,hspace=0.55)# hspace=0.55
    ax_b = plt.subplot(ax_all[0, 0])
    for j,dataset in enumerate(['il__il_validation','us']):
        pheno_df= pd.read_csv(os.path.join(basepath, 'Phenotype-Alpha-Shannon__%s.csv'%dataset),index_col=0)
        pheno_df['alpha-decile']=pd.qcut(pheno_df['alpha'],10,labels=[str(x) for x in range(1,11)])
        pheno_df=pheno_df[['alpha-decile','age', 'bmi', 'hba1c', 'bt__fasting_glucose',
                  'bt__fasting_triglycerides', 'bt__hdl_cholesterol', 'alpha']]
        # ax_all = gridspec.GridSpecFromSubplotSpec(pheno_df.shape[1] - 1, 1, ax_outer,hspace=0.55)
        ax_b = plt.subplot(ax_all[0, 0])
        plt.text(-.4, 1.1, 'b', ha='center', va='center', transform=ax_b.transAxes, fontsize=16)
        phenotype_alpha_correlations=[]
        pvalues=[]
        for i,phenotype in enumerate(pheno_df.columns):
            if phenotype=='hba1c':
                pass
            pheno_sorted_df = pheno_df.sort_values(phenotype).dropna()
            if s*19./10 > len(pheno_sorted_df):
                used_s = 10* int(len(pheno_sorted_df)/19.)
            else:
                used_s = s
            if phenotype == 'alpha-decile' or phenotype=='alpha':
                continue
            ax_p = plt.subplot(ax_all[i - 1, 0])
            window_phenotype = calcMeanWindow(pheno_sorted_df[phenotype], used_s, int(used_s/10))
            window_alpha = calcMeanWindow(pheno_sorted_df['alpha'], used_s, int(used_s/10))
            r,p=pearsonr(pheno_sorted_df[phenotype],pheno_sorted_df['alpha'])
            phenotype_alpha_correlations.append(r)
            pvalues.append(p)
            outlier_high_alpha=pheno_sorted_df['alpha'].mean()+2*pheno_sorted_df['alpha'].std()
            # outlier_low_alpha = pheno_sorted_df['alpha'].mean() - 2 * pheno_sorted_df['alpha'].std()
            # non_outlier_pheno = pheno_sorted_df.loc[(pheno_sorted_df['alpha']<=outlier_high_alpha) & \
            #                                         (pheno_sorted_df['alpha'] >=outlier_low_alpha)]
            # ax_p.scatter(non_outlier_pheno[phenotype],non_outlier_pheno['alpha'],c='gray',marker='.',alpha=0.5)
            x = pheno_sorted_df[phenotype].reset_index()[phenotype].values
            y = pheno_sorted_df['alpha'].reset_index()['alpha'].values
            y = y+np.random.random(y.shape[0])/100
            xy = np.vstack([x, y])
            z = gaussian_kde(xy)(xy)
            print('yey %s'%phenotype)
            idx = z.argsort()
            x, y, z = x[idx], y[idx], z[idx]
            res = ax_p.scatter(x, y, c=z, cmap=plt.cm.get_cmap('gray'), s=1, edgecolor='',alpha=0.3)
            #cbar = plt.colorbar(res, shrink=1.15, pad=0.025)#, ticks=[0.0002, 0.0025]

            ax_p.plot(window_phenotype,window_alpha,c=colors_rgb[j])
            ax_p.autoscale(enable=True, axis='x', tight=True)

            ax_p.scatter(window_phenotype,window_alpha,c=colors_rgb[j], edgecolor='',s=3)
            ax_p.set_ylabel('Alpha\ndiversity',labelpad=2)
            # ax_p.set_yticks([4.9,5.5] )
            # ax_p.set_yticklabels([4.9,5.5] )
            # ax_p.set_ylim([4.9,5.5] )
            #ax_p.set_yticklabels([2,4.5,7])
            ax_p.set_yticks([1,4,7])
            ax_p.set_ylim([1, 7])
            ax_p.set_yticklabels([1,4,7])
            ax_p.set_xlabel(rename[phenotype],labelpad=2)
            params_for_subplots(ax_p,remove_x=False)
            # if phenotype == 'bt__hdl_cholesterol' or phenotype =='age':
            #     proportion_x = ax_p.get_xlim()[0] + (ax_p.get_xlim()[1] - ax_p.get_xlim()[0]) * 0.02
            # else:
            #     proportion_xroportion_x = ax_p.get_xlim()[0] + (ax_p.get_xlim()[1] - ax_p.get_xlim()[0]) * 0.7
            proportion_x = ax_p.get_xlim()[0] + (ax_p.get_xlim()[1] - ax_p.get_xlim()[0])
            proportion_y = ax_p.get_ylim()[0] + (ax_p.get_ylim()[1] - ax_p.get_ylim()[0]) * 0.6

            if dataset=='il__il_validation':
                ax_p.annotate('R=%.2f\n'
                              r'P<$10^{%d}$' % (r, np.log10(p)-1),
                              (proportion_x, proportion_y), fontsize=8)
    ax_p = plt.subplot(ax_all[i - 1, 0])
    for i in range(2):
        plt.plot([], [],color=colors_rgb[i],linewidth=2)
    ax_p.set_axis_off()
    ax_p.legend(labels, loc=10, ncol=1, frameon=False)

def m2inch(value,unit='mm'):
    if unit=='cm':
        return value/2.54
    if unit=='mm':
        return value/10/2.54

def two_columns():
    return m2inch(183)

def full_page():
    return m2inch(247)

def plot_two_parts():
    outer_grid = gridspec.GridSpec(1, 2, width_ratios=[0.5, 0.5])
    outer_grid.update(wspace=0.5)
    plt.figure(figsize=(two_columns(),
                               full_page()), dpi=300)
    ax__a = outer_grid[0, 0]
    ax__b = outer_grid[0, 1]
    plot_alpha_deciles_vs_pheontypes(ax__a)
    plot_pheotype_trend(ax__b)
    plt.savefig(os.path.join(FIGURES_DIR, 'figure2_revision.pdf'), bbox_inches='tight', format='pdf')
    plt.savefig(os.path.join(FIGURES_DIR, 'figure2_revision.png'), bbox_inches='tight', format='png')
if __name__ == '__main__':
    plot_two_parts()
    # plot_pheotype_trend()
    # plot_alpha_deciles_vs_pheontypes()