import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import os
#from Unicorn.Figures import nature_guidline_utils
from scipy.stats import gaussian_kde
from sklearn.metrics import r2_score
from sklearn.metrics._ranking import roc_curve, roc_auc_score
import matplotlib.gridspec as gridspec

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

phenotypes = ['age','bmi','hba1c',"bt__fasting_glucose","bt__fasting_triglycerides",'bowel_movement_frequency',
              "bt__hdl_cholesterol","bt__cholesterol","bt__protein","bt__albumin","bt__sgot","bt__alkaline_phosphatase",
              "bt__tsh","bt__inr",'gender','currently_smokes','ever_smoked','type_2_diabetes']
letters='abcdefghijklmnopqr'
for predictor in ['xgb','Ridge']:
    fig = plt.figure(figsize=(nature_guidline_utils.two_columns(),
                              nature_guidline_utils.full_page()), dpi=300)  # m2inch(165)

    outer_grid = gridspec.GridSpec(6, 3, hspace=0.9, wspace=0.7)#, width_ratios=[1. / 3, 1. / 3, 1. / 3])
    for i,pheno in enumerate(phenotypes):
        print(pheno)
        ax = plt.subplot(outer_grid[i])
        plt.sca(ax)
        plt.text(-.5, 1.2, letters[i], ha='center', va='center', transform=ax.transAxes, fontsize=16)
        scatter_df = pd.read_csv(os.path.join(FIGURES_DIR,'%s_MB_pred_%s.csv'%(predictor,pheno)))
        x = scatter_df['real_val']
        if pheno in ['gender','currently_smokes','ever_smoked','type_2_diabetes']:
            y=scatter_df['pred_prob']
            fpr, tpr, _ = roc_curve(x, y)
            auc = roc_auc_score(x, y)
            # plot the roc curve for the model
            plt.plot(fpr,tpr, linewidth=0.5 , marker='.',markersize=1)
            # axis labels
            plt.xlabel('False positive rate')
            plt.ylabel('True positive rate')
        else:
            y = scatter_df['pred_val']
            xy = np.vstack([x,y])
            z = gaussian_kde(xy)(xy)
            idx = z.argsort()
            x, y, z = x[idx], y[idx], z[idx]
            res = ax.scatter(x,y, c=z, cmap=plt.cm.get_cmap('Blues_r'),s=2, edgecolor='')
            r_square=r2_score(x,y)
            #cbar = plt.colorbar(res,ticks=[0.0002,  0.0025],shrink=1.15,pad=0.025)
            #cbar.ax.tick_params(axis='both', which='major', pad=1)
            #cbar.outline.set_visible(False)
            params = {'mathtext.default': 'regular' }
            plt.rcParams.update(params)
            #cbar.ax.set_yticklabels(['$2x10^{-4}$', '$2.5x10^{-3}$'])
            plt.ylabel('Predicted value' )
            plt.xlabel('Actual value' )

        plt.title(rename[pheno],fontsize=10)
        ax.plot(plt.xlim(), plt.xlim(), 'k', linewidth=0.5)
        if pheno in ['gender','currently_smokes','ever_smoked','type_2_diabetes']:
            ax.annotate('AUC=%.2f'%auc,(plt.xlim()[0]+(plt.xlim()[1]-plt.xlim()[0])/20,
                                        plt.ylim()[0]+(plt.ylim()[1]-plt.ylim()[0])*0.9),fontsize=8)
        else:
            ax.annotate('$R^{2}$=%.2f'%r_square,(plt.xlim()[0]+(plt.xlim()[1]-plt.xlim()[0])/20,
                                                 plt.ylim()[0]+(plt.ylim()[1]-plt.ylim()[0])*0.9),fontsize=8)
        ax.tick_params(top='off',right='off',pad=2,labelsize=fontsize)
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
    #plt.savefig(os.path.join(FIGURES_DIR, 'figureS2-3%s.pdf'%predictor), bbox_inches='tight', format='pdf')
    print('Saving %s'%predictor)
    plt.savefig(os.path.join(FIGURES_DIR, 'figureS2-3%s.png'%predictor), bbox_inches='tight', format='png')
    plt.close()