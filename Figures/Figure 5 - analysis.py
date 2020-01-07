import os
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
from mne.stats import fdr_correction
import numpy as np
FIGURES_DIR = '/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/figures'
threshold =0.05
phenotypes = ['hba1c','bmi','age']
lw = 2
def_data= pd.read_csv(os.path.join(FIGURES_DIR, 'Figures - il_vs_us_sgb_pvals.csv'))
def_data=def_data.set_index('species')
for phenotype in phenotypes:
    pheno_data=def_data[def_data['pheno']==phenotype]
    pheno_data=pheno_data.sort_values('spear_pval')
    # significant_il = def_data['spear_pval']<threshold/len(def_data.index)
    significant_il = fdr_correction(pheno_data['spear_pval'])[1]<threshold
    pheno_data=pheno_data.loc[significant_il]
    significant_us= fdr_correction(pheno_data['spear_pval_us'])[1]<threshold
    print len(pheno_data.loc[significant_us].index)
    print pheno_data.loc[significant_us].index
    cumsum=np.cumsum(significant_us)
    precent_significant=cumsum=1.0*cumsum/range(1,len(significant_us)+1)
    plt.scatter(np.log10(pheno_data['spear_pval']),precent_significant)
    plt.xlabel('P-values train-IL (log10)')
    plt.ylabel('Significant percent test-US')
    plt.title(phenotype)
    plt.xlim(np.log10(pheno_data['spear_pval'])[0]-1,0)
    plt.ylim(0,1.1)
    pd.DataFrame({'species':pheno_data.index,
        'Train-IL-p-value':np.log10(pheno_data['spear_pval'].values),
        'Test-US-significant':significant_us.astype(int),
        'Test-US-significant-percent':precent_significant}).set_index('species')\
        .to_csv(os.path.join(FIGURES_DIR,'Figure5-%s-precent_agreement.csv'%phenotype))
    plt.savefig(os.path.join(FIGURES_DIR,'Figure5-%s-precent_agreement.png'%phenotype), format='png')
    plt.close()
# fpr, tpr, _ = roc_curve(significant_us.values.astype(int), [1]*len(significant_us))
# roc_auc = auc(fpr, tpr)
# plt.plot(fpr[2], tpr[2], color='darkorange',
#          lw=lw, label='ROC curve (area = %0.2f)' % roc_auc[2])
# plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
# plt.xlim([0.0, 1.0])
# plt.ylim([0.0, 1.05])
# plt.xlabel('False Positive Rate')
# plt.ylabel('True Positive Rate')
# plt.title('Receiver operating characteristic example')
# plt.legend(loc="lower right")
# plt.show()

