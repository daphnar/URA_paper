from scipy.spatial.distance import squareform,pdist
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from scipy.stats import mannwhitneyu
def BrayCurtis(X_orig, is_log_abundance=True, zero_min_value=True):
    if is_log_abundance:
        X = 10 ** X_orig
    else:
        X = X_orig.copy()
    if zero_min_value: X[X_orig == np.min(X_orig)] = 0

    D = squareform(pdist(X, metric='braycurtis'))
    return D

def compare_us_il(train_il_path,test_il_path,test_us_path,save_bc_path=None,save_fig_path=None):
    train_il=pd.read_csv(train_il_path).set_index('client_id')
    test_il = pd.read_csv(test_il_path).set_index('client_id')
    test_us = pd.read_csv(test_us_path).set_index('client_id')
    print('Loaded data')
    all_il = pd.concat([train_il,test_il])
    il_ids = all_il.index
    us_ils = test_us.index
    all_data = pd.concat([all_il,test_us])
    print('Concated data')
    all_data_sgb = all_data.loc[:,all_data.columns.str.startswith('k__')]
    bc_all = BrayCurtis(all_data_sgb,False)
    bd_all = pd.DataFrame(data=bc_all,index=all_data.index,columns=all_data.index)
    print('Calculated BC')
    # if save_path is not None:
    #     bd_all.to_csv(save_path)
    #     print('Saved BC')
    ils=bc_all[il_ids,il_ids].values.flatten()
    uss=bc_all[us_ils,us_ils].values.flatten()
    comb=bc_all[il_ids,us_ils].values.flatten()
    plt.boxplot([ils,uss,comb],color='white',fliersize=0,whis=[5, 95],width=0.5)
    plt.xlabel(['IL-IL','US-US','IL-US'])
    plt.ylabel('Bray Curtis diversity')
    plt.savefig(save_fig_path+'.png')
    plt.savefig(save_fig_path+'.pdf',format='pdf')
    print("should not necessarily be significantly different ",mannwhitneyu(ils,uss))
    print("should hopefully be significant ",mannwhitneyu(ils,comb))
    print("should hopefully be significant ",mannwhitneyu(uss,comb))

if __name__=='__main__':
    data_path='/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Cohort_Paper/DataForRevision'
    train_il_path = os.path.join(data_path, 'ura_q_il.csv')
    test_il_path = os.path.join(data_path, 'ura_q_il_validation.csv')
    test_us_path = os.path.join(data_path, 'ura_q_us.csv')
    bc_output_path='/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Cohort_Paper/revision_Analyses/bc_all.csv'
    fig_output_path='/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Cohort_Paper/revision_Analyses/bc_fig'

    compare_us_il(train_il_path, test_il_path, test_us_path,bc_output_path,fig_output_path)


