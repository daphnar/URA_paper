import numpy as np
import os
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt
import pandas as pd
import scipy.linalg as la
import sys

def BrayCurtis(X_orig, is_log_abundance=True, zero_min_value=True):
    if is_log_abundance: X = 10**X_orig
    else: X = X_orig.copy()
    if zero_min_value:
        X[X_orig==np.min(X_orig)]=0
    D = squareform(pdist(X, metric='braycurtis'))
    return D


def e_matrix(distance_matrix):
    """Compute E matrix from a distance matrix.
    Squares and divides by -2 the input elementwise. Eq. 9.20 in
    Legendre & Legendre 1998."""
    return distance_matrix * distance_matrix / -2.0


def f_matrix(E_matrix):
    """Compute F matrix from E matrix.
    Centring step: for each element, the mean of the corresponding
    row and column are substracted, and the mean of the whole
    matrix is added. Eq. 9.21 in Legendre & Legendre 1998."""
    print type(E_matrix)
    row_means = E_matrix.mean(axis=1)  # , keepdims=True)
    col_means = E_matrix.mean(axis=0)  # , keepdims=True)
    matrix_mean = E_matrix.mean()
    return E_matrix - row_means - col_means + matrix_mean


def PCoA(D, suppress_warning=False, return_explained_variance=False):
    E_matrix = e_matrix(D)
    F_matrix = f_matrix(E_matrix)
    eigvals, eigvecs = la.eigh(F_matrix)
    negative_close_to_zero = np.isclose(eigvals, 0)
    eigvals[negative_close_to_zero] = 0
    if (np.any(eigvals < 0) and not suppress_warning):
        print(
            "The result contains negative eigenvalues."
            " Please compare their magnitude with the magnitude of some"
            " of the largest positive eigenvalues. If the negative ones"
            " are smaller, it's probably safe to ignore them, but if they"
            " are large in magnitude, the results won't be useful. See the"
            " Notes section for more details. The smallest eigenvalue is"
            " {0} and the largest is {1}.".format(eigvals.min(),
                                                  eigvals.max()),
        )

    idxs_descending = eigvals.argsort()[::-1]
    eigvals = eigvals[idxs_descending]
    eigvecs = eigvecs[:, idxs_descending]
    num_positive = (eigvals >= 0).sum()
    eigvecs[:, num_positive:] = np.zeros(eigvecs[:, num_positive:].shape)
    eigvals[num_positive:] = np.zeros(eigvals[num_positive:].shape)
    coordinates = eigvecs * np.sqrt(eigvals)

    if return_explained_variance:
        return coordinates, eigvals / eigvals.sum()
    return coordinates


basepath = '/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Analyses/SGB_composition'
cohort_data = {'Train_IL':'sgb__il.csv', 'Test_IL':'sgb__il_validation.csv',
               'Test_US':'sgb__us.csv'}

cohort_df = {}
plt.figure()
illustraitor_il=(34./255,181./255,115./255)
illustraitor_il_validation=(41./255,171./255,226./255)
illustraitor_us=(102./255,45./255,145./255)
color = {'Train_IL':illustraitor_il,'Test_IL':illustraitor_il_validation,'Test_US':illustraitor_us}
order = ['Train_IL','Test_IL','Test_US']
for cohort in cohort_data:
    cohort_df[cohort] = pd.read_csv(os.path.join(basepath,cohort_data[cohort]),index_col=0)
    cohort_df[cohort][cohort_df[cohort]==0]=0.0001
    cohort_df[cohort]=cohort_df[cohort].apply(np.log10)
pd.DataFrame({'Train_IL':cohort_df['Train_IL'].mean(),
              'Test_IL':cohort_df['Test_IL'].mean(),
              'Test_US':cohort_df['Test_US'].mean()}).to_csv(os.path.join(basepath,'mean_species_abundance.csv'))
sys.exit()
all_cohort = pd.concat([cohort_df[x] for x in order])
if False:
    BC = BrayCurtis(all_cohort,False,False)
    allPCs, all_pcoa_explained_var = PCoA(BC, return_explained_variance=True)
    pd.DataFrame(allPCs,index=all_cohort.index)[[0,1]].to_csv(os.path.join(basepath,'PCoA_log_BC.csv'))
    pd.Series(all_pcoa_explained_var)[[0,1]].to_csv(os.path.join(basepath,'PCoA_explained_var_log_BC.csv'))
else:
    allPCs = pd.read_csv(os.path.join(basepath,'PCoA_log_BC.csv'),index_col=0).T.values
    all_pcoa_explained_var = pd.Series.from_csv(os.path.join(basepath,'PCoA_explained_var_log_BC.csv')).values


pos = 0
plt.figure()
for cohort in order:
    plt.scatter(allPCs[0][pos:pos+len(cohort_df[cohort])],allPCs[1][pos:pos+len(cohort_df[cohort])],
                c=color[cohort],alpha=0.5,edgecolor='')
    pos += len(cohort_df[cohort])
plt.legend(order)

plt.savefig(os.path.join(basepath,'PCoA_log.png'))
plt.close('all')
plt.figure()
plt.scatter(cohort_df['Train_IL'].mean(),cohort_df['Test_IL'].mean())
plt.xlabel('Train_IL mean species log abundance')
plt.ylabel('Test_IL mean species log abundance')
plt.plot([-4,-1],[-4,-1],'r')
plt.figure()
plt.scatter(cohort_df['Train_IL'].mean(),cohort_df['Test_US'].mean())
plt.xlabel('Train_IL mean species log abundance')
plt.ylabel('Test_US mean species log abundance')
plt.plot([-4,-1],[-4,-1],'r')
plt.show()



