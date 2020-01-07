import pandas as pd
from sklearn.metrics import r2_score
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import Ridge
from scipy.stats import pearsonr
import numpy
import matplotlib.pyplot as plt

col = 'bmi'
col = 'hba1c'

print(col)

if col == 'bmi':
    Sup ='/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/figures/Submission/MRS-BMI.csv'
    intercept = -44.79
    # IL r2 0.21359447301291157, pearson 0.522744
    # IL validation r2 0.15656481224420127, pearson 0.42498
    # US r2 0.03804002903446224, pearson 0.256346

else:
    Sup = '/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/figures/Submission/MRS-HbA1C.csv'
    intercept = 9.3
    # IL r2 0.19970355770042725, pearson 0.488118
    # IL validation r2 0.10032848377271375, pearson 0.348412
    # US r2 0.1146600578052892, pearson 0.367968


IL_phenotypes='/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/figures/all_phenotypes__il__il_validation.csv'
US_phenotypes='/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/figures/all_phenotypes__us.csv'
Bacteria_IL = '/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/figures/sgb__il.csv'
Bacteria_US = '/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/figures/sgb__us.csv'
Bacteria_ILv = '/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/figures/sgb__il_validation.csv'

IL_bac = pd.read_csv(Bacteria_IL,index_col=0)
IL_bac = IL_bac.apply(numpy.log)
IL_bacv = pd.read_csv(Bacteria_ILv,index_col=0)
IL_bacv = IL_bacv.apply(numpy.log)
US_bac = pd.read_csv(Bacteria_US,index_col=0)
US_bac = US_bac.apply(numpy.log)
IL = pd.read_csv(IL_phenotypes,index_col=0)[col]
IL.dropna(inplace=True)
US = pd.read_csv(US_phenotypes,index_col=0)[col]
US.dropna(inplace=True)
atlas = pd.read_csv(Sup).set_index('Species')
atlas = atlas['Ridge regression coefficient'].dropna()
top=25000
# atlas=atlas.sort_values('Ridge regression coefficient')

inds_IL = list(set(IL_bac.index).intersection(IL.index))
inds_ILv = list(set(IL_bacv.index).intersection(IL.index))
# coef = {}
# for ind in atlas.index[:top]:
#     reg = LinearRegression().fit(IL_bac.loc[inds_IL,[ind]], IL.loc[inds_IL])
#     coef[ind] = reg.coef_[0]

# reg = Ridge().fit(IL_bac.loc[inds_IL,atlas.index], IL.loc[inds_IL])
# for i, ind in enumerate(atlas.index[:top]):
#     coef[ind] = reg.coef_[i]

#atlas = atlas[atlas['Pearson correlation P value']<1]
IL_bac['sum_val'] = 0
IL_bacv['sum_val'] = 0
US_bac['sum_val'] = 0
for ind in atlas.index[:top]:
    IL_bac['sum_val'] += IL_bac[ind] * atlas[ind]
    IL_bacv['sum_val'] += IL_bacv[ind] * atlas[ind]
    US_bac['sum_val'] += US_bac[ind] * atlas[ind]

r_square=r2_score(IL.loc[inds_IL].values,intercept+IL_bac.loc[inds_IL]['sum_val'].values)
rp = pearsonr(IL.loc[inds_IL].values,IL_bac.loc[inds_IL]['sum_val'].values)
plt.figure()
plt.scatter(intercept+IL_bac.loc[inds_IL]['sum_val'].values,IL.loc[inds_IL].values)
plt.title('IL r2 %s, pearson %g' % (r_square, rp[0]))
plt.xlabel('predicted %s'%col)
plt.ylabel('real %s'%col)
print ('IL r2 %s, pearson %g' % (r_square, rp[0]))

r_square=r2_score(IL.loc[inds_ILv].values,intercept+IL_bacv.loc[inds_ILv]['sum_val'].values)
rp = pearsonr(IL.loc[inds_ILv].values,IL_bacv.loc[inds_ILv]['sum_val'].values)
plt.figure()
plt.scatter(intercept+IL_bacv.loc[inds_ILv]['sum_val'].values,IL.loc[inds_ILv].values)
plt.title('IL validation r2 %s, pearson %g' % (r_square, rp[0]))
plt.xlabel('predicted %s'%col)
plt.ylabel('real %s'%col)
print ('IL validation r2 %s, pearson %g' % (r_square, rp[0]))

inds_US = list(set(US_bac.index).intersection(US.index))
r_square=r2_score(US.loc[inds_US].values,intercept+US_bac.loc[inds_US]['sum_val'].values)
rp = pearsonr(US.loc[inds_US].values,US_bac.loc[inds_US]['sum_val'].values)
plt.figure()
plt.scatter(intercept+US_bac.loc[inds_US]['sum_val'].values,US.loc[inds_US].values)
plt.title('US r2 %s, pearson %g' % (r_square, rp[0]))
plt.xlabel('predicted %s'%col)
plt.ylabel('real %s'%col)
print ('US r2 %s, pearson %g' % (r_square, rp[0]))

plt.show()



