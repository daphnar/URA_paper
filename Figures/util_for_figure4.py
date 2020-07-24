import pandas as pd

def merger(pheno):
    obj1 = pd.read_csv('saturation_xgboostMB_%s_vs_ILval.csv'%pheno)
    obj2 = pd.read_csv('saturation_RidgeMB_%s_vs_ILval.csv'%pheno)
    merged = obj1.merge(obj2,right_index=True,left_index=True)
    merged.columns = [0,'cohort_size','mean_pearson','mean_std',1,2,'mean_pearson_linear','mean_std_linear']
    merged = merged[['cohort_size','mean_pearson','mean_std','mean_pearson_linear','mean_std_linear']]
    merged.to_csv('Figures - %s_saturation2.csv'%pheno)

merger('hba1c')
merger('bmi')
