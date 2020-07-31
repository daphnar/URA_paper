import pandas as pd
import os
basepath='/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Cohort_Paper/revision_Analyses/figures'

phenotypes = ['age','bmi','hba1c',"bt__fasting_glucose","bt__fasting_triglycerides",'bowel_movement_frequency',
              "bt__hdl_cholesterol","bt__cholesterol","bt__protein","bt__albumin","bt__sgot","bt__alkaline_phosphatase",
              "bt__tsh","bt__inr",'gender','currently_smokes','ever_smoked','type_2_diabetes']

def model_quantitive_rename(df):
    return df.rename(index=str,columns={'mean_pearson':'pearson_linear','std_pearson':'stdev_linear'})

def model_binary_rename(df):
    res = df.rename(index=str,columns={'auc':'auc_linear','stdev':'stdev_linear'})
    res = res[res['gender']=='all']
    return res[['target','n','auc_linear','stdev_linear']]

def gender_quantitive_rename(gender, df):
    if gender=='all':
        new_df = df.rename(index=str,columns={'cohort_size':'n',
                     'mean_pearson':'pearson','std_pearson':'stdev'})
    elif gender=='female':
        new_df = df.rename(index=str,columns={'cohort_size':'n_female',
                    'mean_pearson':'pearson_female','std_pearson':'stdev_female'})
    elif gender=='male':
        new_df = df.rename(index=str,columns={'cohort_size':'n_male',
                    'mean_pearson':'pearson_male','std_pearson':'stdev_male'})
    return new_df

def gender_binary_rename(gender, df):
    if gender == 'all':
        new_df = df[['target','auc','stdev']]
    elif gender=='female':
        new_df = df.rename(index=str,columns={'auc':'pearson_female','stdev':'stdev_female'})\
        .reset_index()[['target','pearson_female','stdev_female']]
        new_df.loc[new_df.shape[0]]=['gender',0,0]
    elif gender=='male':
        new_df = df.rename(index=str,columns={'auc':'pearson_male','stdev':'stdev_male'}) \
            .reset_index()[['target','pearson_male', 'stdev_male']]
        new_df.loc[new_df.shape[0]]=['gender',0,0]
    return new_df

def merger_quantitive():
    obj_ridge = pd.read_csv(os.path.join(basepath,'prediction_Ridge_linear.csv'),index_col=0)
    dfs=[]
    obj_xgboost = pd.read_csv(os.path.join(basepath,'prediction_xgboostMB_gender.csv'),index_col=0)
    dfs.append(model_quantitive_rename(obj_ridge))
    for gender,df in obj_xgboost.groupby('gender'):
        dfs.append(gender_quantitive_rename(gender,df))
    merge_on_target(dfs,os.path.join(basepath, 'Figures - prediction_other_phenotypes.csv'))

def merger_binary():
    obj_ridge = pd.read_csv(os.path.join(basepath,'classification_sgdMB_phenotypes.csv'),index_col=0)
    dfs=[]
    obj_xgboost = pd.read_csv(os.path.join(basepath,'classification_xgboostMB_phenotypes.csv'),index_col=0)
    for gender, df in obj_ridge.groupby('gender'):
        if gender=='all':
            dfs.append(model_binary_rename(obj_ridge))
    for gender,df in obj_xgboost.groupby('gender'):
        dfs.append(gender_binary_rename(gender,df))
    merge_on_target(dfs,os.path.join(basepath, 'Figures - classification.csv'))


def merge_on_target(dfs,output):
    merged = pd.DataFrame({'target':phenotypes})
    for df in dfs:
        print(merged)
        merged = merged.merge(df,on='target')
    merged.to_csv(output)

#merger_quantitive()
merger_binary()
