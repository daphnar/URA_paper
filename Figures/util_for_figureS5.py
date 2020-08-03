import pandas as pd
import os
basepath='/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Cohort_Paper/revision_Analyses/figures'

def merger(basepath):
    obj1 = pd.read_csv(os.path.join(basepath,'IL_atlas.csv'))
    obj2 = pd.read_csv(os.path.join(basepath,'US_atlas.csv'))
    merged = obj1.merge(obj2,suffixes=['','_us'],on=['pheno','species'])
    merged.to_csv(os.path.join(basepath,'Figures - supplementary_il_vs_us_sgb_pvals.csv'))

merger(basepath)

