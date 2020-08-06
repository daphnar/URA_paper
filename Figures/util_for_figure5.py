import pandas as pd
import os
basepath='/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Cohort_Paper/revision_Analyses/figures'

data=pd.read_csv(os.path.join(basepath,'il_us_spear_saturation.csv'),index_col=0)
data.columns=["target","spearman_r","size","std"]
data.set_index('target').to_csv(os.path.join(basepath,'Figures - il_vs_us_spearman_saturation.csv'))