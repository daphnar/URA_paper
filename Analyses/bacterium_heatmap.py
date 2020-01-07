import numpy as np
import glob
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from Mappings.config import prepare_alternative_mapping_dataframes
output_path='/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Analyses/'+\
            'MetaPhlan100Uniq/plots/bacteria'
# bac_outputs=prepare_alternative_mapping_dataframe.output_bacteria_path
if not os.path.exists(output_path):
    os.makedirs(output_path)

bacs=["s__eubacterium_hallii_GCF_000173975_genome_51000000.csv",
"s__roseburia_inulinivorans_GCF_000174195_genome_39000000.csv",
"s__coprococcus_catus_GCF_000210555_genome_88000000.csv",
"s__roseburia_hominis_GCF_000225345_genome_37000000.csv",
"s__coprococcus_comes_GCF_000155875_genome_56000000.csv",
"s__bacteroides_faecis_GCF_000226135_genome_46000000.csv",
"s__ruminococcus_torques_GCF_000210035_genome_33010000.csv",
"s__bacteroidales_bacterium_ph8_GCF_000311925_genome_41000000.csv",
"s__faecalibacterium_prausnitzii_GCF_000166035_genome_4000000.csv"]
for filename in bacs:
    print filename
    bac=filename.replace('.csv','')
    bac_outputs=os.path.join(prepare_alternative_mapping_dataframes.output_bacteria_path,
                             filename)
    obj = pd.read_csv(bac_outputs,index_col=0)
    percentileCap=99
    x=obj.values.flatten()
    x = x[~np.isnan(x)]
    cap=np.percentile(x,percentileCap)
    objCapped=obj.applymap(lambda x: min(x,cap))
    objCapped.columns = [int(col)%10000 for col in objCapped.columns]
    plt.figure()
    sns.heatmap(objCapped,xticklabels=500, yticklabels=False)
    plt.ylabel('Samples')
    plt.xlabel('Unique part number')
    plt.title(bac+' capped at %.2f'%cap)
    plt.savefig(os.path.join(output_path,bac+".png"))
    plt.close('all')