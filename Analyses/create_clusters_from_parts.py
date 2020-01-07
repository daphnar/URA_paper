import numpy as np
import os
import pandas
import glob
import time
import seaborn as sns
from scipy.cluster.hierarchy import linkage
import matplotlib.pyplot as plt
from Analyses.getStrainClustersFromCoverage import getRangeGroups

def xorsum( x, y ):
    return (x ^ y).sum()

def run_bac_clusterring( filename, bac, plot_path, data_path, k = 0 ):
    print ("Start %s" % bac, time.ctime())
    th_bac_obj = pandas.read_csv( filename, index_col = 0)
    res = {}
    for i,c in enumerate(th_bac_obj.index):
        if (i % 100) == 0:
              print (i, time.ctime())
        res[c] = th_bac_obj.apply( xorsum, 1, args = (th_bac_obj.loc[c],))
    res = pandas.DataFrame(res)
    res = res[res.index]
    if len(res)==1:
        print ("For %s got %d * %d (%d)" % (bac, len(th_bac_obj),len(th_bac_obj.columns), len(res)))
        return 0
    print ( bac, len(res))
    z = linkage(res, method='average', metric='euclidean')
    sns.clustermap(res, row_linkage=z, col_linkage=z)
#    z2 = sns.clustermap(res)
    plt.suptitle("Manhatten distance between samples,\nof %s parts existance" % bac )
    plt.savefig( os.path.join(plot_path, bac + "_grps.png"))
    plt.close('all')
    grps = getRangeGroups( z=z, k=k )
    if len(grps) == 0:
        return 0
    grps = grps.set_index(res.index)
    if k == 0:
        grps.to_csv( os.path.join(data_path, bac + "_hierarchical_grps.csv") )
    else:
        grps.to_csv(os.path.join(data_path, bac + "_hierarchical_%d_grps.csv" % k ))

    plt.hist( th_bac_obj.sum(1), bins=max(int(len(th_bac_obj)/5),2))
    plt.title( "Number of non-del parts, of %d,\nfor %s" % ( len(th_bac_obj.columns), bac ))
    plt.ylabel("num samples")
    plt.xlabel("num parts")
    plt.savefig( os.path.join(plot_path, bac + "_hist_num_parts.png"))
    plt.close('all')
    print ("End %s" % bac, time.ctime())
    return len(grps.columns)

def run ():
#    inout_path = "/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Analyses/MetaPhlan100Uniq/bacteria/data"
    inout_path = "/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Analyses/MetaPhlan100Uniq/bacteria_tmp/data_non0_cov3"
    out_path = "/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Analyses/MetaPhlan100Uniq/bacteria_tmp/plots_non0_cov3"
    if not os.path.isdir(out_path):
        os.makedirs(out_path)

    divide_to = 0 # 0 - by knee algorithm
    cnt = 0
    num_clusts = 0
    for filename in glob.glob(os.path.join(inout_path, "*_7_call_part_existance.csv")):
        th_bac_obj = pandas.read_csv( filename, index_col = 0)
        bac = filename.split('/')[-1].split('_7_call_part_existance')[0]
        print ("For %s %d participants and %d representative parts" % ( bac, len(th_bac_obj), len(th_bac_obj.columns)))
        if len(th_bac_obj.columns) == 0:
            continue
        num_clusts += run_bac_clusterring( filename, bac, out_path, inout_path, divide_to )
        cnt += 1
    print ("done. Got %d clusters for %d bacterias" % ( num_clusts, cnt ), time.ctime())


if __name__=='__main__':
    run()