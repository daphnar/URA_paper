import numpy as np
from scipy.cluster.hierarchy import linkage,cut_tree,dendrogram
import matplotlib.pyplot as plt
import pandas
import os
bacteria_path='/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Analyses/'+\
            'MetaPhlan100Uniq/bacteria/data'
bacteria_output_plot='/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Analyses/'+\
            'MetaPhlan100Uniq/bacteria/'
bacterium_coverage_file='s__anaerostipes_hadrus_GCF_000332875_genome_100000000_corr95_cols.csv'

def getGroups(z,k,out_groups = True):
    n_participants = len(z) + 1
    groups = []
    cluster_ids = cut_tree(z, k).flatten()
    if not out_groups:
        return cluster_ids
    all_IDs = pandas.Series(range(n_participants))
    for group in range(k):
        groups.append(all_IDs[cluster_ids == group].values)
    return np.array(groups)

def getRangeGroups(z, k = 0):
    # if k = 0 decide on number of groups by knee algorithm
    # else create k groups
    if k == 0:
        k = getNumberOfGroups(z)
    groups = {}
    for i in range(2,k+1):
        groups[i] = getGroups(z, i, False)
    return pandas.DataFrame( groups )

def createLinkage( bacCoverage_file,plot=False):
    obj  = pandas.read_csv(bacCoverage_file,index_col=0)
    z = linkage(obj,method='average',metric='euclidean')
    if plot:
        plt.figure()
        dendrogram(z)
        plt.show()
    return z


def getNumberOfGroups(z, plot=False):
    last = z[:, 2]
    last_rev = last[::-1]
    idxs = np.arange(1, len(last) + 1)
    acceleration = np.diff(last, 2)  # 2nd derivative of the distances
    acceleration_rev = acceleration[::-1]
    if plot:
        plt.figure()
        plt.plot(idxs, last_rev)
        plt.plot(idxs[:-2] + 1, acceleration_rev)
        plt.xlabel('Difference number')
        plt.ylabel('Difference')
        plt.show()
    try:
        k = acceleration_rev.argmax() + 2
    except:
        print ("Failed to cluster")
        return 0
    print "clusters:", k
    return k

def heirarchicalCluster( bacCoverage_file,plot=False ):
    z = createLinkage( bacCoverage_file,plot)
    k = getNumberOfGroups(z, plot)
    return getGroups(z, k)

if __name__=="__main__":
    print getGroups(os.path.join(bacteria_path, bacterium_coverage_file))