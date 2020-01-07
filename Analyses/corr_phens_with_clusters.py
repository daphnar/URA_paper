import numpy as np
import os
import pandas
import glob
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from scipy.stats import ranksums
from scipy.stats import ks_2samp
from scipy.stats import ttest_ind
from scipy.stats import kruskal

CORR_TH = 5*10**-3
MIN_DATA_POINTS = 10

def corr_phens( th_bac_obj, bac, phen, out_path, th ):
    if not os.path.isfile(os.path.join(out_path, bac + "_corrs_phen.csv")):
        full_set = set(phen.StoolSample.values).intersection(th_bac_obj.index)
        tmp_phen = phen.set_index('StoolSample')
        tmp_phen = tmp_phen.loc[list(full_set)]
        print ("%d phen samples pass threshold coverage" % ( len(full_set) ))
        tmp_obj = th_bac_obj.loc[tmp_phen.index]
        d_r = {}
        inds_c = []
        inds_p = []
        for c in tmp_phen.columns:
            r = [{},{}]
            for col in tmp_obj.columns:
                try:
                    res = spearmanr(tmp_obj[col], tmp_phen[c], nan_policy='omit' )
                    r[0][col] = res[0]
                    r[1][col] = res[1]
                except:
                    r[0][col] = np.nan
                    r[1][col] = np.nan

            #        r = tmp_obj.apply(lambda col: spearmanr(col, tmp_phen[c], nan_policy='omit')[0], 0 )
    # pearson
    #        r = tmp_obj.corrwith(tmp_phen[c])
            d_r[c+"_cor"] = pandas.Series(r[0])
            d_r[c+"_pval"] = pandas.Series(r[1])
            inds_c.append( c+"_cor" )
            inds_p.append( c+"_pval")
        d_r = pandas.DataFrame(d_r).T
        d_r.to_csv(os.path.join(out_path, bac + "_corrs_phen.csv"))
    else:
        d_r = pandas.read_csv( os.path.join(out_path, bac + "_corrs_phen.csv"), index_col = 0 )
        inds_c = d_r.index[::2]
        inds_p = d_r.index[1::2]

    x = d_r.loc[inds_c].values.flatten()
    h = plt.hist(x, bins=100, range=[-1, 1])
    plt.savefig(os.path.join(out_path, bac + "_corrs_phen_hist.png") )

    res05 = {}
    for i,ind in enumerate(inds_p):
        tmp = d_r.loc[ind]
        tmp = tmp[tmp<th]
        for j in tmp.index:
            res05[(ind[:-5], j)] = [ tmp[j], d_r.loc[inds_c[i]][j]]
    if len(res05)==0:
        print("No corrs passed threshold")
        return
    res05=pandas.DataFrame(res05).T
    res05.columns = ['p_val', 'corr']
    res05.sort_values('p_val',inplace=True)
    if len(res05) < 20:
        print ("Top spearman corrs:")
        print ( res05 )
    else:
        print ("To many correlation (%d) abouve p_value %g, error in p-values of %d samples on %d parts" % \
               ( len(res05), th, len(th_bac_obj), len(th_bac_obj.columns )))
    res05.to_csv( os.path.join(out_path, bac + "_top_corrs_phen.csv") )
    return

def clust_2_phens( clusts, bac, phen_obj, out_path, th, stat_test ):
    if not os.path.isfile(os.path.join(out_path, bac + "_clust_2_phen.csv")):
        full_set = set(phen_obj.index).intersection(clusts.index)
        tmp_phen = phen_obj.loc[list(full_set)]
        print ("%d samples: 2 clusters vs %d phens" % ( len(full_set),  len(tmp_phen.columns )))
        d_r = {}
        inds_p = []
        no_check = []
        num_tests = 0
        for phen in tmp_phen.columns:
            phen_vals = tmp_phen[phen][~tmp_phen[phen].isna()]
            tmp_clusts = clusts.loc[phen_vals.index]
            if ( len(phen_vals[tmp_clusts==0]) < MIN_DATA_POINTS ) or ( len(phen_vals[tmp_clusts==1]) < MIN_DATA_POINTS ):
                no_check.append( phen )
                continue
            try:
                if stat_test == ttest_ind:
                    res = stat_test(phen_vals[tmp_clusts == 0], phen_vals[tmp_clusts == 1], equal_var=False)[1]
                else:
                    res=stat_test(phen_vals[tmp_clusts==0], phen_vals[tmp_clusts==1])[1]
                d_r[phen + "_grp0_size"] = len(phen_vals[tmp_clusts == 0])
                d_r[phen + "_grp1_size"] = len(phen_vals[tmp_clusts == 1])
                d_r[phen + "_pval"] = res
                inds_p.append(phen + "_pval")
                num_tests += 1
            except:
                print ("Exception for %s" % phen )
        print ("Didn't check %s for lack of data" % str(no_check))
        if len(d_r) != 0:
            d_r = pandas.Series(d_r)
            d_r.to_csv(os.path.join(out_path, bac + "_clust_2_phen.csv"))
        else:
            pandas.Series().to_csv(os.path.join(out_path, bac + "_clust_2_phen.csv"))
            print ("Nothing checked")
            return 0
    else:
        try:
            d_r = pandas.Series.from_csv( os.path.join(out_path, bac + "_clust_2_phen.csv"))
        except pandas.errors.EmptyDataError:
            print ("Nothing checked")
            return 0
        inds_p = d_r.index[2::3]
        num_tests = len(inds_p)

    pass_pval = {}
    tmp = d_r.loc[inds_p]
    tmp = tmp.astype('float')
    print ("Minimal p_val was %g" % tmp.min())
    tmp = tmp[tmp<th]
    for j in tmp.index:
        pass_pval[j[:-5]] = [ tmp[j], d_r.loc[j[:-4]+"grp0_size"], d_r.loc[j[:-4]+"grp1_size"]]
    if len(pass_pval)==0:
        print("No clust2 %s passed threshold" % stat_test.__name__)
        return num_tests
    pass_pval=pandas.DataFrame(pass_pval).T
    pass_pval.columns = ['p_val', 'size_g0', 'size_g1']
    pass_pval.sort_values('p_val',inplace=True)
    if len(pass_pval) < 5:
        print ("Top %s res:" % stat_test.__name__)
        print ( pass_pval )
    else:
        print ("To many correlation (%d) below p_value %g, for %s 2 clusters" % \
               ( len(pass_pval), th, bac ))
    pass_pval.to_csv( os.path.join(out_path, bac + "_top_clust_2_phen.csv") )
    return num_tests

def clust_n_phens( clusts, bac, phen_obj, out_path, th, stat_test ):
    if not os.path.isfile(os.path.join(out_path, bac + "_clust_n_phen.csv")):
        full_set = set(phen_obj.index).intersection(clusts.index)
        tmp_phen = phen_obj.loc[list(full_set)]
        min_cols = clusts.columns.astype(int).min()
        max_cols = clusts.columns.astype(int).max()
        print ("%d samples: %d-%d clusters vs %d phens" % ( len(full_set),  min_cols, max_cols, len(tmp_phen.columns)))
        d_r = {}
        inds_p = []
        no_check = {}
        num_tests = 0
        for phen in tmp_phen.columns:
            phen_vals = tmp_phen[phen][~tmp_phen[phen].isna()]
            tmp_clusts = clusts.loc[phen_vals.index]
            for c in range(min_cols,max_cols+1):
                phen_grps = []
                lens = []
                for i in range(c):
                    if (len(phen_vals[tmp_clusts[str(c)] == i]) >= MIN_DATA_POINTS):
                        phen_grps.append(phen_vals[tmp_clusts[str(c)] == i])
                        lens.append(len(phen_grps[-1]))
                if len(phen_grps) < 2:
                    if phen in no_check.keys():
                        no_check[phen].append( c )
                    else:
                        no_check[phen] = [c]
                    continue
                try:
                    res=stat_test( *phen_grps )[1]
                    d_r[phen + "_%d_grp_sizes" % c] = lens
                    d_r[phen + "_%d_pval" % c ] = res
                    inds_p.append(phen + "_%d_pval" % c )
                    num_tests += 1
                except:
                    print ("Exception for %s" % phen )
        print ("Didn't check %s for lack of data" % str(no_check))
        if len(d_r) != 0:
            d_r = pandas.Series(d_r)
            d_r.to_csv(os.path.join(out_path, bac + "_clust_n_phen.csv"))
        else:
            pandas.Series().to_csv(os.path.join(out_path, bac + "_clust_n_phen.csv"))
            print ("Nothing checked")
            return 0
    else:
        try:
            d_r = pandas.Series.from_csv(os.path.join(out_path, bac + "_clust_n_phen.csv"))
        except pandas.errors.EmptyDataError:
            print ("Nothing checked")
            return 0
        inds_p = d_r.index[1::2]
        num_tests = len(inds_p)

    pass_pval = {}
    tmp = d_r.loc[inds_p]
    tmp = tmp.astype('float')
    print ("Minimal p_val was %g" % tmp.min())
    tmp = tmp[tmp<th]
    for j in tmp.index:
        pass_pval[j[:-5]] = [ tmp[j], d_r.loc[j[:-4]+"grp_sizes"]]
    if len(pass_pval)==0:
        print("No clust_n %s passed threshold" % stat_test.__name__)
        return num_tests
    pass_pval=pandas.DataFrame(pass_pval).T
    pass_pval.columns = ['p_val', 'grp_sizes']
    pass_pval.sort_values('p_val',inplace=True)
    if len(pass_pval) < 5:
        print ("Top %s res:" % stat_test.__name__)
        print ( pass_pval )
    else:
        print ("To many correlation (%d) below p_value %g, for %s n clusters" % \
               ( len(pass_pval), th, bac ))
    pass_pval.to_csv( os.path.join(out_path, bac + "_top_clust_n_phen.csv") )
    return num_tests

def run ():
    stat_test = ttest_ind #ranksums # ks_2samp
    print ("Running test %s" % stat_test.__name__ )

    divide_to = 0  # 0 - by knee algorithm

    in_path = "/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Analyses/MetaPhlan100Uniq/bacteria/data"
    if divide_to == 0:
        out_path = "/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Analyses/MetaPhlan100Uniq/bacteria/cluster_vs_phen_%s" \
                   % stat_test.__name__
    else:
        out_path = "/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Analyses/MetaPhlan100Uniq/bacteria/cluster%d_vs_phen_%s" \
                   % ( divide_to, stat_test.__name__ )
    if not os.path.isdir(out_path):
        os.makedirs(out_path)

    phen = pandas.read_csv("/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/datasets/" + \
                       "SamplesAndPhenotypes_Genotek_study1_throw_dups_of_participants.csv")
    phen['Weight'] = phen['BMI'] * phen['Height'] * phen['Height'] / 100 / 100
    phen = phen.set_index('StoolSample')

    num_tests = 0
    if divide_to == 0:
        name = "_hierarchical_grps"
    else:
        name = "_hierarchical_%d_grps" % divide_to
    for filename in glob.glob(os.path.join(in_path, "*%s.csv" % name )):
        clst_bac_obj = pandas.read_csv( filename, index_col = 0)
        bac = filename.split('/')[-1].split(name)[0]
        print ("For %s %d participants and %d clustering divisions" % ( bac, len(clst_bac_obj), len(clst_bac_obj.columns)))
        if len(clst_bac_obj.columns) == 0:
            continue
        num_tests += clust_2_phens( clst_bac_obj['2'], bac, phen, out_path, CORR_TH, stat_test )
        cols = clst_bac_obj.columns.drop(['2'])
        if len(cols) > 0:
            num_tests += clust_n_phens( clst_bac_obj[cols], bac, phen, out_path, CORR_TH, kruskal )
    print ("Performed %d tests" % num_tests )

if __name__=='__main__':
    run()