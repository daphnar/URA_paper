import numpy as np
import os
import pandas
import glob
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from scipy.stats import ranksums
from scipy.stats import ks_2samp
from scipy.stats import ttest_ind

CORR_TH = 5*10**-5

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

def grp_exist_phens( th_bac_obj, bac, phen, out_path, th, stat_test ):
    if not os.path.isfile(os.path.join(out_path, bac + "_grp_exist_phen.csv")):
        full_set = set(phen.StoolSample.values).intersection(th_bac_obj.index)
        tmp_phen = phen.set_index('StoolSample')
        tmp_phen = tmp_phen.loc[list(full_set)]
        print ("%d samples: %d representative parts vs %d phens" % ( len(full_set),  len(th_bac_obj.columns), len(tmp_phen.columns )))
        tmp_obj = th_bac_obj.loc[tmp_phen.index]
        d_r = {}
        inds_p = []
        for phen in tmp_phen.columns:
            phen_vals = tmp_phen[phen][~tmp_phen[phen].isna()]
            r = [{},{},{}]
            for part_exist in tmp_obj.columns:
                try:
                    res=stat_test(phen_vals[tmp_obj[part_exist]==0], phen_vals[tmp_obj[part_exist]==1])
                    r[0][part_exist] = len(phen_vals[tmp_obj[part_exist]==0])
                    r[1][part_exist] = len(phen_vals[tmp_obj[part_exist]==1])
                    r[2][part_exist] = res[1]
                except:
                    r[0][part_exist] = np.nan
                    r[1][part_exist] = np.nan
                    r[2][part_exist] = np.nan
            d_r[phen+"_grp0_size"] = pandas.Series(r[0])
            d_r[phen+"_grp1_size"] = pandas.Series(r[1])
            d_r[phen+"_pval"] = pandas.Series(r[2])
            inds_p.append( phen+"_pval")
        d_r = pandas.DataFrame(d_r).T
        d_r.to_csv(os.path.join(out_path, bac + "_grp_exist_phen.csv"))
    else:
        d_r = pandas.read_csv( os.path.join(out_path, bac + "_grp_exist_phen.csv"), index_col = 0 )
        inds_p = d_r.index[2::3]

    x = d_r.loc[inds_p].values.flatten()
    x1 = x[~np.isnan(x)]
    if len(x1)<len(x):
        print ("Had %d nans in test" % (len(x)-len(x1)))
    h = plt.hist(x1, bins=100,log=True)
    plt.title("p_values")
    plt.savefig(os.path.join(out_path, bac + "_grp_exist_phen_p_values.png") )
    print ("save png")
    plt.close('all')

    pass_pval = {}
    for i,ind in enumerate(inds_p):
        tmp = d_r.loc[ind]
        tmp = tmp[tmp<th]
        for j in tmp.index:
            pass_pval[(ind[:-5], j)] = [ tmp[j], d_r.loc[ind[:-4]+"grp0_size"][j], d_r.loc[ind[:-4]+"grp1_size"][j]]
    if len(pass_pval)==0:
        print("No %s passed threshold" % stat_test.__name__)
        return
    pass_pval=pandas.DataFrame(pass_pval).T
    pass_pval.columns = ['p_val', 'size_g0', 'size_g1']
    pass_pval.sort_values('p_val',inplace=True)
    if len(pass_pval) < 20:
        print ("Top %s res:" % stat_test.__name__)
        print ( pass_pval )
    else:
        print ("To many correlation (%d) above p_value %g, for %s" % \
               ( len(pass_pval), th, bac ))
    pass_pval.to_csv( os.path.join(out_path, bac + "_top_grp_exist_phen.csv") )
    return

def run ():
    stat_test = ranksums # ttest_ind ks_2samp ranksums
    print ("Running test %s" % stat_test.__name__ )

#    in_path = "/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Analyses/MetaPhlan100Uniq/bacteria_genotek0_and_1/data"
    in_path = "/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Analyses/MetaPhlan100Uniq/bacteria/data"
    out_path = "/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Analyses/MetaPhlan100Uniq/bacteria/exist_vs_phen_%s" \
               % stat_test.__name__
    if not os.path.isdir(out_path):
        os.makedirs(out_path)

    phen = pandas.read_csv("/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/datasets/" + \
                       "SamplesAndPhenotypes_throw_dups_of_participants.csv")
    phen['Weight'] = phen['BMI'] * phen['Height'] * phen['Height'] / 100 / 100

    for filename in glob.glob(os.path.join(in_path, "*_8_for_corr_with_phen.csv")):
        th_bac_obj = pandas.read_csv( filename, index_col = 0)
        bac = filename.split('/')[-1].split('_8_for_corr_with_phen')[0]
        print ("For %s %d participants and %d representative parts" % ( bac, len(th_bac_obj), len(th_bac_obj.columns)))
        if len(th_bac_obj.columns) == 0:
            continue
        grp_exist_phens( th_bac_obj, bac, phen, out_path, CORR_TH, stat_test )

if __name__=='__main__':
    run()