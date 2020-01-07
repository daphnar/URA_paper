import pandas
import os
import matplotlib.pyplot as plt
import glob
import numpy

def create_names( df_tree, check ):
    names={}
    for i in check:
        tmp_tree = df_tree[df_tree.parent==i]
        for ind in tmp_tree.index:
            names[df_tree.loc[ind].node] = [ i, tmp_tree.loc[ind]['name'] ]
    return names

cnt_exp = 0
def add_exp(x):
    global cnt_exp
    cnt_exp += x.expected_score_window
    if numpy.isnan(x['count']):
        return 0
    else:
        r = cnt_exp
        cnt_exp = 0
        return r

def check_scores( gnms, in_path ):
    res = {}
    for gnm in gnms:
        fs = glob.glob( os.path.join( in_path, "*_%d.csv" % gnm ))
        if len(fs) != 1:
            print ( "WTF. %d files (not 1) for genome %d, skipping" % ( len(fs), gnm ))
            continue
        f = fs[0]

        df=pandas.read_csv(f, index_col = 0)
        m = df.mean(1)
        print ( "for %d %s" % ( gnm, f.split("/")[0] ))
        if len(m[m>500]) != 2:
            print ("Too many genomes (%d) got above 500" %  len(m[m>500]) )
            for ind in m[m>500].index:
                print ( "%s %g (%s)" % ( ind.split('_')[3], m.loc[ind], ind ))
            continue
        if ( int(m[m>500].index[0].split('_')[3]) != gnm ) or ( int(m[m>500].index[1].split('_')[3]) != gnm ):
            print ( "Wrong genomes?!?" )
            continue
        else:
            if 'noise' in m[m>500].index[0]:
                print ("Real genome got mean %g with noise mean %g" % ( m[m>500].iloc[1], m[m>500].iloc[0] ))
                res[gnm] = [ m[m>500].iloc[1], m[m>500].iloc[0] ]
            else:
                print ("Real genome got mean %g with noise mean %g" % (m[m > 500].iloc[0], m[m > 500].iloc[1]))
                res[gnm] = [m[m > 500].iloc[0], m[m > 500].iloc[1]]
        other = m[m<500].index
        print ( "Max of mean of all others %g on %s" % ( m.loc[other].max(),  m.loc[other].idxmax() ))
        res[gnm] += [ m.loc[other].max(),  m.loc[other].idxmax() ]
        print ( "Max of max of all others %g on %s" % ( df.max(1).loc[other].max(), df.max(1).loc[other].idxmax() ))
        res[gnm] += [ df.max(1).loc[other].max(), df.max(1).loc[other].idxmax() ]
        x = df.loc[other].iloc[0]
        df['num_non_0'] = df.apply( count_non0, 1 )
        res[gnm] += [ df.loc[other]['num_non_0'].max(), df.loc[other]['num_non_0'].idxmax(), len(x[x>=0]) ]
        print ()
    res = pandas.DataFrame(res).T
    res.columns = [ 'mean_clean', 'mean_noise', 'max_mean', 'max_mean_on', 'max_max', 'max_max_on', 'max_non_0', \
                    'max_non_0_on', 'of_parts']
    res.sort_index( inplace = True )
    return res

def count_non0(x):
    return len(x[x>0])

def run():
    path = "/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/databases/MetaPhlan100ANI"
    read_len = 75
    df_tree = pandas.read_csv( os.path.join( path, "Part",  "df_tree.csv" ))

    out_path = "/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Mappings/MetaPhlan100ANIUniq_debug"
    in_path = "/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Mappings/MetaPhlan100ANIUniq_debug/bacteria"
    check = list(range(10,1001,10))
    check.remove(310)
    check.remove(970)

    names = create_names( df_tree, check )
    gnms = names.keys()
    gnms.sort()
    res = check_scores( gnms, in_path )
    res.to_csv( os.path.join( out_path, "res_debug.csv" ))

if __name__=='__main__':
    run()
