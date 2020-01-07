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

def plot_scores( gnms, norm, in_path, out_path ):
    for gnm in gnms:
        fs = glob.glob( os.path.join( in_path, "*_%d_*.csv" % gnm ))
        if len(fs) != 2:
            print ( "Note: %d files (not 2) for genome %d, skipping" % ( len(fs), gnm ))
            continue
        if ('noise' in fs[0]):
            fs = [ fs[1], fs[0] ]

        df1=pandas.read_csv(fs[0])
        fname = fs[0].split('/')[-1].split('.')[0]
        if os.path.isfile( os.path.join( out_path, fname + ".png" )):
            print ("Outfile for genome %d exists, skipping" %  gnm)
            continue
        grp = int(fname.split('_')[2])*10
        tmp = df1[(df1.genome!=gnm)&(df1['count']>0)]
        tmp = tmp[(tmp.result!=1)&(tmp.result!=grp)]
        print ("For %s %d out fo genome/grp/root (sum %g)" % ( fname, len(tmp), tmp['count'].sum()))
        df1=df1[(df1.genome==gnm)][['part','count']]
        df2=norm[norm.Strain==gnm][['PartNum','expected_score_window','num_unique']]
        df_m=pandas.merge( df1, df2, left_on='part',right_on='PartNum')
        cnt_exp = 0
        df_m['exp']=df_m.apply(add_exp,1)
        df1=pandas.read_csv(fs[1])
        fname2 = fs[1].split('/')[-1]
        tmp = df1[(df1.genome!=gnm)&(df1['count']>0)]
        tmp = tmp[(tmp.result!=1)&(tmp.result!=grp)]
        print ("For %s %d out fo genome/grp/root (sum %g)" % ( fname2, len(tmp), tmp['count'].sum()))
        df1=df1[(df1.genome==gnm)][['part','count']]
        df_m=pandas.merge( df1, df_m, on='part')
        df_m.dropna(inplace=True)
        df_m['m_x']=df_m['count_x']*df_m['exp']/1000
        df_m['m_y']=df_m['count_y']*df_m['exp']/1000
#        plt.figure()
        plt.scatter(df_m['exp'].values,df_m['m_x'].values,c='b',label='clean')
        plt.scatter(df_m['exp'].values,df_m['m_y'].values,c='r',label='noisy')
        plt.plot([df_m['exp'].min(),df_m['exp'].max()],[df_m['exp'].min(),df_m['exp'].max()],c='g')
        plt.legend(loc=2)
        plt.xlabel('expected score')
        plt.ylabel('actual score')
        plt.title(fname + "\nparts scoring")
        plt.savefig( os.path.join( out_path, fname + ".png" ))
        plt.close('all')

def run():
    path = "/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/databases/MetaPhlan100ANI"
    read_len = 75
    df_tree = pandas.read_csv( os.path.join( path, "Part",  "df_tree.csv" ))
    norm = pandas.read_csv(  os.path.join( path, "Uniq", "united_normalized_score.csv"))

    in_path = "/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Mappings/MetaPhlan100ANIUniq_debug/rawMappings"
    out_path = os.path.join( in_path, "plots" )
    if not os.path.isdir( out_path ):
        os.makedirs( out_path )
    check = list(range(10,1001,10))
    check.remove(310)
    check.remove(970)

    names = create_names( df_tree, check )
    gnms = names.keys()
    gnms.sort()
    plot_scores( gnms, norm, in_path, out_path )

if __name__=='__main__':
    run()
