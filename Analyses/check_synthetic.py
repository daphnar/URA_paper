import pandas
import os
import glob

def run():
    inout_path = "/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Analyses/MetaPhlan100ANIUniq/bacteria_synthetic/data"
    in_file_ext = "_meannon0_abnd"
    out_path = "/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Analyses/MetaPhlan100ANIUniq/bacteria_synthetic"
    if not os.path.isdir(out_path):
        os.makedirs(out_path)

    dfs = {}
    for filename in glob.glob(os.path.join(inout_path, "*%s.csv" % in_file_ext )):
        bac = filename.split('/')[-1].split(in_file_ext)[0]
        tmp = pandas.read_csv( filename, index_col = 0,header=None).T
        if tmp.loc[1].sum() > 0:
            dfs[bac] = tmp.loc[1]
    abnd = pandas.concat(dfs)
    abnd = abnd.unstack(-1)
    rel_abnd = abnd/abnd.sum()
    rel_abnd = rel_abnd.sort_index()
    rel_abnd.to_csv ( os.path.join( out_path, "abnd_synthetic.csv"))
    print ("done. Got %d clients for %d bacterias" % ( len(rel_abnd.columns), len(rel_abnd) ))


if __name__ == '__main__':
    run()