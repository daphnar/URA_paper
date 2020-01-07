import pandas
import os
import time
import random
from numpy.random import choice
from addloglevels import sethandlers
from SegalQueue.qp import qp,fakeqp

err_ps = None


def create_names( df_tree, check ):
    names={}
    for i in check:
        tmp_tree = df_tree[df_tree.parent==i]
        for ind in tmp_tree.index:
            names[df_tree.loc[ind].node] = [ i, tmp_tree.loc[ind]['name'] ]
    return names

def create_reads_files( fict_path, singles_path, names, read_len, p_noise = None ):
    print ("Initializing")
    basedir = "/net/mraid08/export/genie/Lab/Unicorn"

    sethandlers()
    os.chdir(basedir)
    print ("Starting")

    if os.path.isdir( fict_path ):
        print ( "Fict reads directory exists. Writing to it")
    else:
        os.mkdir( fict_path )
    i = [0,0]
    with qp(jobname='make_fict', q=['himem7.q'], tryrerun=True, mem_def ='2G',trds_def=2, \
                    qworker='~/Develop/Python/lib/SegalQueue/qworker.py') as q:
        q.startpermanentrun()
        waiton = []
        for dirpath, dnames, fnames in os.walk( singles_path ):
            # GCF_000281435_97_985260000_s__klebsiella_pneumoniae.fa
            for fname in fnames:
                i[1] += 1
                ps = fname.split("_")
                if not ( int(ps[3]) in names.keys()):
                    continue
                i[0] += 1
                waiton.append(q.method(run_create, (dirpath, fict_path, fname, names, read_len, p_noise)))
            print ("sent %d jobs (of %d possible)" % (i[0], i[1]))
            q.wait(waiton)
            break
    return


def run_create( dirpath, fict_path, fname, names, read_len, p_noise ):
    print ("Started genome %s" % ( fname ), time.ctime())
    fin = open( os.path.join( dirpath, fname ),"r")
    fs = []
    fs.append( open( os.path.join( fict_path, fname), "w"))
    if p_noise != None:
        if (len(p_noise) != 2 ) or ( sum(p_noise) > 0.5 ) or ( sum(p_noise) <= 0 ):
            raise ("WTF. noise values are wrong")
        fs.append(open(os.path.join(fict_path, fname[:-3] + "_noise_%g_%g.fa" %
                                    ( p_noise[0], p_noise[1] )), "w"))
    out = False
    seq = ""
    while True:
        l = fin.readline()
        if (len(l) == 0) or (l[0] == '>'):
            if seq != "":
                out_all_read_to_file( seq, fs, p_noise, save_name, read_len )
            if len(l) == 0:
                break
            taxid = int(l.split('|')[1])
            if taxid in names.keys():
                out = True
                pos = l.find('|')
                pos = l.find('|',pos+1)
                save_name = l[pos+1:-1]
                seq = ""
            else:
                out = False
                seq = ""
        elif out:
            seq += l[:-1]
    fin.close()
    for f in fs:
        f.close()
    print ("Finished genome %s" % ( fname ), time.ctime())
    return

# def add_noise( seq, l_seq, ps ):
#     global err_ps
#     if err_ps == None:
#         err_ps = {}
#         p = ps[1]
#         err_ps['indel']  = [ (1-p)**l_seq, (1-p)**(l_seq-1)*p*l_seq, \
#                                      (1-p)**(l_seq-2)*p**2*l_seq*(l_seq-1)/2]
#         err_ps['indel'].append( 1- sum(err_ps['indel']))
#         # so err could be to self, for efficiency
#         p = ps[0]*4./3
#         err_ps['err'] = [ (1-p)**l_seq, (1-p)**(l_seq-1)*p*l_seq, \
#                         (1-p)**(l_seq-2)*p**2*l_seq*(l_seq-1)/2, \
#                         (1 - p) ** (l_seq - 3) * p ** 3 * l_seq * (l_seq - 1) * (l_seq - 2)/ 6 ]
#         err_ps['err'].append( 1- sum(err_ps['err']))
#     n_seq = seq
#     try:
#         pr = random.random()
#         num_errs = 0
#         if pr > err_ps['err'][0]:
#             pr -= err_ps['err'][0]
#             if pr > err_ps['err'][1]:
#                 pr -= err_ps['err'][1]
#                 if pr > err_ps['err'][2]:
#                     pr -= err_ps['err'][2]
#                     if pr > err_ps['err'][3]:
#                         num_errs = 4
#                     else:
#                         num_errs = 3
#                 else:
#                     num_errs = 2
#             else:
#                 num_errs = 1
#         locs = choice(range(l_seq), num_errs, False)
#         locs.sort()
#         for l in locs[::-1]:
#             opts = ["A", "G", "C", "T"]
#             n_seq = n_seq[:l] + choice(opts, 1)[0] + n_seq[l + 1:]
#         pr = random.random()
#         num_indel = 0
#         if pr > err_ps['indel'][0]:
#             pr -= err_ps['indel'][0]
#             if pr > err_ps['indel'][1]:
#                 pr -= err_ps['indel'][1]
#                 if pr > err_ps['indel'][2]:
#                     num_indel = 3
#                 else:
#                     num_indel = 2
#             else:
#                 num_indel = 1
#         locs = choice(range(l_seq), num_indel, False)
#         locs.sort()
#         for l in locs[::-1]:
#             if random.random() < 0.5:
#                 n_seq = n_seq[:l] + n_seq[l + 1:]
#             else:
#                 n_seq = n_seq[:l] + choice(["A", "G", "C", "T"], 1)[0] + n_seq[l:]
#         return n_seq[:l_seq]
#     except:
#         print ("exception occured on %s %d %d" % ( seq, len(seq), len(n_seq)))
#         return seq[:l_seq]

def add_noise( seq, l_seq, ps ):
    n_seq = ""
    for aa in seq:
        pr = random.random()
        if pr < ps[1]/2:
            continue
        if pr < ps[1]:
            n_seq += choice(["A", "G", "C", "T"], 1)[0]
        pr -= ps[1]
        if pr < ps[0]*4./3:
            n_seq += choice(["A", "G", "C", "T"], 1)[0]
        else:
            n_seq += aa
    return n_seq[:l_seq]

def out_all_read_to_file( seq, fs, p_noise, save_name, len_read ):
    cnt = 0
    for i in range(len(seq)-len_read+1):
        for ind in range(len(fs)):
            fs[ind].write((">%s|pos%d\n" % (save_name, i)))
            if ind == 0:
                fs[ind].write( seq[i:i+len_read] +'\n')
            else:
                tmp = add_noise(seq[i:i + len_read + 10], len_read, p_noise )
                if tmp != seq[i:i + len_read]:
                    cnt += 1
                fs[ind].write( tmp + '\n')
#    print ("Done", i, cnt, time.ctime() )

def run():
    path = "/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/databases/MetaPhlan100ANI"
    in_path = os.path.join( path, "single_genomes" )
    read_len = 75
    df_tree = pandas.read_csv( os.path.join( path, "Part",  "df_tree.csv" ))
    out_path = os.path.join( path, "fict_reads_%d" % read_len )
#    check = [ 10, 20, 30 ]
    check = list(range(40,1001,10))
    check.remove(310)
    check.remove(970)

    p_noise = [0.005, 0.002] # None

    names = create_names( df_tree, check )
    create_reads_files( out_path, in_path, names, read_len, p_noise )

if __name__=='__main__':
    run()


