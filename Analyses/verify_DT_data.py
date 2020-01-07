import pandas
import os
import glob
import time
import pickle

path = "/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/datasets"

all_files = glob.glob( os.path.join( path, "SSUnitedFastq", "DT*.fastq" ))
all_bars = list(map(lambda x: x.split('DT')[-1].split('_')[0], all_files))

bars_df = []
no_bars_all = []
full_path = os.path.join(path, "DayTwo", "pheno" )
f = os.path.join( full_path, "for_wz_barcode_client_mapping.csv" )
df = pandas.read_csv(f)
bars = {}
no_bars = []
for i, bar in enumerate(df.barcode):
    bar = str(bar)
    if (i % 1000) == 0:
        print ("%d (%d found) of %d" % (i, len(bars), len(df)), time.ctime())
    if bar in all_bars:
        bars[bar] = all_files[all_bars.index(bar)]
    else:
        no_bars.append(bar)
print ("got: %d found analysed, %d not" % ( len(bars), len(no_bars)))
bars_df = pandas.Series(bars)
#    bars_df[-1].to_csv(os.path.join(path, "DayTwo", f.split('/')[-1][:-4] + "_found_bars.csv"))
bars_df.to_csv(os.path.join( full_path, "found_bars.csv"))

all_files = glob.glob( "/net/mraid08/export/genie/Data/NextSeq/*/FastqOutput/Unaligned_fastq/Sample_*" )
all_bars = list(map(lambda x: x.split('Sample_')[-1].split('_')[0], all_files))

not_in = []
found = {}
for i, bar in enumerate(no_bars):
    if (i % 1000) == 0:
        print ("%d (%d found, not analysed) of %d" % (i, len(found), len(no_bars)), time.ctime())
    bar = str(bar)
    if not(bar in all_bars):
        not_in.append(bar)
    else:
        found[bar] = all_files[all_bars.index(bar)]
print ("got: %d found analysed, %d fount not analysed, %d not found" % ( len(bars), len(found), len(not_in)))
found_bars = pandas.Series(found)
found_bars.to_csv( os.path.join( full_path, "found_non_analized_bars.csv" ))
tmp=pandas.Series(not_in)
tmp.to_csv(os.path.join( full_path, "not_found.csv"))
print ("done")