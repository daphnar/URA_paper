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
for f in glob.glob(os.path.join(path, "DayTwo", "*_barcode_list.csv")):
    f = os.path.join( full_path, "for_wz_barcode_client_mapping.csv" )
    df = pandas.read_csv(f, index_col=0)
    bars = {}
    no_bars = []
    for i, bar in enumerate(df.barcode):
        bar = str(bar)
        if (i % 1000) == 0:
            print ("%d (%d good) of %d" % (i, len(bars), len(df)), time.ctime())
        if bar in all_bars:
            bars[bar] = all_files[all_bars.index(bar)]
        else:
            no_bars.append(bar)
    print ("for %s got: %d good, %d bad" % (f.split('/')[-1], len(bars), len(no_bars)))
    bars_df.append(pandas.Series(bars))
    bars_df[-1].to_csv(os.path.join(path, "DayTwo", f.split('/')[-1][:-4] + "_found_bars.csv"))
    no_bars_all.append(no_bars)
df = pandas.concat(bars_df)
df.to_csv(os.path.join(path, "DayTwo", "all_found_bars.csv"))
pickle.dump(no_bars_all, open(os.path.join(path, "DayTwo", "no_bars.pkl"), "wb"))

# all_samples = pandas.read_csv( os.path.join( path, "DayTwo", "samples follow-up table - v3 automated.csv"))
# all_samples=all_samples[['Sample ID \n(barcode)','Collection date']]
# all_samples.columns=['barcode','date']
#
# not_in = []
# dates = {}
# for no_bars in no_bars_all:
#     for bar in no_bars:
#         tmp = all_samples[all_samples['barcode'] == bar]
#         if len(tmp)==0:
#             not_in.append(bar)
#         elif len(tmp)>1:
#             print ("WTF %s" % bar )
#             continue
#         else:
#             date = tmp.date.iloc[0]
#             if type(date) == str:
#                 dates[bar] = [ date, date[::-1] ]
#             else:
#                 dates[bar] = [ date, date ]
# dates2 = pandas.DataFrame(dates).T

all_files = glob.glob( "/net/mraid08/export/genie/Data/NextSeq/*/FastqOutput/Unaligned_fastq/Sample_*" )
all_bars = list(map(lambda x: x.split('Sample_')[-1].split('_')[0], all_files))

not_in = []
found = {}
for no_bars in no_bars_all:
    for bar in no_bars:
        bar = str(bar)
        if not(bar in all_bars):
            not_in.append(bar)
        else:
            found[bar] = all_files[all_bars.index(bar)]
found_bars = pandas.Series(found)
found_bars.to_csv( os.path.join( path, "DayTwo", "found_non_analized_bars.csv" ))
tmp=pandas.Series(not_in)
tmp.to_csv(os.path.join( path, "DayTwo", "not_found.csv"))
print ("done")