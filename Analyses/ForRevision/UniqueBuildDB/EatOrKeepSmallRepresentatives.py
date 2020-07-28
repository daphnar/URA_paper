import configparser
import os
import pandas as pd
import numpy

unicorn_path = '/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/'

# def eatOrKeepSGBs(rep_to_sgb_path,distanceFromLargeThreshold=0.075,distanceFromSmallThreshold=0.075):
#     allGenomes=pd.read_csv(os.path.join(unicorn_path, 'Segata/SupplementaryTable8-SGBsDescription.csv'),index_col=0)
#     SGBdescription=pd.read_csv(build_small_representatives.SGBdescription).set_index('SGB ID')
#     mash_path = os.path.join(unicorn_path, 'Bowtie_Segata/mash_4930/top_rep_4930.csv')
#
#     rep_to_sgb = pd.read_csv(rep_to_sgb_path, delim_whitespace=True, index_col=0, header=None)[1]
#     mash_df = pd.read_csv(mash_path, index_col=0)
#     new_names = [rep_to_sgb.loc[rep.replace('.fa', '')] for rep in mash_df.index]
#     mash_df.columns = new_names
#     mash_df.index = new_names
#     sgbCounts=allGenomes['# Reconstructed genomes']
#     largeSGBs=sgbCounts[sgbCounts>4].index.values
#     print ("We have %s large SGBs"%len(largeSGBs))
#     smallSGBcounts=sgbCounts.drop(largeSGBs)
#     tooClose=[]
#     for sgb in smallSGBcounts.index:
#         closestSGB=mash_df.loc[sgb,largeSGBs].min()
#         if closestSGB<distanceFromLargeThreshold:
#             tooClose.append(sgb)
#     smallSGBcounts = smallSGBcounts.drop(tooClose)
#     groupedSGBs={k:[] for k in smallSGBcounts.index}
#     annotation={k:[SGBdescription.loc[k,'Estimated taxonomy'].split('|')[-1]] for k in smallSGBcounts.index}
#     res = pd.DataFrame(index=smallSGBcounts.index,columns=['Assigned genomes','SGBs','annotation'])
#     res['Assigned genomes']=smallSGBcounts.values
#     mashIter=mash_df.loc[res.index, res.index].min().copy()
#     while mashIter.min()<distanceFromSmallThreshold:
#         smallestID1 = mashIter.idxmin()
#         mashCouples = mash_df.loc[res.index, res.index].idxmin()
#         smallestID2 = mashCouples.loc[smallestID1]
#         n1=smallSGBcounts.loc[smallestID1]
#         n2=smallSGBcounts.loc[smallestID2]
#         if n1>n2:
#             res.loc[smallestID1,'Assigned genomes']+=n2
#             groupedSGBs[smallestID1]=groupedSGBs[smallestID1]+groupedSGBs[smallestID2]+[smallestID2]
#             annotation[smallestID1]=annotation[smallestID1]+annotation[smallestID2]
#             res=res.drop(smallestID2)
#         elif n2>n1:
#             res.loc[smallestID2, 'Assigned genomes'] += n1
#             groupedSGBs[smallestID2] = groupedSGBs[smallestID1] + groupedSGBs[smallestID2] + [smallestID1]
#             annotation[smallestID2]=annotation[smallestID1]+annotation[smallestID2]
#             res = res.drop(smallestID1)
#         else:
#             if mash_df.loc[smallestID1,largeSGBs].min()<mash_df.loc[smallestID2,largeSGBs].min():
#                 res.loc[smallestID2, 'Assigned genomes'] += n1
#                 groupedSGBs[smallestID2] = groupedSGBs[smallestID1] + groupedSGBs[smallestID2] + [smallestID1]
#                 annotation[smallestID2] = annotation[smallestID1] + annotation[smallestID2]
#                 res = res.drop(smallestID1)
#             else:
#                 res.loc[smallestID1, 'Assigned genomes'] += n2
#                 groupedSGBs[smallestID1] = groupedSGBs[smallestID1] + groupedSGBs[smallestID2] + [smallestID2]
#                 annotation[smallestID1] = annotation[smallestID1] + annotation[smallestID2]
#                 res = res.drop(smallestID2)
#         mashIter = mash_df.loc[res.index, res.index].min().copy()
#     res['SGBs']=[groupedSGBs[sgb] for sgb in res.index]
#     res['annotation']=[annotation[sgb] for sgb in res.index]
#     return res,list(res.index)+list(largeSGBs)

def eatOrKeepDistantGenus(rep_to_sgb_path,SGBdescription_path):
    allGenomes=pd.read_csv(SGBdescription_path,index_col=0)
    GenusInfo = pd.read_csv(os.path.join(unicorn_path,'Segata/strain_taxonomy.csv'),index_col=0)
    #SGBdescription=pd.read_csv(build_small_representatives.SGBdescription).set_index('SGB ID')
    mash_path = os.path.join(unicorn_path, 'Bowtie_Segata/mash_4930/top_rep_4930.csv')
    #rep_to_sgb_path = os.path.join(unicorn_path, 'Segata/reconstructed_REPRESENTATIVE.txt')
    rep_to_sgb = pd.read_csv(rep_to_sgb_path, delim_whitespace=True, index_col=0, header=None)[1]
    mash_df = pd.read_csv(mash_path, index_col=0)
    new_names = [rep.replace('.fa', '') for rep in mash_df.index]
    mash_df.columns = new_names
    mash_df.index = new_names
    mash_df=mash_df.loc[rep_to_sgb.index,rep_to_sgb.index]
    new_names = [rep_to_sgb.loc[rep] for rep in mash_df.index]
    mash_df.columns = new_names
    mash_df.index = new_names
    sgbCounts=allGenomes.loc[new_names,'# Reconstructed genomes']
    largeSGBs=sgbCounts[sgbCounts>4].index.values
    print ("We have %s large SGBs"%len(largeSGBs))
    genusOfSGBs=GenusInfo.loc[rep_to_sgb.index]
    genusOfSGBs=genusOfSGBs.set_index('SGB')
    generaOfLargeSGBs=genusOfSGBs.loc[largeSGBs,'genus']
    smallSGBcounts=sgbCounts.drop(largeSGBs)
    tooClose=[]
    count=0
    for sgb in smallSGBcounts.index:
        if genusOfSGBs.loc[sgb,'genus'] in generaOfLargeSGBs.values:
            tooClose.append(sgb)
        else:
            count+=1
    print ("Added %s SGBs as they come from new genera"%count)
    smallSGBcounts = smallSGBcounts.drop(tooClose)
    annotation={k:[allGenomes.loc[k,'Estimated taxonomy'].split('|')[-1]] for k in smallSGBcounts.index}
    res = pd.DataFrame(index=smallSGBcounts.index,columns=['Assigned genomes','genus','annotation','minMashFromLarge','groupedSGBs','Grouped annotation'])
    res['Assigned genomes']=smallSGBcounts.values
    res['genus']=genusOfSGBs.loc[smallSGBcounts.index,'genus']
    res['annotation'] = [annotation[sgb] for sgb in res.index]
    res['minMashFromLarge']=mash_df.loc[smallSGBcounts.index, largeSGBs].min(axis=1)
    keepSGBs=[]
    keepSGBsAnnotations=[]
    groupedSGBs=[]
    for genus,genus_df in res.groupby('genus'):
        maxAssignedGenomes=genus_df['Assigned genomes'].max()
        maxOptions=genus_df[genus_df['Assigned genomes']==maxAssignedGenomes].index.values
        location = res.loc[maxOptions, 'minMashFromLarge'].argmax()
        keepID = res.loc[maxOptions,'minMashFromLarge'].index[location]
        keepSGBsAnnotations.append(genus_df['annotation'].values)
        keepSGBs.append(keepID)
        groupedSGBs.append(genus_df.index.values)
    res = res.loc[keepSGBs]
    res.loc[keepSGBs,'Grouped annotation'] = keepSGBsAnnotations
    res.loc[keepSGBs,'groupedSGBs'] = groupedSGBs
    return res,list(res.index)+list(largeSGBs)

def run(configFile):
    config = configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation())
    config.read(configFile)
    LargeOrNew = config['LargeOrNewGenus']
    SelectedSGBs,allSGBs = eatOrKeepDistantGenus(LargeOrNew['representatives'],LargeOrNew['SGBdescription'])
    pd.Series(allSGBs).to_csv(LargeOrNew['all_large_or_new_sgbs'],header=False,index=False)

if __name__=='__main__':
    run()