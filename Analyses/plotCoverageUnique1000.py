import matplotlib.pyplot as plt
import numpy  as np
import pandas
import os
import glob

#samplesPath='/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Analyses/mappingUnique1000Samples'
#plotDir='/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Analyses/plots/mappedSamplesToUnique1000'

#samplesPath='/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Analyses/mappingUnique1000Samples_MetaphlanByParts'
samplesPath='/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Analyses/mappingUnique1000Samples_MetaphlanBySinglePart_thresh200/'
plotDir='/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Analyses/plots/mappedSamplesToUnique1000_MetaphlanBySingleParts_thresh200'
unicornCoverageMedianDir='/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Analyses/mappingUnique1000Samples_MetaphlanBySinglePart_thresh200/unicornMedianCoverage'
unicornCoverageMedianWithZeroDir='/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Analyses/mappingUnique1000Samples_MetaphlanBySinglePart_thresh200/unicornMedianWithZeroCoverage'
unicornCoverageMeanDir='/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Analyses/mappingUnique1000Samples_MetaphlanBySinglePart_thresh200/unicornMeanCoverage'
unicornCoverageMeanWithZeroDir='/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Analyses/mappingUnique1000Samples_MetaphlanBySinglePart_thresh200/unicornMeanWithZeroCoverage'
unicornCoverageByGroup='/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Analyses/mappingUnique1000Samples_MetaphlanBySinglePart_thresh200/unicornMedianCoverageByGroup'
hierarchyPath="/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/krakenDBs/"+\
    "DB100MetaPhlanUnique/hierarchyDic.dat"
hierarchy=pandas.read_pickle(hierarchyPath)
toPlot=True
for sample in glob.glob(os.path.join(samplesPath,'*.dat')):
    unicornCoverageMedian={}
    unicornCoverageMedianWithZero={}
    unicornCoverageMean={}
    unicornCoverageMeanWithZero={}
    sampleName = sample.split('/')[-1].strip('.dat')
    obj = pandas.read_csv(sample)#pickle(sample)
    if 'count' in obj.columns:
        mappedToPart=obj.dropna(subset=['part'])
        groups = mappedToPart.groupby('group').sum()['count']
        groups=groups.sort_values(ascending=False)
    else:
        groups=obj.group.value_counts()

    print sample
    for i,grp in enumerate(groups.index):
        if groups[grp] == 0: #added by Sigal
            continue
        group=str(int(grp))
        if not os.path.exists(os.path.join(plotDir,group)):
            os.makedirs(os.path.join(plotDir,group))
        if not os.path.exists(os.path.join(plotDir, group,sampleName)):
            os.makedirs(os.path.join(plotDir, group,sampleName))

        objgroup=obj[obj.group==grp]
        if 'count' in obj.columns:
            genomes = objgroup.groupby('genome').sum()['count']
        else:
            genomes = objgroup.genome.value_counts()
        for idx,gnm in enumerate(genomes.index):
#        gnm=genomes.index[0]
            if genomes[gnm] == 0: # added by Sigal
                continue
            objgenome=objgroup[(objgroup.genome==gnm) & (~ np.isnan(objgroup.part))]
            if 'count' in objgenome:
                distibution=objgenome.set_index('result')['count']
            else:
                distibution=objgenome.part.value_counts()
            first = int(hierarchy[hierarchy.father==gnm].index.min())
            last=int(hierarchy[hierarchy.father==gnm].index.max())
            distibutionWithEmpty=pandas.Series(index=range(first,last+1))
            distibutionWithEmpty.loc[distibution.index]=distibution.values
            #distibutionWithEmpty.fillna(0,inplace=True)
            genome=str(int(gnm))
            distributionNoZero=distibution[distibution>0]
            medianCoverage=np.nanmedian(distributionNoZero.values)
            medianCoverageWithZero=np.nanmedian(distibution.values)
            meanCoverage = np.nanmean(distributionNoZero.values)
            meanCoverageWithZero = np.nanmean(distibution.values)
            unicornCoverageMedian[genome]=medianCoverage
            unicornCoverageMedianWithZero[genome]=medianCoverageWithZero
            unicornCoverageMean[genome]=meanCoverage
            unicornCoverageMeanWithZero[genome]=meanCoverageWithZero
            if toPlot:
                if i<5:
                    plt.figure()
                    plt.scatter(distibutionWithEmpty.index.values,1+distibutionWithEmpty.values,
                            label='window coverage')
                    plt.yscale('log',basey=10)
                    plt.xlabel('Genome position - unique window n=1000')
                    plt.ylabel('Count')
                    plt.ylim( 0, (max(distibutionWithEmpty.values)*2) ) #added by sigal
                    group_size=str(len(hierarchy[hierarchy.father==grp]))
                    relativeAbundance=pandas.read_csv(os.path.join(unicornCoverageByGroup,sampleName+'.csv')).set_index('group')
                    gRA=relativeAbundance.loc[int(group),'relative abundance']
                    sRA=medianCoverage/relativeAbundance.loc[int(group),'coverage']
                    plt.plot(distibutionWithEmpty.index.values,[medianCoverage+1]*len(distibutionWithEmpty.index.values),
                         'r',label='coverage %.2f RA %.2f'%(medianCoverage,sRA))
                    plt.legend(loc=3,bbox_to_anchor=(0.,0.925,1.,.102),ncol=2,borderaxespad=0.,mode='expand')
                    if 'result' in obj.columns and len(obj[obj['result']==gnm]['count'])>0:
                        plt.title("S:"+sampleName+'  G:'+genome+ '  group size:'+group_size+ \
                                  '   group RA: %.2f' % gRA)
                        #'  UpperLevel:'+str(int(obj[obj['result']==gnm]['count'])))
                    else:
                        plt.title("S:" + sampleName + '  G:' + genome + '  group size:' + group_size+\
                                  '   group RA: %.2f'%gRA)
                    plt.savefig(os.path.join(plotDir, group,sampleName,'%s_%sOf%s_%s.png'%(sampleName,idx,group_size,genome)))
                    plt.close('all')
    pandas.Series(unicornCoverageMedianWithZero).to_csv(os.path.join(unicornCoverageMedianWithZeroDir,sampleName+'.csv'))
    pandas.Series(unicornCoverageMedian).to_csv(os.path.join(unicornCoverageMedianDir,sampleName+'.csv'))
    pandas.Series(unicornCoverageMeanWithZero).to_csv(os.path.join(unicornCoverageMeanWithZeroDir, sampleName + '.csv'))
    pandas.Series(unicornCoverageMean).to_csv(os.path.join(unicornCoverageMeanDir, sampleName + '.csv'))
