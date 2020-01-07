import os
import pandas

base_path = "/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/datasets/list_and_phens_for_FD_samples"
opt = 2

#from PNPChip.PNPChip.ForPaper import extract_phenotypes
#extract_phenotypes.extract('include_allPNP','antropo','blood','covars').to_csv( os.path.join( base_path, 'Antropo_blood_pnp.csv'))
#extract_phenotypes.extract('s').to_csv(os.path.join( base_path, 'MetaPhlan_s_pnp.csv'))

withRegs = pandas.read_pickle('/home/daphnar/Microbiome/Analyses/AllSeqProjects/DFOut/StoolMetadataSpidDF.dat').reset_index()
withPhenotype=pandas.read_csv( os.path.join( base_path, 'Antropo_blood_pnp.csv'))
withRegs['RegistrantNumber']=withRegs['RegistrantNumber'].astype(int)
withRegs=withRegs.reset_index().merge(withPhenotype, left_on='RegistrantNumber',right_on='IID')
withRegs=withRegs[withRegs['PostSubSamp']==10e6]
withRegs['SPID']=withRegs['SPID'].astype(str).apply(lambda x: x[:-2] if x[-2:]=='.0' else x)
withRegs['StoolSample']=withRegs[['FD', 'SPID']].apply(lambda x: '_'.join(x), axis=1)
withRegs=withRegs[withRegs['StoolSample'].apply(lambda x: not (x[-1]=='_' or ',' in x))]
withRegs=withRegs[withRegs['IsBlacklisted'] !=1]
drop_cols = ['index','SPID','FD','RegistrantNumber','PostQC','RawReads',
             'HGMapped','PostHGF','Exp','PostSubSamp','GeneSetMap_Count',
             'GeneSetMap_Perc','SE/PE','StudyTypeID','IsBlacklisted',
             'ProductionDate','ConnectionID','dT','IID','Protain_g',
             'Calories_kcal','Fat_g','Carbs_g']

if opt == 1:
    # throw every participant that has more then one stool sample
    vals_to_throw = withRegs['RegistrantNumber'].value_counts()
    vals_to_throw = vals_to_throw[vals_to_throw>1].index
    print ("%d in merge" % len(withRegs))
    print ("throwing %d participants" % len(vals_to_throw))
    withRegs=withRegs[~withRegs['RegistrantNumber'].isin(vals_to_throw)]
    print ("%d after drop" % len(withRegs))
    withRegs=withRegs.drop(columns=drop_cols).set_index('StoolSample')
    withRegs.to_csv( os.path.join( base_path, 'SamplesAndPhenotypes_throw_participants_with_dup.csv'))
if opt == 2:
    # makeshift solution for keeping oldest stool sample for each participant
    withRegs.sort_values(['RegistrantNumber','ProductionDate'],inplace=True)
    withRegs['Reg_min1']=[-1] + list(withRegs.RegistrantNumber.values[:-1])
    withRegs=withRegs[withRegs.RegistrantNumber!=withRegs.Reg_min1]
    withRegs=withRegs.drop(columns=drop_cols+['Reg_min1']).set_index('StoolSample')
    withRegs.to_csv( os.path.join( base_path, 'SamplesAndPhenotypes_throw_dups_of_participants.csv'))
if opt == 3:
    # keep all samples
    withRegs=withRegs.drop(columns=drop_cols).set_index('StoolSample')
    withRegs.to_csv( os.path.join( base_path, 'SamplesAndPhenotypes_all.csv'))
if opt == 4:
    withRegs = withRegs[withRegs.StudyTypeID==1]
    withRegs = withRegs[withRegs.IsGenotec == 1]
    withRegs.sort_values(['RegistrantNumber','ProductionDate'],inplace=True)
    withRegs['Reg_min1']=[-1] + list(withRegs.RegistrantNumber.values[:-1])
    withRegs=withRegs[withRegs.RegistrantNumber!=withRegs.Reg_min1]
    withRegs=withRegs.drop(columns=drop_cols+['Reg_min1']).set_index('StoolSample')
    withRegs.to_csv( os.path.join( base_path, 'SamplesAndPhenotypes_Genotek_study1_throw_dups_of_participants.csv'))
if opt == 5:
    withRegs = withRegs[withRegs.StudyTypeID==1]
    withRegs = withRegs[withRegs.IsGenotec == 0]
    withRegs.sort_values(['RegistrantNumber','ProductionDate'],inplace=True)
    withRegs['Reg_min1']=[-1] + list(withRegs.RegistrantNumber.values[:-1])
    withRegs=withRegs[withRegs.RegistrantNumber!=withRegs.Reg_min1]
    withRegs=withRegs.drop(columns=drop_cols+['Reg_min1']).set_index('StoolSample')
    withRegs.to_csv( os.path.join( base_path, 'SamplesAndPhenotypes_Swab_study1_throw_dups_of_participants.csv'))

print ("Outputted %d phens for %d samples" % ( len(withRegs.columns), len(withRegs)))