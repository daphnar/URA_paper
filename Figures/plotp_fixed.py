import numpy as np
import scipy as sp
import scipy.stats as st
import scipy.interpolate
import scipy.linalg as la
import pylab as pl
from fastlmm.association.tests import Cv
import fastlmm.util.util as ut
import matplotlib.pyplot as plt
import seaborn
seaborn.set_style("white")
def pnames():
    '''
    The order here is important. The nearer the start of the list takes precedence if multiple ones
    are present.
    '''
    return ['P-value_adjusted',"P-value","pValue","p","P-value(50/50)","Pvalue"]

def rownames():
    return ["SetId","SNPId","setId","SNP","SnpId"]

def excludefiles():
    return ['work_directory','info','lrtperm','synthPhen.txt']

def threshlist():
    return sp.array([1e-10, 1e-9, 1e-8,1e-7,1e-6, 1e-5,1e-4,1e-3,1e-2,1e-1])


def getfiles(dirin, filepattern):
    import glob   
    filepatternorig=filepattern 
    filepattern=dirin + r"/" + filepattern    
    myfilestmp=glob.glob(filepattern)        
    myfiles=[]
    f=-1
    for f in myfilestmp: 
        keep=True
        for st in excludefiles():
            if st in f: keep=False
        if keep: myfiles.append(f)
      
    print str(len(myfiles)) + " files found"
    
    assert f>0, "no files found"
    return myfiles

def recalibrate(dirin,filepattern='*.txt', lrtpermfile=None, pnames=["P-value(50/50)"], rownames=rownames(), nullfit="qq",qmax=0.1, postfix="RECALIBRATED"):
    '''
    Read in each results file, use the null stats in lrtpermfile to re-calibrate the null distribution, 
    and then add a column to each file
    '''
    assert lrtpermfile is not None, "must provide lrtpermfile (output by FastLmmSet.py with lrt)"

    myfiles = getfiles(dirin, filepattern)
        
    ii=0
    for f in myfiles:       
        ii=ii+1
        print str(ii) + ") " + f
        pv,rowids,llnull,llalt = extractpvals(f,pnames,rownames)
        lrt = -2*(llnull-llalt)  
        alteqnull = (lrt==0)
        pv_adj = Cv.pv_adj_and_ind(nperm=0, pv_adj=None, nullfit=nullfit, lrt=lrt, lrtperm=None,
                                 alteqnull=alteqnull, alteqnullperm=None, qmax=qmax, 
                                 nullfitfile=lrtpermfile, nlocalperm=0, sort=False)[0]
        outfile = ut.appendtofilename(f,postfix)        
        np.savetxt(outfile, pv_adj)
            

def _qqplot_bar(M=1000000, alphalevel = 0.05,distr = 'log10'):
    '''
    calculate error bars for a QQ-plot
    --------------------------------------------------------------------
    Input:
    -------------   ----------------------------------------------------
    M               number of points to compute error bars
    alphalevel      significance level for the error bars (default 0.05)
    distr           space in which the error bars are implemented
                    Note only log10 is implemented (default 'log10')
    --------------------------------------------------------------------
    Returns:
    -------------   ----------------------------------------------------
    betaUp          upper error bars
    betaDown        lower error bars
    theoreticalPvals    theoretical P-values under uniform
    --------------------------------------------------------------------
    '''


    #assumes 'log10'

    mRange=10**(sp.arange(sp.log10(0.5),sp.log10(M-0.5)+0.1,0.1));#should be exp or 10**?
    numPts=len(mRange);
    betaalphaLevel=sp.zeros(numPts);#down in the plot
    betaOneMinusalphaLevel=sp.zeros(numPts);#up in the plot
    betaInvHalf=sp.zeros(numPts);
    for n in xrange(numPts):
        m=mRange[n]; #numplessThanThresh=m;
        betaInvHalf[n]=st.beta.ppf(0.5,m,M-m);
        betaalphaLevel[n]=st.beta.ppf(alphalevel,m,M-m);
        betaOneMinusalphaLevel[n]=st.beta.ppf(1-alphalevel,m,M-m);
        pass
    betaDown=betaInvHalf-betaalphaLevel;
    betaUp=betaOneMinusalphaLevel-betaInvHalf;

    theoreticalPvals=mRange/M;
    return betaUp, betaDown, theoreticalPvals

def pairedpvalsin_onefile(dirin,filepattern='*.txt',pnames=['P-value(50/50)','P-value_adjusted'], rownames=rownames(),logspace=True):
    '''
    Used to compare two sets of p-values in one file (e.g. 50/50 to adjusted).
    '''
    import os.path
    fs=10

    import glob
    filepattern=dirin + r'/' + filepattern    
    myfilestmp=glob.glob(filepattern)  
    myfiles=[]    
    for f in myfilestmp: 
        keep=True
        for str in excludefiles():
            if str in f: keep=False
        if keep: myfiles.append(f)
    indrange=sp.arange(0,len(myfiles))
    pv0={}
    pv1={}
    rowids={}
    label={}
    for j in indrange:      
        pv0[j],rowids[j],garb1,garb2=extractpvals(myfiles[j],[pnames[0]],rownames,sort=True)        
        if logspace: pv0[j]=-sp.log10(pv0[j])              
        label[j]=os.path.basename(myfiles[j])        
        
        pv1[j],rowids[j],garb1,garb2=extractpvals(myfiles[j],[pnames[1]],rownames,sort=True)        
        if logspace: pv1[j]=-sp.log10(pv1[j])
        
    import pylab as pl
    pl.ion()
    for j1 in indrange:
        assert len(pv0[j1]) == len(pv1[j1]), "different # of pvals"           
        pl.figure()
        pl.plot(pv0[j1],pv1[j1],'.')
        pl.xlabel(pnames[0])
        pl.ylabel(pnames[1])
        maxval=max(pl.xlim()[1],pl.ylim()[1])
        pl.plot([0,maxval],[0,maxval],'r--',linewidth=1)
        pl.title(label[j],fontsize=fs)
        fix_axes()

def pairedpvalsplot(dirin,filepattern='*.txt',pnames=pnames(), rownames=rownames(),minpval=1e-30,newplot=True,heatmap=False,extent=None,ms=5,fs=8,fliporder=False):
    '''
    Paired p-value plot in -log_10 space among all pairs of files in the directory that match the filepattern.
    Reorders all p-values so that rownames match up between files.
    Returns list of p-value lists, one per file
    '''
    import os.path
    import glob
    filepattern=dirin + r'/' + filepattern    
    myfilestmp=glob.glob(filepattern)  
    myfiles=[]       
    for f in myfilestmp: 
        keep=True
        for strn in excludefiles():
            if strn in f: keep=False
        if keep: 
            myfiles.append(f)
            print "keeping: " + f
        else:
            print "not including: " + f
    indrange=sp.arange(0,len(myfiles))
    if len(myfiles)==0: raise Exception("no files found")
    pv={};
    pvorig={}
    lambdas={}
    rowids={}
    label={}
    for j in indrange:      
        pv[j],rowids[j],dummy1,dummy2=extractpvals(myfiles[j],pnames,rownames,sort=True);
        pvorig[j]=pv[j]
        lambdas[j]=lambda_gc=estimate_lambda(pv[j])      
        #pv[j]=-sp.log10(pv[j])
        #pv[j][pv[j]<minpval]=minpval
        label[j]=os.path.basename(myfiles[j]) + ", $\mathrm{\mathsf{\lambda$=%1.3f}}" % lambdas[j]
        if j==0: 
            idorder=rowids[j]
        else:
            if not all(idorder==idorder):
                raise Exception("ids do not match up across file:" + label[j])

    import pylab as pl
    pl.ion()
    #these loops are to find out which are "notokany"
    notokany = sp.zeros_like(pv[0],dtype = bool)
    for j1 in indrange:      
        for j2 in sp.arange(j1+1,len(myfiles)):  
                assert len(pv[j1]) == len(pv[j2]), "different # of pvals in each file"                
                imag = (pv[j1]<=0.0) | (pv[j2]<=0.0)
                one = (pv[j1]>1.0) | (pv[j2]>1.0)
                iok = (~imag) & (~one)
                notokany = notokany | (~iok)                          
                
    #these are the standard plotting loops
    for j1 in indrange:      
        for j2 in sp.arange(j1+1,len(myfiles)):
                if fliporder: 
                    tmp=j1; j1=j2; j2=tmp;
                assert len(pv[j1]) == len(pv[j2]), "different # of pvals in each file"                
                if newplot: pl.figure()
                print "%i, %i" % (j1,j2)
                imag1 =(pv[j1]<=0.0)
                imag2 = (pv[j2]<=0.0)
                imag = imag1 | imag2
                one = (pv[j1]>1.0) | (pv[j2]>1.0)
                
                #print pv[j1][imag]
                #print pv[j2][imag]
                iok = (~imag) & (~one)       
                
                if not heatmap:                             
                    tmpj1=sp.copy(pv[j1])
                    tmpj2=sp.copy(pv[j2])
                    minpval1=min(tmpj1[~imag1])
                    minpval2=min(tmpj2[~imag2])
                    tmpj1[imag1]=minpval1
                    tmpj2[imag2]=minpval2
                                                            
                    pl.plot(-sp.log10(tmpj1[iok]),-sp.log10(tmpj2[iok]),'.k',markersize=ms)
                    #pl.plot(-sp.log10(tmpj1[notokany]),-sp.log10(tmpj2[notokany]),'b.',markersize=ms)
                    #maxval=max(pl.xlim()[1],pl.ylim()[1])
                    pl.plot(-sp.log10(tmpj1[imag1]),-sp.log10(tmpj2[imag1]),'g.',markersize=ms)
                    pl.plot(-sp.log10(tmpj1[imag2]),-sp.log10(tmpj2[imag2]),'r.',markersize=ms)

                    #pl.plot(-sp.log10(pv[j1][~iok]),-sp.log10(pv[j2][~iok]),'g.')
                    maxval=max(pl.xlim()[1],pl.ylim()[1])
                    pl.plot([0,maxval+1],[0,maxval+1],'b--',linewidth=1)                                  
                    #fix_axes()                     
                    pl.xlim([0,maxval+1])
                    pl.ylim([0,maxval+1])
                    pl.xlabel(label[j1],fontsize=fs)
                    pl.ylabel(label[j2],fontsize=fs)   
                else:     
                    heatmap, xedges, yedges = np.histogram2d(-sp.log10(pv[j1]),-sp.log10(pv[j2]), bins=100)
                    if extent is None: extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]                                   #heatmap=heatmap[-1:heatmap.shape[0]:-1,:]
                    pl.imshow(sp.log(heatmap+1), extent=extent,cmap=pl.cm.Greys)                   
                    pl.colorbar() 
                    import ipdb; ipdb.set_trace()         
    return pvorig                  
              
                
               

def type1errdir_agg(dirin,filepattern=['*.txt'],savefiles=False, pnames=pnames(), rownames=rownames(), alphalevel = 0.05,legend=None,xlim=None,ylim=None,plotsize="436x355",plot=False,dopower=True,auditstring=".S{0}.",auditrange=None,threshInd=threshlist()):
    '''
    filepattern is a list of filepatterns, each of which is processed seperate as in type1errdir
    '''
    assert dopower, "need to add in clause if false"
    S=len(filepattern)
    obs_count=[None for i in range(S)]
    obs_count_tot=sp.zeros([len(threshlist())])
    N_tot=0
    Norig=None        
    for s in sp.arange(0,S):
        if filepattern[s] is not None:
            thresh, obs_error, obs_count[s], p_twot, N, nfiles, missingind, duplicates=type1errdir(dirin, filepattern[s],dopower=True,auditrange=range(1,6),verbose=False)
            obs_count_tot = obs_count_tot + obs_count[s]
            N_tot+=N        
            if Norig is None:
                Norig=N
            else:
                if not (N==Norig):
                    type1errdir(dirin, filepattern[s],dopower=True,auditrange=range(1,6),verbose=True);
                    raise Exception("number of tests has changed, look for missing files")

    thresh=thresh[threshInd]
    obs_count_tot=obs_count_tot[threshInd]
    print "\n\n**********AGGREGATE************************"
    print "filepattern="
    for f in filepattern: print f
    padl=15
    strfm="%1.2e"

    print "Threshold".ljust(padl," ") + "\t" + "Obs. Count".ljust(padl," ") 
    print "--------------------------------------------------------"    
    for t in sp.arange(0,len(thresh)):
        print str(strfm % thresh[t]).ljust(padl," ") + "\t" + (str(int(obs_count_tot[t])) + "/" + str(N_tot)).ljust(padl," ")

    return thresh,obs_count_tot,N_tot

def type1errdir(dirin,filepattern='*.txt',savefiles=False, pnames=pnames(), rownames=rownames(), alphalevel = 0.05,legend=None,xlim=None,ylim=None,plotsize="436x355",plot=False,dopower=False,auditstring=".S{0}.",auditrange=None,verbose=True,thresh=threshlist()):
    '''
    Aggregate for type I error, and print out the observed type I error
    filepattern can use python regular expressions. For e.g. "wtcbf*P[0-2]*.txt" allows any digit from 0 to 2 after P
    Use filepattern (reg expression) and auditrange to make sure you get all the values you expect, will print the missing ones to screen
    auditrange=range(1,501) with filepattern=".S*." will make sure that each file contains SX, where X goes from 1-500    
    Output: filepattern, thresh, obs_error, obs_count, pv, nfiles, missingind, duplicates
    '''
    import re
    import glob    
    filepatternorig=filepattern
    filepattern=dirin + r"/" + filepattern    
    myfilestmp=glob.glob(filepattern)        
    myfiles=[]    
    for f in myfilestmp: 
        keep=True
        for st in excludefiles():
            if st in f: keep=False
        if keep: myfiles.append(f)
  
    if verbose: print str(len(myfiles)) + " files found"    
    ycoord=0    
  
    auditdic={}; duplicates=[]; extras=[]; auditind=None
    
    #put all p-values in to one qqplot, and also report the type 1 errors 
    #for 1e-1 through to smallest order of magnitude
    pv=[]    
    obs_count=-1*sp.zeros([len(thresh)])        
    obs_error=-1*sp.zeros([len(thresh)])  
    p_onet=sp.NaN*sp.zeros([len(thresh)])     
    p_twot=sp.NaN*sp.zeros([len(thresh)])     
    ii=0        
    for f in myfiles:        
        if "info.txt" not in f:                 
            import os
            shortf=os.path.split(f)[1]   
            #if verbose: print str(ii) + ") " + shortf            
            if auditrange is not None: 
                searchstring=str.format(auditstring,"([0-9]*)") #".S([0-9]*)."            
                auditind= float(re.search(searchstring,shortf).groups()[0])
            if auditdic.has_key(auditind): 
                if auditrange is not None: duplicates.append(int(auditind))
            elif auditrange is not None and auditind not in auditrange: 
                extras.append(int(auditind))
            else:
                if auditrange is not None: auditdic[auditind]=True
                pvtmp,rowids,garb1,garb2=extractpvals(f,pnames,rownames,includefilename=False,verbose=False)            
                pvtmp.sort()
                pv.extend(pvtmp)
                # tally up type I error
                ct=0   
                t=0 
                j=0;           
                while True:
                    if j>=len(pvtmp): break
                    if pvtmp[j]<=thresh[t]:
                        ct=ct+1
                        j=j+1
                    else: 
                        obs_count[t]=obs_count[t]+ct
                        t=t+1
                        if t>=len(thresh): break                                
                ii=ii+1   
    N=len(pv) 
    
    for t in sp.arange(0,len(thresh)):
        p_onet[t]=1 - sp.stats.binom.sf(obs_count[t]-1, N, thresh[t]) # -1 is there to include 51 as well ;-) (one-tailed test)
        p_twot[t]=sp.stats.binom_test(obs_count[t], N, thresh[t])# (two-tailed test)
        
    nfiles=ii
    missingind=[]    
    if verbose:        
        if auditrange is not None:        
            nummis=0
            msg="missinginds=["
            for a in auditrange:        
                if not auditdic.has_key(a):
                    nummis+=1
                    msg+=str(a) + ", "
            msg+="]"
            print "--------------------------------------------------------"
            print msg
            print "found " + str(nummis) + " missing entries"
        print "--------------------------------------------------------"
        if len(duplicates)>0:
            msg="["
            print "ignored duplicates for ["
            for a in duplicates:
                msg+= str(a) + ","
            msg+= "]"
            print msg
        else: print "found no duplicates"
        print "--------------------------------------------------------"
        if len(extras)>0:
            msg="["
            print "ignored extras ["
            for a in extras:
                msg+= str(a) + ","
            msg+= "]"
            print msg
        else: print "found no extras"
        print "--------------------------------------------------------"
                
        print "--------------------------------------------------------"
        print "found total of " + str(N) + " p-values over " + str(nfiles) + " files with pattern:\n " + filepatternorig
        print "from dir= " + dirin
        print "--------------------------------------------------------"
        padl=15
    
        if not dopower:
            print "Threshold".ljust(padl," ") + "\t" + "Obs. Error".ljust(padl," ") + "\t" + "Obs. Count".ljust(padl," ")  +"\t" + "two-tail P".ljust(padl," ") #+ "\t" + "one-tail P".ljust(padl," ")
            print "--------------------------------------------------------"
            strfm="%1.2e"
            for t in sp.arange(0,len(thresh)):
                print str(strfm % thresh[t]).ljust(padl," ") + "\t" + str(strfm % (obs_count[t]/float(N))).ljust(padl," ") + "\t" + (str(int(obs_count[t])) + "/" + str(N)).ljust(padl," ")  + "\t" + str(strfm % (p_twot[t])).ljust(padl," ")
        else:
            print "Threshold".ljust(padl," ") + "\t" + "Obs. Count".ljust(padl," ") 
            print "--------------------------------------------------------"
            strfm="%1.2e"
            for t in sp.arange(0,len(thresh)):
                print str(strfm % thresh[t]).ljust(padl," ") + "\t" + (str(int(obs_count[t])) + "/" + str(N)).ljust(padl," ")
    if plot:
        fileout=None
        qqplot(pv, fileout, alphalevel,legend,xlim,ylim)     
        
    return thresh, obs_error, obs_count, p_twot, N, nfiles, missingind, duplicates



def qqplotdir(dirin,filepattern='*.txt',savefiles=False, pnames=pnames(), rownames=rownames(), alphalevel = 0.05,legend=None,xlim=None,ylim=None,plotsize="436x355",aggregate=False,dohist=True):
    '''
    Qqplot and histogram of p-values in pname. Does not print any filename containing in "info.txt"
    Currently only saves the qqplots to file, not the histograms.
    filepattern can use python regular expressions. For e.g. "wtcbf*P[0-2]*.txt" allows any digit from 0 to 2 after P
    aggregate=True puts all the p-values in one plot
    '''
    import glob   
    filepatternorig=filepattern 
    filepattern=dirin + r"/" + filepattern    
    myfilestmp=glob.glob(filepattern)        
    myfiles=[]
    f=-1
    for f in myfilestmp: 
        keep=True
        for st in excludefiles():
            if st in f: keep=False
        if keep: myfiles.append(f)
      
    print str(len(myfiles)) + " files found"
    maxycoord=800
    ycoordshift=int(plotsize.split("x")[1])
    ycoord=-ycoordshift
    ii=0 
    
    allp=[]   

    assert f>0, "no files found"
        
    for f in myfiles:       
        ii=ii+1
        print str(ii) + ") " + f
        ycoord=ycoord+ycoordshift
        if ycoord>maxycoord: ycoord=ycoord-ycoordshift+10   
        if savefiles: fileout = f[0:-4]+".out.png" 
        else: fileout=None   
        if (not aggregate):  
            qqplotfile(f, fileout, pnames,rownames, alphalevel,legend,xlim,ylim,str(ycoord),plotsize,dohist=dohist)        
        else:
            pv,rowids,garb1,garb2=extractpvals(f,pnames,rownames)
            #allp.append(pv) #this gives all seperate plots on one (makes several series)
            allp=allp+pv.tolist() #this just makes one series
    if aggregate:
        title="agg over " + filepatternorig + " (" + str(len(allp)) + " p-values)"
        qqplotp(sp.array(allp),fileout = None, pnames=pnames, rownames=rownames,alphalevel =alphalevel,legend=None,title=title,dohist=dohist)
   

def extractpvals(filein,pnames=pnames(), rownames=rownames(),sort=False,includefilename=True,verbose=True,altname="LogLikeAlt",nullname="LogLikeNull"):
    import pandas as pd          
    header = pd.read_csv(filein,delimiter = '\t',usecols=[])
    pname=None
    rowname=None
    for p in pnames:
       if p in header.keys():
           pname=p
           break
    for r in rownames:
       if r in header.keys():
           rowname=r
           break   
    if pname is None: 
       print "header: " +  str(header.keys())
       raise Exception(str(pnames) + " not found")    
    if rowname is None:
       print "header: " +  str(header.keys())
       raise Exception(str(rownames) + " not found")
    try:
        data=pd.read_csv(filein,delimiter = '\t',dtype={pname:np.float64,nullname:np.float64,altname:np.float64},usecols=[pname,rowname,nullname,altname])
    except:
        data=pd.read_csv(filein,delimiter = '\t',dtype={pname: np.float64},usecols=[pname,rowname])
    if sort:        
        data=data.sort([rowname],ascending=[True])
    pv=data[pname].values      
    try:
        llalt = data[altname]
        llnull  =data[nullname]
    except:
        llalt = None
        llnull= None
    rowids=data[rowname].values
    import os    
    if verbose: 
        print "M=" + str(len(pv)) + ", pname=" + pname + ", rowname=" + rowname 
        if includefilename: print " (" + os.path.split(filein)[1] + ")"
    return pv,rowids,llnull,llalt

def qqplotfile(filein,fileout = None, pnames=pnames(), rownames=rownames(),alphalevel = 0.05,legend=None,xlim=None,ylim=None,ycoord=10,plotsize="652x526",dohist=True):
     pv,rowids,garb1,garb2=extractpvals(filein,pnames,rownames)
     import os.path
     title=os.path.basename(filein)
     qqplotp(pv,fileout, pnames, rownames,alphalevel,legend,xlim,ylim,ycoord,plotsize,title=title,dohist=dohist)

def qqplotp(pv,fileout = None, pnames=pnames(), rownames=rownames(),alphalevel = 0.05,legend=None,xlim=None,ylim=None,ycoord=10,plotsize="652x526",title=None,dohist=True):
     '''
     Read in p-values from filein and make a qqplot adn histogram.
     If fileout is provided, saves the qqplot only at present.
     Searches through p until one is found.   '''       
     
     import pylab as pl     
     pl.ion()     
          
     fs=8     
     h1=qqplot(pv, fileout, alphalevel,legend,xlim,ylim,addlambda=True)
     #lambda_gc=estimate_lambda(pv)
     #pl.legend(["gc="+ '%1.3f' % lambda_gc],loc=2)     
     pl.title(title,fontsize=fs)
     #wm=pl.get_current_fig_manager()
     #e.g. "652x526+100+10
     xcoord=100
     #wm.window.wm_geometry(plotsize + "+" + str(xcoord) + "+" + str(ycoord))

     if dohist:
         h2=pvalhist(pv)
         pl.title(title,fontsize=fs)
         #wm=pl.get_current_fig_manager()
         width_height=plotsize.split("x")
         buffer=10
         xcoord=int(xcoord + float(width_height[0])+buffer)
         #wm.window.wm_geometry(plotsize + "+" + str(xcoord) + "+" + str(ycoord))
     else: h2=None

     return h1,h2
     
def pvalhist(pv,numbins=50,linewidth=3.0,linespec='--r'):    
    '''
    Plots normalized histogram, plus theoretical null-only line.
    '''
    import pylab as pl   
    h2=pl.figure()  
    [nn,bins,patches]=pl.hist(pv,numbins,normed=1)    
    pl.plot([0, 1],[1,1],linespec,linewidth=linewidth)
        
def tmp():

    import scipy as sp
    import fastlmm.util.stats.plotp as pt    

    pall=[]
    for j in sp.arange(0,100):
        p=sp.rand(1,1000)
        p100=sp.repeat(p,100)
        pall.append(p100)
    qqplotavg(pall)


def qqplotavg(pvallist, fileout = None, alphalevel = 0.05,legend=None,xlim=None,ylim=None,fixaxes=True):
    '''
    Average the qqplots in pvallist (assumes qq theoretical is same for all, and averages qq empirical)
    '''    
    L=len(pvallist)
    qemp=None
    for i in sp.arange(0,L):
        pval =pvallist[i].flatten()
        M = pval.shape[0]
        pnull = (0.5 + sp.arange(M))/M
        qnull = -sp.log10(pnull)
        qemptmp = -sp.log10(sp.sort(pval))   
        if qemp is None:
            qemp=qemptmp
        else:
            qemp=qemp+qemptmp
    qavg=qemp/L               
    qqplotfromq(qnull,qavg)    
    
def qqplotfromq(qnull,qemp):
    '''
    Given uniform quartile values, and emprirical ones, make a qq plot
    '''    
    pl.ion()
    pl.plot(qnull, qemp, '.',markersize = 2)                
    addqqplotinfo(qnull,qnull.flatten().shape[0])
    
def addqqplotinfo(qnull,M,xl='P value observed (-log10)',yl='P value expected (-log10)',xlim=None,ylim=None,alphalevel=0.05,legendlist=None,fixaxes=False, fontsize=12):    
    distr='log10'
    if ylim is None:
        pl.plot([0,qnull.max()], [0,qnull.max()],'k')
    else:
        pl.plot([0,ylim[1]], [0,ylim[1]],'k')
    pl.ylabel(xl, fontsize=fontsize)
    pl.xlabel(yl, fontsize=fontsize)
    if xlim is not None:
        pl.xlim(xlim)
    if ylim is not None:
        pl.ylim(ylim)
    if alphalevel is not None:
        if distr == 'log10':
            betaUp, betaDown, theoreticalPvals = _qqplot_bar(M=M,alphalevel=alphalevel,distr=distr)
#             lower = -sp.log10(theoreticalPvals-betaDown)
#             upper = -sp.log10(theoreticalPvals+betaUp)
#             pl.fill_between(-sp.log10(theoreticalPvals),lower,upper,color="grey",alpha=0.5)
            #pl.plot(-sp.log10(theoreticalPvals),lower,'g-.')
            #pl.plot(-sp.log10(theoreticalPvals),upper,'g-.')
    if legendlist is not None:
#         leg = pl.legend(legendlist, loc=4, numpoints=1, frameon=False, fontsize=fontsize)
        plt.gca().annotate(legendlist[0], xy=(0.9, 0.25), xycoords='axes fraction', fontsize=fontsize,
                xytext=(0, -15), textcoords='offset points',
                ha='right', va='top')
        # set the markersize for the legend
#         for lo in leg.legendHandles:
#             lo.set_markersize(0)
        
#         plt.legend(markerscale=0)
        plt.gca().grid(False)
    if fixaxes:
        fix_axes()        


def qqplot(pvals, fileout = None, alphalevel = 0.05,
           legend=None,xlim=None,ylim=None,fixaxes=True,
           addlambda=True,minpval=1e-20,title=None,h1=None,
           figsize=[5,5],grid=True,fontsize=12,markersize=2,
           color=None,group1=None,group2=None):
    '''
    performs a P-value QQ-plot in -log10(P-value) space
    -----------------------------------------------------------------------
    Args:
        pvals       P-values, for multiple methods this should be a list (each element will be flattened)
        fileout    if specified, the plot will be saved to the file (optional)
        alphalevel  significance level for the error bars (default 0.05)
                    if None: no error bars are plotted
        legend      legend string. For multiple methods this should be a list
        xlim        X-axis limits for the QQ-plot (unit: -log10)
        ylim        Y-axis limits for the QQ-plot (unit: -log10)
        fixaxes    Makes xlim=0, and ylim=max of the two ylimits, so that plot is square
        addlambda   Compute and add genomic control to the plot, bool
        title       plot title, string (default: empty)
        h1          figure handle (default None)
        figsize     size of the figure. (default: [5,5])
        grid        boolean: use a grid? (default: True)
    Returns:   fighandle, qnull, qemp
    -----------------------------------------------------------------------
    '''    
    distr = 'log10'
    import pylab as pl
    if type(pvals)==list:
        pvallist=pvals
    else:
        pvallist = [pvals]
    if type(legend)==list:
        legendlist=legend
    elif legend is not None:
        legendlist = [legend]
    else: legendlist = None
    
    if h1 is None:
        h1=pl.figure(figsize=figsize)
         
    maxval = 0

    for i in xrange(len(pvallist)):        
        pval =pvallist[i].flatten()
        M = pval.shape[0]
        pnull = (0.5 + sp.arange(M))/M
        # pnull = np.sort(np.random.uniform(size = tests))
                
        pval[pval<minpval]=minpval
        pval[pval>=1]=1

        if distr == 'chi2':
            qnull = st.chi2.isf(pnull, 1)
            qemp = (st.chi2.isf(sp.sort(pval),1))
            xl = 'LOD scores'
            yl = '$\chi^2$ quantiles'

        if distr == 'log10':
            qnull = -sp.log10(pnull)            
            qemp = -sp.log10(sp.sort(pval)) #sorts the object, returns nothing
            xl = 'Observed P value (-log10)'
            yl = 'Expected P value (-log10)'
        if not (sp.isreal(qemp)).all(): raise Exception("imaginary qemp found")
        # if qnull.max>maxval:
        #     maxval = qnull.max()


        # plt.axhline(7.5,linestyle="--",color='gray')
        if color is None:
            color=np.array([69.0,117.0,180.0])/255.0
        if group1 is not None and group2 is not None:
            pl.plot(qnull[group1], qemp[group1], '+', markersize=markersize,c=color[0])
            pl.plot(qnull[group2], qemp[group2], '_', markersize=markersize,c=color[1])
        else:
            pl.plot(qnull, qemp, 'o', markersize=markersize,color=color)
        #pl.plot([0,qemp.max()], [0,qemp.max()],'r')        
        if addlambda:
            lambda_gc = estimate_lambda(pval)
            print "lambda=%1.4f" % lambda_gc
            #pl.legend(["gc="+ '%1.3f' % lambda_gc],loc=2)   
            # if there's only one method, just print the lambda
            if len(pvallist) == 1:
                legendlist=["$\mathrm{\mathsf{\lambda_{GC}=}}$%1.4f" % lambda_gc]   
            # otherwise add it at the end of the name
            else:
                legendlist[i] = legendlist[i] + " ($\mathrm{\mathsf{\lambda_{GC}}}=$%1.4f)" % lambda_gc

    addqqplotinfo(qnull,M,xl,yl,xlim,ylim,alphalevel,legendlist,fixaxes, fontsize=fontsize)
    
    if title is not None:
        pl.title(title)            
    
    if fileout is not None:
        pl.savefig(fileout)

    return h1,qnull, qemp,

def fix_axes(buffer=0.1):
    '''
    Makes x and y max the same, and the lower limits 0.
    '''    
    maxlim=max(pl.xlim()[1],pl.ylim()[1])    
    pl.xlim([0-buffer,maxlim+buffer])
    pl.ylim([0-buffer,maxlim+buffer])

def estimate_lambda(pv):
    '''
    estimate the lambda for a given array of P-values
    ------------------------------------------------------------------
    pv          numpy array containing the P-values
    ------------------------------------------------------------------
    L           lambda value
    ------------------------------------------------------------------
    '''
    LOD2 = np.median(st.chi2(1).isf(pv))
    L = (LOD2/0.456)
    return L








def _compute_x_positions_snps(positions, chromosome_starts):
    cumulative_pos = np.zeros(positions.shape[0])
    for i, chromosome_start in enumerate(chromosome_starts):
        idx_chr = positions[:,0]==chromosome_start[0]
        cumulative_pos[idx_chr] = positions[idx_chr][:,1] + chromosome_start[1]
    return cumulative_pos
def _rel_to_midpoint(rle):
    previous_count = 0
    for item, count in rle:
        yield previous_count + count // 2
        previous_count += count

def _color_list(chr_list,rle):
    chr_to_index = dict((chr,index) for index,(chr,count) in enumerate(rle))
    index_to_color = {1:np.array([69,117,180])/255.0,0:np.array([145,191,219])/255.0} #{1:"purple",0:"pink"}
    result = [index_to_color[chr_to_index[chr]%len(index_to_color)] for chr in chr_list]
    return result
def _compute_x_positions_chrom(positions, offset=1e5):
    chromosomes = np.unique(positions[:,0])
    chromosomes.sort()
    chromosome_starts = np.zeros((chromosomes.shape[0],3),dtype="object")
    chr_start_next = 0
    for i, chromosome in enumerate(chromosomes):
        pos_chr = positions[positions[:,0]==chromosome]
        chromosome_starts[i,0] = chromosome                     #the chromosome
        chromosome_starts[i,1] = chr_start_next                 #start of the chromosome
        chromosome_starts[i,2] = chr_start_next + pos_chr.max() #end of the chromosome
        chr_start_next = chromosome_starts[i,2] + offset
    return chromosome_starts
def _run_length_encode(seq):
    count = 0
    previous = None
    for item in seq:
        if count == 0:
            count = 1
            previous = item
        elif item == previous:
            count += 1
        else:
            yield previous, count
            previous = item
            count =1
    if count > 0:
        yield previous, count
def manhattan_plot(chr_pos_pvalue_array,pvalue_line=None,plot_threshold=1.0,vline_significant=False,
                   marker="o", chromosome_starts=None, xaxis_unit_bp=True, alpha=0.5,ticks=True,s=None,extraTicks=[]):
    """
    Function to create a Manhattan plot.  See http://en.wikipedia.org/wiki/Manhattan_plot.

    Args:
        chr_pos_pvalue_array:   an n x 3 numpy array. The three columns are the chrom number 
                                (as a number), the position, and pvalue.
                                :type chr_pos_pvalue_array: numpy array
        pvalue_line:            (Default: None). If given, draws a line at that PValue.
                                :type pvalue_line: a 'pheno dictionary' or a string
        plot_threshold:         plot only SNPs that achieve a P-value smaller than pvalue_threshold
                                to speed up plotting
        vline_significant:      boolean. Draw a vertical line at each significant Pvalue?
                                :rtype: none, but changes the global current figure.
        marker:                 marker for the scatter plot. default: "o"
        chromosome_starts:      [Nchrom x 3] ndarray: chromosome, cumulative start position, cumulative stop position
                                cumulative chromosome starts, for plotting. If None (default), this is estimated from data
        xaxis_unit_bp:          plot cumulative position in basepair units on x axis? If False, only 
                                use rank of SNP positions. (default: True)
        alpha:                  alpha (opaquness) for P-value markers in scatterplot (default 0.5)

    Returns:
        chromosome_starts       [Nchrom x 3] ndarray: chromosome, cumulative start position, cumulative stop position
                                cumulative chromosome starts used in plotting.

    :Example:

    >>> from fastlmm.association import single_snp
    >>> from pysnptools.snpreader import Bed
    >>> import matplotlib.pyplot as plt
    >>> import fastlmm.util.util as flutil
    >>> pheno_fn = "../feature_selection/examples/toydata.phe"
    >>> results_dataframe = single_snp(test_snps="../feature_selection/examples/toydata.5chrom", pheno=pheno_fn, h2=.2)
    >>> chromosome_starts = flutil.manhattan_plot(results_dataframe.as_matrix(["Chr", "ChrPos", "PValue"]),pvalue_line=1e-7)
    >>> #plt.show()

    """

    # create a copy of the data and sort it by chrom and then position
    array = np.array(chr_pos_pvalue_array)
    if plot_threshold:
        array = array[array[:,2]<=plot_threshold]
    else:
        plot_threshold = 1.0
    array=array[np.argsort(array[:,1]),:] #sort by ChrPos
    array=array[np.argsort(array[:,0],kind='mergesort'),:] #Finally, sort by Chr (but keep ChrPos in case of ties)
    rle = list(_run_length_encode(array[:,0]))
        
    if xaxis_unit_bp:   #compute and use cumulative basepair positions for x-axis
        if chromosome_starts is None:
            chromosome_starts = _compute_x_positions_chrom(array)
        chr_pos_list = _compute_x_positions_snps(array, chromosome_starts)
        plt.xlim([0,chromosome_starts[-1,2]+1])
        if ticks:
            plt.xticks(chromosome_starts[:,1:3].mean(1),chromosome_starts[:,0])
    else:               #use rank indices for x-axis
        chr_pos_list = np.arange(array.shape[0])
        xTickMarks = [str(int(item)) for item,count in rle]
        plt.xlim([0,array.shape[0]])
        if ticks:
            plt.xticks(list(_rel_to_midpoint(rle)), xTickMarks)
    y = -np.log10(array[:,2])
    max_y = y.max()

    if pvalue_line and vline_significant:   #mark significant associations (ones that pass the pvalue_line) by a red vertical line:
        idx_significant = array[:,2]<pvalue_line
        if np.any(idx_significant):
            y_significant = y[idx_significant]
            chr_pos_list_significant = chr_pos_list[idx_significant]
            for i in xrange(len(chr_pos_list_significant)):
                plt.axvline(x=chr_pos_list_significant[i],ymin = 0.0, ymax = y_significant[i], color = 'r',alpha=0.8)
    
    plt.scatter(chr_pos_list,y,marker=marker,c=_color_list(array[:,0],rle),edgecolor='none',s=(y/max_y*20+0.5 if s is None else s), alpha=alpha)
    plt.xlabel("chromosome")
    plt.ylabel("-log10(P value)")
    print "Max pvalue:",max(array[:,2])
    if pvalue_line:
        printTicks=['0','5']+[str(extra) for extra in extraTicks]
        plt.yticks([0,5]+extraTicks,printTicks,fontsize=12)
        for val in pvalue_line:
            plt.axhline(-np.log10(val),linestyle="--",color='gray')
#             ###############Add pvalue###############
# #             printTicks.append(r"$%0de^{-%02d}$"%(-np.log10(val),-np.log10(val)+1))
#             b=int(np.ceil(-np.log10(val)))
#             a = int(round(val*10**b))
#             #printTicks.append('$%0.1fe^{-%02d}$'%(-np.log10(val),-np.log10(val)+1))
#             printTicks.append('%.1f'%-np.log10(val))#'$\mathrm{\mathsf{%1d\cdot10^{-%02d}}}$'%(a, b))
#         plt.yticks([0,5]+\
#                    list(-np.log10(np.array(pvalue_line))),printTicks,fontsize=12)
    if len(extraTicks)==0:
        plt.ylim([0,None])#[-np.log10(plot_threshold),None]
    else:
        plt.ylim([0,max(extraTicks)])
    return chromosome_starts