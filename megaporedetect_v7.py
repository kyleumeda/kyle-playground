#Classification:
# 1 pore to bilayer
# 2 pore to short
# 3 pore/megapore to bilayer
# 4 pore/megapore to short
# 5 megapore to bilayer
# 6 megapore to short

# 9 pore lasted all run - never deactivated
# number > 10 number of cycle it was a Megapore(c) - never deactivated
# -2 shorted
# -1 is bilayer

#export-cycle-stats -a P_0_2015-09-01-153518_KB-server4_ac-analysis_v2.2.0/annotations.h5 --primary raw/multi_ubf.h5 --oc-calibration raw/oc_calibration.h5


import matplotlib;
#matplotlib.use("Agg")
import numpy as np, os,sys,pandas as pd, matplotlib.pyplot as plt, h5py,pdb
#import matplotlib;
#matplotlib.use("Agg")
matplotlib.verbose.set_level('debug') 
from matplotlib import animation
import datetime,glob
starttime = datetime.datetime.now().time()


def contiguous_regions(condition,N):
    try:
        idxT=np.where(np.diff(np.where(np.concatenate(([condition[0]],condition[:-1] != condition[1:],[True])))[0])[::2]>N)[0][0]
        idx= np.where(np.concatenate(([condition[0]],condition[:-1] != condition[1:],[True])))[0][2*idxT]
    except:
        idx=-1
    return idx


def mad_based_outlier(points, thresh=3.5):
    if len(points.shape) == 1:
        points = points[:,None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)

    modified_z_score = 0.6745 * diff / med_abs_deviation

    return modified_z_score > thresh

subsample=1
datapath=os.getcwd()
#datapath=sys.argv[1]
P_dir = glob.glob("P*")[0]
cell_anno = P_dir + '/cell_annotations.csv'
cell_anno_df = pd.DataFrame; cell_anno_df = cell_anno_df.from_csv(cell_anno)

subsample=int(subsample)
print "Subsample to every: " +str(subsample)
#pdb.set_trace()
#path='/Users/wahbaa/Data/deacphase/150922_SIG-D_01_vaporeon_WL21R13C07/raw/multi_ubf'
path=datapath + '/raw/multi_ubf'
#path='sweep_2'
f=h5py.File(path+'.h5')

#f=h5py.File('sweep_2.h5')

#checking for dark cycle reduction

if np.where((f['other_traces/averaged_frame_mask'].value)>0):
    dark_cycle_reduc = True
    print 'Dark Cycle Reduction enabled'
else:
    dark_cycle_reduc = False
    print 'Dark Cycle Reduction disabled'


#getting voltages and time traces
liq=f['other_traces']['liq'].value;liq=(liq/(2.0**16-1.))*(7./8.)*(2.048);
v_pre=f['other_traces']['v_pre'].value;v_pre=(v_pre/(2.0**16-1.))*(7./8.)*(2.048);
ts=f['other_traces/timestamp'].value;ts=ts-ts[0];ts=ts/1e6

vapp=(v_pre-liq)*1000.;vapp=vapp.astype(np.float64)

ts_orig=ts

tsjump=np.diff(ts)
tsjumplocation=np.where(tsjump>0.5)
print tsjumplocation[0]

print tsjumplocation[0].shape
starts=np.zeros(1+len(tsjumplocation[0]),'int')
stops=np.zeros(1+len(tsjumplocation[0]),'int')
starts[0]=0
for section in range(1,1+len(tsjumplocation[0]),1):
    stops[section-1]=tsjumplocation[0][section-1] 
    starts[section]=stops[section-1]+1
stops[-1]= ts.size  
print starts
print stops
print starts.size


posvol=np.where(vapp>0)[0] #positive regions
negvol=np.where(vapp<0)[0] #negative regions


#getting the positive regions
firstPosptslocs=np.where(np.diff(posvol)>1)[0]+1
firsttPospts=ts[posvol[firstPosptslocs]]
firstvappPospts=vapp[posvol[firstPosptslocs]]
#firstdataPospts=data[posvol[firstPosptslocs]]
#getting the last positive regions
lastPosptslocs=np.where(np.diff(posvol)>1)[0]
lasttPospts=ts[posvol[lastPosptslocs]]
lastvappPospts=vapp[posvol[lastPosptslocs]]
#firstdataPospts=data[posvol[firstPosptslocs]]


#getting the first negative regions
firstNegptslocs=np.where(np.diff(negvol)>1)[0]+1
firsttNegpts=ts[negvol[firstNegptslocs]]
firstvappNegpts=vapp[negvol[firstNegptslocs]]
#firstdataNegpts=data[negvol[firstNegptslocs]]

#getting the last negative regions
lastNegptslocs=np.where(np.diff(negvol)>1)[0]
lasttNegpts=ts[negvol[lastNegptslocs]]
lastvappNegpts=vapp[negvol[lastNegptslocs]]
#firstdataNegpts=data[negvol[firstNegptslocs]]


#pdb.set_trace()
# trying a new scheme to make the array sizes the same
firstptMinlen=np.min([len(firstPosptslocs),len(firstNegptslocs)])
firstPosptslocs=firstPosptslocs[:firstptMinlen]
firstNegptslocs=firstNegptslocs[:firstptMinlen]

#getting the last locs the same size
lastptMinlen=np.min([len(lastPosptslocs),len(lastNegptslocs)])
lastPosptslocs=lastPosptslocs[:lastptMinlen]
lastNegptslocs=lastNegptslocs[:lastptMinlen]

df=pd.DataFrame();
dfrmsbright=pd.DataFrame();
dfoutlierbright=pd.DataFrame()
dfrmsdark=pd.DataFrame();
dfoutlierdark=pd.DataFrame()
dffpd=pd.DataFrame();
dfvliq=pd.DataFrame();
dfactivetime=pd.DataFrame();
df=pd.DataFrame();


dflast=pd.DataFrame();

dfdeactivationtime=pd.DataFrame()
dfshorttime=pd.DataFrame()
dfdeactivationcause=pd.DataFrame()
dfcharge=pd.DataFrame()
dfMEGAPORE=pd.DataFrame()
dfPORE=pd.DataFrame()
dfSHORTED=pd.DataFrame()
dfBILAYER=pd.DataFrame()
dfDEACTIVATED=pd.DataFrame()
dfSINGLEPORE=pd.DataFrame()
dfMULTIPORE=pd.DataFrame()
dfEnthoughtSP=pd.DataFrame()
dfEnthoughtSeq=pd.DataFrame()

#df['cells']=f['cells'].keys()


cellcount=0
#celllist=[cellnameP,cellnamePtobilayer,cellnameMP,cellnameMPtoshort,cellnameMPtobilayer,cellnamePMP,cellnamePMPtobilayer,cellnamecrap]
#celllist=['b01c7019', 'b05c6111','b09c7698','b00c4099','b08c7982','b14c0885','b09c1445','b09c1442','b10c3476']
#celllist=['b00c4099']
#celllist=['b01c7019']
#celllist=['b09c1445']
#celllist=['b10c3476']
#celllist=['b09c7698']
#celllist=['b13c7185']
#celllist=['b01c4842']
#celllist=['b09c4921']
#celllist=['b11c5814']
#celllist=['b12c5512']
#for cellOI in celllist:

for cellOI in f['cells'].keys():
    if cellcount>200000:
        break
    cellcount=cellcount+1
    #outliersbright=[]
    #outliersdark=[]
    #rmsbrightlist=[]
    #rmsdarklist=[]
    vliqlist=[]
    fpdlist=[]
    activetimelist=[]
    deactivationtimelist=[]
    shorttimelist=[]
    deactivationcauselist=[]
    chargelist=[]
    MEGAPORE = []
    PORE = []
    SINGLEPORE = []
    MULTIPORE = []
    DEACTIVATED = []
    BILAYER = []
    SHORTED = []
    
    # get first point difference
    data=f['cells'][cellOI].value
    data=data.astype(np.float)

    # loop over sections
#    for section in range(50,65,1):
    #for section in range(len(starts)):
    for section in range(1,2,1):    
        sectiondata=data[starts[section]:stops[section]]
        sectionvapp=vapp[starts[section]:stops[section]]
        sectionts=ts[starts[section]:stops[section]]
        sectiont0=sectionts[0]
        
        posvol=np.where(sectionvapp>0)[0] #positive regions
        negvol=np.where(sectionvapp<0)[0] #negative regions

        bright=sectiondata[posvol]
        dark=sectiondata[negvol]

        bwidth=30
        dwidth=46    
        brightsq=np.reshape(bright[:bright.size/bwidth*bwidth],(bright.size/bwidth,bwidth))

        if dark_cycle_reduc == True:
            
            # deleting fake data points from dark cycles 
            dark2=np.delete(np.array(dark), np.arange(46,np.array(dark).size,49))
            dark3=np.delete(np.array(dark2), np.arange(46,np.array(dark2).size,48))
            dark4=np.delete(np.array(dark3), np.arange(46,np.array(dark3).size,47))
            dark=dark4

            # to compensate for dark cycle reduction, bright cycle must be reduced in length
            brightsq_reduced=brightsq[0::4]

        darksq=np.reshape(dark[:dark.size/dwidth*dwidth],(dark.size/dwidth,dwidth))



 #dark2=np.delete(np.array(dark), np.arange(46,np.array(dark).size,49))
 #dark3=np.delete(np.array(dark2), np.arange(46,np.array(dark2).size,48))
 #dark4=np.delete(np.array(dark3), np.arange(46,np.array(dark3).size,47))

 #darksq2=np.reshape(dark2[:dark2.size/dwidth*dwidth],(dark2.size/dwidth,dwidth))
 #darksq3=np.reshape(dark3[:dark3.size/dwidth*dwidth],(dark3.size/dwidth,dwidth))
 #darksq4=np.reshape(dark4[:dark4.size/dwidth*dwidth],(dark4.size/dwidth,dwidth))
       	
        mincycles=min((bright.size/4)/bwidth,dark.size/dwidth)

        if dark_cycle_reduc == True:
            brightpoint=(brightsq_reduced[:mincycles,4]+brightsq_reduced[:mincycles,5]+brightsq_reduced[:mincycles,6])/3
            lastpoint=brightsq_reduced[:mincycles,-1]
        else:
            brightpoint=(brightsq[:mincycles,4]+brightsq[:mincycles,5]+brightsq[:mincycles,6])/3
            lastpoint=brightsq[:mincycles,-1]

        brightslope = (brightpoint - lastpoint)/24
        megapore = [brightslope > 0.4]
        megapore = np.asarray(megapore)
        megapore = contiguous_regions(megapore,3)

        highFPD = np.ndarray.flatten(brightsq)
        highFPD = contiguous_regions((highFPD>254),8)


    #    avgbright = np.asarray([np.mean(x) for x in brightsq[:mincycles,4:]])
     #   avgbrightbincount = plt.hist(avgbright)[0]
      #  avgbrightbinvalues = plt.hist(avgbright)[1]
       # diffcount = np.diff(avgbrightbincount)
        #diffbool = [np.where(diffcount>0)]
        
         
        #actualbright = np.asarray(brightsq[:mincycles,4:])


        #find pores,shorts:
        if dark_cycle_reduc == True:
            lpd = brightsq_reduced[:mincycles,-1]-darksq[:mincycles,-1]
        else:
            lpd = brightsq[:mincycles,-1]-darksq[:mincycles,-1]
        #charge=np.cumsum(lpd)
        chargelist.append(np.sum(lpd))
        fpd= brightsq_reduced[:mincycles,0]-darksq[:mincycles,0]
        sdd=darksq[:mincycles,9]-darksq[:mincycles,4]
        deac=np.where(brightsq_reduced[:mincycles,0]==0)
        deacn=brightsq_reduced[:mincycles,0]==0
        deactest=contiguous_regions(deacn,3)
        short =contiguous_regions(sdd>5,3)
        pore = contiguous_regions(((lpd>10)&(lpd<200)),3)
        singlepore = contiguous_regions(((lpd>=55)&(lpd<=85)),3)
        multipore = contiguous_regions(((lpd>90)&(lpd<200)),3)
        blm = contiguous_regions((lpd<5),3)
#        megapore=contiguous_regions((lpd>100),11)
        #if deac:
            #get time
            #get cause


#did it deactivated?
        if (deactest>-1): #yes, deactivated
            DEACTIVATED.append(1)

#what was it before it deactivated?
            if (pore>-1): #was it ever a pore?
                PORE.append(1)
                if (singlepore>-1): #was it ever a single pore?
                    SINGLEPORE.append(1)
                elif (singlepore==-1):
                    SINGLEPORE.append(0)
                if (multipore>-1): #was it ever a multipore?
                    MULTIPORE.append(1)
                elif (multipore==-1):
                    MULTIPORE.append(0)
            elif (pore==-1):
                PORE.append(0)
                SINGLEPORE.append(0)
                MULTIPORE.append(0)
            if (megapore>-1) or (highFPD>-1): #was it ever a megapore?
                MEGAPORE.append(1)
            else:
                MEGAPORE.append(0)
            if (short>-1): #did it short?
                SHORTED.append(1)
                BILAYER.append(0)
            else:
                SHORTED.append(0)
                BILAYER.append(1)

        else: #no, did not deactivate
            DEACTIVATED.append(0)
            SHORTED.append(0)
            BILAYER.append(0)
            if (pore>-1): #was it ever a pore?
                PORE.append(1)
                if (singlepore>-1): #was it ever a single pore?
                    SINGLEPORE.append(1)
                elif (singlepore==-1):
                    SINGLEPORE.append(0)
                if (multipore>-1): #was it ever a multipore?
                    MULTIPORE.append(1)
                elif (multipore==-1):
                    MULTIPORE.append(0)
            elif (pore==-1):
                PORE.append(0)
                SINGLEPORE.append(0)
                MULTIPORE.append(0)
            if (megapore>-1) or (highFPD>-1): #was it ever a megapore?
                MEGAPORE.append(1)
            else:
                MEGAPORE.append(0)



        if (deactest>-1):
            deactivationtime=deactest
         #   DEACTIVATED.append(1)            
            # check if it was ever a pore:
            if (pore>-1):
                #PORE.append(1)
                if (megapore>-1) or (highFPD>-1): # if megapore
                    #MEGAPORE.append(1)
                    if (short>-1) :#short
                        shorttime =short
                        deactivationcauselist.append(4)
                       # SHORTED.append(1)
                        #BILAYER.append(0)
                    else: #bilayer   
                        shorttime=-1
                        deactivationcauselist.append(3)
                        #SHORTED.append(0)
                        #BILAYER.append(1)
                else: #no megapore in pore in deac
                    #MEGAPORE.append(0)
                    if (short>-1) :#short
                        shorttime =short
                     #   SHORTED.append(1)
                      #  BILAYER.append(0)
                        deactivationcauselist.append(2)
                    else: #bilayer
                        shorttime=-1
                       # SHORTED.append(0)
                        #BILAYER.append(1)
                        deactivationcauselist.append(1)
                       
                        #fig=plt.figure()
                        #plt.plot(sectionts[(deactivationtime*76)-3000:(deactivationtime*76)+3000],sectiondata[(deactivationtime*76)-3000:(deactivationtime*76)+3000],linestyle='None', marker='.')
                        #plt.plot(sectionts[(deactivationtime*76)-3000:(deactivationtime*76)+3000],sectionvapp[(deactivationtime*76)-3000:(deactivationtime*76)+3000],linestyle='None', marker='x')
                        ##plt.plot(sectionts[(megapore*76)-3000:(megapore*76)+3000],sectiondata[(megapore*76)-3000:(megapore*76)+3000],linestyle='None', marker='.')
                        #
                        ##plt.plot(sectionts[(megapore*76)-3000:(megapore*76)+3000],sectionvapp[(megapore*76)-3000:(megapore*76)+3000],linestyle='None', marker='x')
                        #
                        ##plt.plot(sectionts,sectiondata,linestyle='None', marker='.')
                        ##plt.plot(sectionts,sectionvapp,linestyle='None', marker='x')
                        #
                        #plt.title("1D trace for %s"%(cellOI))
                        #figManager = plt.get_current_fig_manager()
                        #figManager.window.showMaximized()    
                        #plt.show()
                        #fig=plt.figure()
                        ##plt.plot(sectionts[(deactivationtime*76)-3000:(deactivationtime*76)+3000],sectiondata[(deactivationtime*76)-3000:(deactivationtime*76)+3000],linestyle='None', marker='.')
                        ##plt.plot(sectionts[(deactivationtime*76)-3000:(deactivationtime*76)+3000],sectionvapp[(deactivationtime*76)-3000:(deactivationtime*76)+3000],linestyle='None', marker='x')
                        #plt.plot(sectionts[(pore*76):(pore*76)+3000],sectiondata[(pore*76):(pore*76)+3000],linestyle='None', marker='.')
                        #plt.plot(sectionts[(pore*76):(pore*76)+3000],sectionvapp[(pore*76):(pore*76)+3000],linestyle='None', marker='x')
                        #
                        ##plt.plot(sectionts,sectiondata,linestyle='None', marker='.')
                        ##plt.plot(sectionts,sectionvapp,linestyle='None', marker='x')
                        #
                        #plt.title("1D trace for %s"%(cellOI))
                        #figManager = plt.get_current_fig_manager()
                        #figManager.window.showMaximized()    
                        #plt.show()
                        #
                        #fig=plt.figure()
                        #plt.subplot(1,2,1)
                        ##plt.title("Vliq=%d,Outliers:%d"%(np.max(sectionvapp),np.sum(brightoutliers)))
                        #plt.imshow(brightsq[deactivationtime-50:deactivationtime+50,:],interpolation='nearest')
                        #plt.colorbar() 
                        #plt.subplot(1,2,2)
                        #plt.title("Vliq=%d"%(np.min(sectionvapp)))
                        #plt.imshow(darksq[deactivationtime-50:deactivationtime+50,:],interpolation='nearest')
                        #plt.colorbar()
                        #figManager = plt.get_current_fig_manager()
                        #figManager.window.showMaximized()    
                        #plt.show()
                        #break
                        
            else: # no pore found
                #but megapore
              #  PORE.append(0)
               # SINGLEPORE.append(0)
                #MULTIPORE.append(0)
                if (megapore>-1) or (highFPD>-1): # if megapore
                 #   MEGAPORE.append(1)
                    if (short>-1) :#short
                        shorttime =short
                  #      SHORTED.append(1)
                   #     BILAYER.append(0)
                        deactivationcauselist.append(6)
                    else: #bilayer   
                        shorttime=-1
                    #    SHORTED.append(0)
                     #   BILAYER.append(1)
                        deactivationcauselist.append(5)
                
                else: #no pore or megapore
                  #  MEGAPORE.append(0)
                    if (short>-1) :#short
                        shorttime =short
                   #     SHORTED.append(1)
                    #    BILAYER.append(0)
                        deactivationcauselist.append(-2)
                    else: #bilayer
                        shorttime=-1
                     #   SHORTED.append(0)
                      #  BILAYER.append(1)
                        deactivationcauselist.append(-1)
        else:  #pore or multi pore
            deactivationtime=-1
            shorttime=-1
            #SHORTED.append(0)
            #BILAYER.append(0)
            #DEACTIVATED.append(0)
            if (megapore>-1) or (highFPD>-1):
                deactivationcauselist.append(10)
             #   MEGAPORE.append(1)
              #  PORE.append(0)
            else:
                deactivationcauselist.append(9)
               # MEGAPORE.append(0)
                #PORE.append(1)
        sectiondata=None
        deactivationtimelist.append(deactivationtime)
        shorttimelist.append(shorttime)
        vliqlist.append(np.max(sectionvapp))
        
        print cellOI + ' DEAC:' + str(deactest) + ' Short:' +str(short) +' Pore:' +str(pore) + ' Megapore:' +str(megapore) +' Class:' +str(deactivationcauselist[-1])
        print cellOI,DEACTIVATED,BILAYER,SHORTED,PORE,MEGAPORE,SINGLEPORE,MULTIPORE
        
    dfvliq[cellOI]=vliqlist
    dfactivetime[cellOI]=activetimelist
    dfdeactivationtime[cellOI]=deactivationtimelist
    dfshorttime[cellOI]=shorttimelist
    dfdeactivationcause[cellOI]=deactivationcauselist
    dfcharge[cellOI]=chargelist
    dfMEGAPORE[cellOI] = MEGAPORE
    dfPORE[cellOI] = PORE
    dfSHORTED[cellOI] = SHORTED
    dfBILAYER[cellOI] = BILAYER
    dfDEACTIVATED[cellOI] = DEACTIVATED
    dfSINGLEPORE[cellOI] = SINGLEPORE
    dfMULTIPORE[cellOI] = MULTIPORE
    cell_idx = np.where(cell_anno_df['cell_id']==cellOI)[0][0]
    dfEnthoughtSP[cellOI] = [cell_anno_df['single_pore (unitless)'][cell_idx]]
    dfEnthoughtSeq[cellOI] = [cell_anno_df['align_good_sequencing_pore_150909 (unitless)'][cell_idx]]


    data=None

path = datapath.split('/')[-1]
frames=[dfdeactivationtime,dfdeactivationcause,dfDEACTIVATED,dfSHORTED,dfBILAYER,dfPORE,dfSINGLEPORE,dfMULTIPORE,dfMEGAPORE,dfEnthoughtSP, dfEnthoughtSeq]
dfsummary=pd.concat(frames,keys=['DeacTime','Class','deactivated','shorted','bilayer','pore','singlepore','multipore','megapore','Enthought SP', 'Enthought Seq'])   
rowID = []
columnID  = []
bankID = []
for X in dfsummary.keys().values:
    rowID.append(int(X.split('c')[-1])/128)
    columnID.append(int(X.split('c')[-1])%128)
    bankID.append(X.split('c')[0].split('b')[-1])              
dfvliq=dfvliq.transpose()
#dfvliq.to_csv(path+'_vliq.csv')
dfactivetime=dfactivetime.transpose()
#dfactivetime.to_csv(path+'_activetime.csv')
dffpd=dffpd.transpose()        
#dffpd.to_csv(path+'_fpd.csv')

dfdeactivationtime=dfdeactivationtime.transpose()
#dfdeactivationtime.to_csv(path+'_DeactivationTime.csv')
dfdeactivationcause=dfdeactivationcause.transpose()
#dfdeactivationcause.to_csv(path+'_DeactivationCause.csv')
dfcharge=dfcharge.transpose()
#dfcharge.to_csv(path+'_charge.csv')
date = path.split('_')[0]
station = path.split('_')[-2]
group = path.split('_')[1]
run = path.split('_')[2]
chip = path.split('_')[-1]
dfsummary=dfsummary.transpose()
dfsummary['row'] = pd.Series(np.asarray(rowID),index=dfsummary.index)
dfsummary['column'] = pd.Series(np.asarray(columnID),index=dfsummary.index)
dfsummary['bank'] = pd.Series(np.asarray(bankID),index=dfsummary.index)
dfsummary['date'] = pd.Series(date,index=dfsummary.index)
dfsummary['group'] = pd.Series(group,index=dfsummary.index)
dfsummary['run'] = pd.Series(run,index=dfsummary.index)
dfsummary['station'] = pd.Series(station,index=dfsummary.index)
dfsummary['chipID'] = pd.Series(chip,index=dfsummary.index)
dfsummary.to_csv(path+'_multiubf_summary.csv')
print dfsummary

#Counting Classes:
one = [dfsummary['Class']==1]
one = len(np.where(one[0]==True)[0])
two = [dfsummary['Class']==2]
two = len(np.where(two[0]==True)[0]) 
three = [dfsummary['Class']==3]
three = len(np.where(three[0]==True)[0]) 
four = [dfsummary['Class']==4]
four = len(np.where(four[0]==True)[0]) 
five = [dfsummary['Class']==5]
five = len(np.where(five[0]==True)[0]) 
six = [dfsummary['Class']==6]
six = len(np.where(six[0]==True)[0]) 
nine = [dfsummary['Class']==9]
nine = len(np.where(nine[0]==True)[0]) 
ten = [dfsummary['Class']==10]
ten = len(np.where(ten[0]==True)[0]) 
negone = [dfsummary['Class']==-1]
negone = len(np.where(negone[0]==True)[0]) 
negtwo = [dfsummary['Class']==-2]
negtwo = len(np.where(negtwo[0]==True)[0]) 

print "1: "+str(one)
print "2: "+str(two)
print "3: "+str(three)
print "4: "+str(four)
print "5: "+str(five)
print "6: "+str(six)
print "9: "+str(nine)
print "10: "+str(ten)
print "-1: "+str(negone)
print "-2: "+str(negtwo)

os.system('touch ' + path + '_multiubf_megapore_metrics.txt')
with open(path + '_multiubf_megapore_metrics.txt','w') as w:
    w.write('---- Megapore Metrics ----\n')
    w.write('Class 1 (pore to bilayer): ' + str(one))
    w.write('\nClass 2 (pore to short): ' + str(two))
    w.write('\nClass 3 (pore/megapore to bilayer): ' + str(three))
    w.write('\nClass 4 (pore/megapore to short): ' + str(four))
    w.write('\nClass 5 (megapore to bilayer): ' + str(five))
    w.write('\nClass 6 (megapore to short): ' + str(six))
    w.write('\nClass 9 (pore, no deac): ' + str(nine))
    w.write('\nClass 10 (megapore, no deac): ' + str(ten))
    w.write('\nClass -1 (bilayer only): ' + str(negone))
    w.write('\nClass -2 (short only): ' +str(negtwo))

#got the cell list function
def cell_list_to_df(cellsOI):
    cellsOI=[x.upper() for x in cellsOI]
    celldf=pd.DataFrame()
    celldf['cells']=cellsOI
    celldf['bank']=celldf.cells.apply(lambda x: int(x.split('B')[1].split('C')[0]))
    celldf['c']    =celldf.cells.apply(lambda x: int(x.split('C')[1]))
    celldf['row']  =(celldf['bank']//4)*64 + (celldf['c']//128)
    celldf['col']  =(celldf['bank']%4)*128 + (celldf['c']%128)
    return celldf

#converting the first point deltas to a chipmap
cdf=cell_list_to_df(f['cells'].keys())
cdf.set_index(['cells'],drop=True,inplace=True)
tdf=df.transpose()
df2=pd.concat([cdf,tdf],join='outer',axis=1)

#converting the last point delta to a chipmap
cdflast=cell_list_to_df(f['cells'].keys())
cdflast.set_index(['cells'],drop=True,inplace=True)
tdflast=dflast.transpose()
df2last=pd.concat([cdflast,tdflast],join='outer',axis=1)
#pdb.set_trace()

chipmap=np.zeros([256,512])*np.nan
chipmaplast=np.zeros([256,512])*np.nan

endtime = datetime.datetime.now().time()

print "Started @: "+str(starttime)
print "Ended @: " + str(endtime)

