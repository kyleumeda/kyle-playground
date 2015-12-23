import numpy as np, pandas as pd, matplotlib.pyplot as plt, h5py,pdb
import glob
import os, sys,datetime
import pdb
import multiprocessing


def beforedeac(data,cellOI,df2,vapp):
    vappidx=np.where(vapp>0)[0] # getting indices of all the positive values
    vappidxdiff=np.diff(vappidx) # this is getting the diff of the locations so it can be used to find where it goes from negative to positive
    vappidxdiffgrtone=vappidx[np.where(vappidxdiff>1)[0]+1] # this finds the spot where it switches from negative to positive

    temp=vappidxdiffgrtone<df2.ix[cellOI,'zero_idx']
    negvappidx=np.where(vapp<0)[0] #get negvoltage indices
    #print len(vappidxdiffgrton[temp[-3
    if (len(vappidxdiffgrtone[temp])==0):
        #pdb.set_trace()
        df2.loc[cellOI,'before_deac_start']=np.nan
        return df2
    findneg=np.where((vappidxdiffgrtone[temp][-9:][-1]<negvappidx)&(df2.ix[cellOI,'zero_idx']>negvappidx))[0] 
    if (len(findneg)==0):
        try:
            poslocs_before_deac=vappidxdiffgrtone[temp][-6]
        except IndexError:
            poslocs_before_deac=np.nan
    else:
        try:
            poslocs_before_deac=vappidxdiffgrtone[temp][-5]
        except IndexError:
            poslocs_before_deac=np.nan
    #cellOI='b5c3594'
    #if(cellOI=='b08c4345'):
    #    pdb.set_trace()

    df2.loc[cellOI,'before_deac_start']=poslocs_before_deac
    return df2

#this function just gets the locations of the first positive and negative data pts given a certain time window
def getfirstptdelta(vapp,ts,tupper,tlower):
 #   print "getfirstptdelta"
    subtsidx=np.where((ts>=tlower)&(ts<=tupper))[0]
    subvapp=vapp[subtsidx]
    
    subposvapp=np.where((vapp>0)&(ts>=tlower)&(ts<=tupper))[0]
    subnegvapp=np.where((vapp<0)&(ts>=tlower)&(ts<=tupper))[0]
    '''keeping the old code for record keeping sake
    posvappidx=np.where(vapp>0)[0];
    negvappidx=np.where(vapp<0)[0];
    subposvapp=np.intersect1d(posvappidx,subtsidx)
    subnegvapp=np.intersect1d(negvappidx,subtsidx);
    pdb.set_trace()'''
    # this plot is to make sure your code is working
    #plt.plot(ts,vapp,'b.');plt.plot(ts[subtsidx],vapp[subtsidx],'ro');plt.plot(ts[subposvapp],vapp[subposvapp],'kx',ms=10);plt.plot(ts[subnegvapp],vapp[subnegvapp],'cx',ms=10);plt.plot(ts,data,'.');plt.show();
    posvappdiff=np.diff(subposvapp)
    negvappdiff=np.diff(subnegvapp)
    firstpospts=subposvapp[np.where(posvappdiff>1)[0]+1]
    firstnegpts=subnegvapp[np.where(negvappdiff>1)[0]+1];
    #plt.plot(ts,vapp,'b.');plt.plot(ts[subtsidx],vapp[subtsidx],'rx');plt.plot(ts[firstpospts],vapp[firstpospts],'ko',ms=10);plt.plot(ts[firstnegpts],vapp[firstnegpts],'co',ms=10);plt.plot(ts,data,'.');plt.show();
    #pdb.set_trace()
    if ((len(firstnegpts)==0)|(len(firstpospts)==0)):
        return np.nan,np.nan
    if (firstnegpts[0]<firstpospts[0]):
        firstnegpts=firstnegpts[1:]
    elif (firstpospts[-1]>firstnegpts[-1]):
        firstpospts=firstpospts[:-1]
    #pdb.set_trace()
    return firstpospts,firstnegpts

#this function looks back 3 seconds before cell gets deactivated
def look3secs_before_deac(data,cellOI,dfOI,df2,ts,vapp):
    #pdb.set_trace()
    tupper=dfOI[1].zero_time-3
    tlower=dfOI[1].zero_time-5
    firstpospts,firstnegpts=getfirstptdelta(vapp,ts,tupper,tlower);
    try:
        test = len(firstpospts)
        if len(firstpospts) == len(firstnegpts):
            try:
                #pdb.set_trace();print('in try');
                deltas=data[firstpospts]-data[firstnegpts]
            except IndexError:
                #pdb.set_trace();print('in except');
                deltas=np.array([np.nan])
        elif len(firstpospts) > len(firstnegpts):
            len_diff = len(firstpospts)-len(firstnegpts)
            for x in range(len_diff):
                firstpospts = np.delete(firstpospts,-1)
            try:
                deltas=data[firstpospts]-data[firstnegpts]
            except IndexError:
                #pdb.set_trace();print('in except');
                deltas=np.array([np.nan])

        elif len(firstnegpts) > len(firstpospts):
            len_diff = len(firstnegpts) - len(firstpospts)
            for x in range(len_diff):
                firstnegpts = np.delete(firstnegpts,-1)
            try:
                deltas=data[firstpospts]-data[firstnegpts]
            except IndexError:
                #pdb.set_trace();print('in except');
                deltas=np.array([np.nan])
    except TypeError:
        deltas = np.array([np.nan])
        
    return deltas;


def fstptdelta_in_beginning(data,cellOI,dfOI,df2,ts,vapp):
    tupper=2;
    tlower=0;
    firstpospts,firstnegpts=getfirstptdelta(vapp,ts,tupper,tlower);
    deltas=data[firstpospts]-data[firstnegpts]
    #plt.plot(ts,vapp,'k.');plt.plot(ts,data,'b.');plt.plot(ts[firstpospts],vapp[firstpospts],'ko',ms=10);plt.plot(ts[firstnegpts],vapp[firstnegpts],'co',ms=10);plt.plot(ts[firstpospts],data[firstpospts],'rx',ms=10);plt.plot(ts[firstnegpts],data[firstnegpts],'kx',ms=10);plt.show();
    return deltas;
#this for calculating negdroop (taken from Roger's ocf calculator)
#def getocf(cellOI,df2,data):
#    pdb.set_trace()
#import pdb;pdb.set_trace()

def getlastptdelta(vapp,ts,tupper,tlower):
#    print "getlastptdelta"
    subtsidx=np.where((ts>=tlower)&(ts<=tupper))[0]
    subvapp=vapp[subtsidx]
    
    subposvapp=np.where((vapp>0)&(ts>=tlower)&(ts<=tupper))[0]
    subnegvapp=np.where((vapp<0)&(ts>=tlower)&(ts<=tupper))[0]
    '''keeping the old code for record keeping sake
    posvappidx=np.where(vapp>0)[0];
    negvappidx=np.where(vapp<0)[0];
    subposvapp=np.intersect1d(posvappidx,subtsidx)
    subnegvapp=np.intersect1d(negvappidx,subtsidx);
    pdb.set_trace()'''
    # this plot is to make sure your code is working
    #plt.plot(ts,vapp,'b.');plt.plot(ts[subtsidx],vapp[subtsidx],'ro');plt.plot(ts[subposvapp],vapp[subposvapp],'kx',ms=10);plt.plot(ts[subnegvapp],vapp[subnegvapp],'cx',ms=10);plt.plot(ts,data,'.');plt.show();
    posvappdiff=np.diff(subposvapp)
    negvappdiff=np.diff(subnegvapp)
    lastpospts=subposvapp[np.where(posvappdiff>1)[0]]
    lastnegpts=subnegvapp[np.where(negvappdiff>1)[0]];
    #plt.plot(ts,vapp,'b.');plt.plot(ts[subtsidx],vapp[subtsidx],'rx');plt.plot(ts[firstpospts],vapp[firstpospts],'ko',ms=10);plt.plot(ts[firstnegpts],vapp[firstnegpts],'co',ms=10);plt.plot(ts,data,'.');plt.show();
    #pdb.set_trace()
    if ((len(lastnegpts)==0)|(len(lastpospts)==0)):
        return np.nan,np.nan
    if (lastnegpts[0]<lastpospts[0]):
        lastnegpts=lastnegpts[1:]
    elif (lastpospts[-1]>lastnegpts[-1]):
        lastpospts=lastpospts[:-1]
    #pdb.set_trace()
    return lastpospts,lastnegpts    

def getlstpointdelta_end(data,cellOI,dfOI,df2,ts,vapp):
    tupper=dfOI[1]['zero_time']
    tlower=dfOI[1]['zero_time']-0.4
    lastpospoints,lastnegpts=getlastptdelta(vapp,ts,tupper,tlower)
    deltas=data[lastpospoints]-data[lastnegpts]
    return deltas

def getshortdarkdecay(vapp,ts,tupper,tlower):
#    print "get dark decay"
    subtsidx=np.where((ts>=tlower)&(ts<=tupper))[0]
    subvapp=vapp[subtsidx]
    
    subposvapp=np.where((vapp>0)&(ts>=tlower)&(ts<=tupper))[0]
    subnegvapp=np.where((vapp<0)&(ts>=tlower)&(ts<=tupper))[0]

    posvappdiff=np.diff(subposvapp)
    negvappdiff=np.diff(subnegvapp)
 #   lastpt=subnegvapp[np.where(negvappdiff>1)[0]]
    fourth_negpt=subnegvapp[np.where(negvappdiff>1)[0]]-42
    ninth_negpt=subnegvapp[np.where(negvappdiff>1)[0]]-40

    return fourth_negpt,ninth_negpt

def shortdarkdecay(data,cellOI,dfOI,df2,ts,vapp):
    tupper=dfOI[1]['zero_time']
    tlower=dfOI[1]['zero_time']-0.4
    fourth_negpt,ninth_negpt=getshortdarkdecay(vapp,ts,tupper,tlower)
    Deltas=data[ninth_negpt]-data[fourth_negpt]
    return Deltas

def execute_getpore_analysis(dataDir,anno_dir,raw_dir,bank_id):

    target=raw_dir
    runName=dataDir.split('/')[-1]
    date=runName.split('_')[0]
    group=runName.split('_')[1]
    station=runName.split('_')[-2]
    chipid=runName.split('_')[-1]
    f=h5py.File(target,'r');
    if os.path.isfile(anno_dir)==True:
        p=h5py.File(anno_dir,'r')
        analysis_complete=True
        cell_list=p['cells'].value
        DF=pd.DataFrame(index=cell_list)
        singlePores=p['experiments/primary/cell_anno/single_pore'].value
        seq_pores=p['experiments/primary/cell_anno/align_homo_good_sequencing_pore_150909'].value
        DF['single pore']=pd.Series(singlePores,index=cell_list)
        DF['seq pore']=pd.Series(seq_pores,index=cell_list)
    else:
        analysis_complete=False
        print 'No annotations.h5 found'


    #getting timestamp and vapp values
    liq=f['other_traces']['liq'].value;liq=(liq/(2.0**16-1.))*(7./8.)*(2.048);
    ts=f['other_traces/timestamp'].value;ts=ts-ts[0];ts=ts/1e6
    v_pre=f['other_traces']['v_pre'].value;v_pre=(v_pre/(2.0**16-1.))*(7./8.)*(2.048);
    vapp=(v_pre-liq)*1000.;vapp=vapp.astype(np.float64)
    df=pd.DataFrame()
    #import pdb;pdb.set_trace()
    tempzeros=np.zeros(len(ts),dtype=bool)
    rolltempzero=np.zeros(len(ts),dtype=bool)
    idx=np.arange(len(ts))

    #CELL_LIST = ['b14c5106','b14c3935']
    #CELL_LIST = ['b14c3935']
    #CELL_LIST = ['b11c1437']
    #CELL_LIST = ['b07c6243']
    #CELL_LIST = ['b12c1953']
    #CELL_LIST = ['b09c1445']
    cellcount=0
    #for cellOI in CELL_LIST:
    for cellOI in f['cells'].keys():
        cellcount=cellcount+1
        if cellcount > 200000:
            break
        print cellOI, cellcount;
        tempdf=pd.DataFrame();
        import time;
        t0=time.time()
        tempval=f['cells'][cellOI].value;
        t0=time.time()
        tempzeros[:]=tempval==0
        fst=np.where((tempzeros & ~np.hstack([False,tempzeros[:-1]]))==True)[0]
        lst=np.where((tempzeros & ~np.hstack([tempzeros[1:],False]))==True)[0]
        pr=[(i,j) for i,j in zip(fst,lst) if j>i+4]

        if (len(pr)==0):
            tempdf2=pd.DataFrame([np.nan],columns=['zero_idx'],index=[cellOI])
            tempdf2['len_zeros']=len(pr)
        else:
            tempdf2=pd.DataFrame([pr[-1][0]],columns=['zero_idx'],index=[cellOI]);
            tempdf2['len_zeros']=len(pr)
        df=pd.concat([tempdf2,df],axis=0,join='outer');

    df2=df.dropna(how='any');
    df2['zero_time']=ts[df2['zero_idx']]

    for dfOI in df2.iterrows():
        import time;t0=time.time();
        cellOI=dfOI[0]
        rowID=int(cellOI.split('c')[-1])/128
        columnID=int(cellOI.split('c')[-1])%128
        data=f['cells'][dfOI[0]].value.astype(np.float64)
        zero_idx=dfOI[1][0]
        temp_idx=np.where(vapp[zero_idx-100:zero_idx]>0)[0][-1]
        elec_voltage=vapp[zero_idx-100:zero_idx][temp_idx]
        sdd=shortdarkdecay(data,cellOI,dfOI,df2,ts,vapp); sdd=np.median(sdd)
        print cellOI
        #print dfOI[0]
        if (dfOI[1]['zero_time']>7):
            df2=beforedeac(data,cellOI,df2,vapp);
            three_sec_beforefstpt_deltas=look3secs_before_deac(data,cellOI,dfOI,df2,ts,vapp);
            df2.loc[cellOI,'row']=rowID;
            df2.loc[cellOI,'column']=columnID;
            df2.loc[cellOI,'med_delta_3_secs_before_deac']=np.median(three_sec_beforefstpt_deltas);
            fstptdeltas_in_beginning=fstptdelta_in_beginning(data,cellOI,dfOI,df2,ts,vapp);
            df2.loc[cellOI,'med_delta_in_beginning']=np.median(fstptdeltas_in_beginning);
            lstptdeltas_end=getlstpointdelta_end(data,cellOI,dfOI,df2,ts,vapp);
            df2.loc[cellOI,'med_lst_pt_delta_end']=np.median(lstptdeltas_end);
            df2.loc[cellOI,'mean_lst_pt_delta_end']=np.mean(lstptdeltas_end);
            df2.loc[cellOI,'short dark decay']=sdd
            tempsubdata=data[df2.ix[cellOI,'before_deac_start']:df2.ix[cellOI,'zero_idx']]#this will be used for calculating percentiles (used in ocf analysis)
            df2.loc[cellOI,'negdroop']=np.percentile(tempsubdata,55)-np.percentile(tempsubdata,15);
            df2.loc[cellOI,'littledelta']=np.percentile(tempsubdata,85)-np.percentile(tempsubdata,15);
            #getocf(cellOI,df2,data)
        else:
            df2=beforedeac(data,cellOI,df2,vapp);
            df2.loc[cellOI,'row']=rowID;
            df2.loc[cellOI,'column']=columnID;
            #df2.loc[cellOI,'med_delta_3_secs_before_deac']=np.nan
            #df2.loc[cellOI,'before_deac_start']=np.nan
            fstptdeltas_in_beginning=fstptdelta_in_beginning(data,cellOI,dfOI,df2,ts,vapp);
            #print fstptdeltas_in_beginning
            #if cellOI=='b7c4709':
            #    pdb.set_trace()
            lstptdeltas_end=getlstpointdelta_end(data,cellOI,dfOI,df2,ts,vapp)
            df2.loc[cellOI,'med_delta_in_beginning']=np.median(fstptdeltas_in_beginning);
            df2.loc[cellOI,'med_lst_pt_delta_end']=np.median(lstptdeltas_end);
            df2.loc[cellOI,'mean_lst_pt_delta_end']=np.mean(lstptdeltas_end);
            #tempsubdata=data[df2.ix[cellOI,'before_deac_start']:df2.ix[cellOI,'zero_idx']]#this will be used for calculating percentiles (used in ocf analysis)
            #df2.loc[cellOI,'negdroop']=np.percentile(tempsubdata,55)-np.percentile(tempsubdata,15);
            #df2.loc[cellOI,'littledelta']=np.percentile(tempsubdata,85)-np.percentile(tempsubdata,15
        if sdd > 5:
            deac_cause='shorted'
        elif df2.loc[cellOI]['med_lst_pt_delta_end']>50:
            deac_cause='megapore'
        elif 0<df2.loc[cellOI]['med_lst_pt_delta_end']<=50:
            deac_cause='pore'
        else:
            deac_cause='None'
        df2.loc[cellOI,'deactivated as']=deac_cause
        df2.loc[cellOI,'deac volt']=elec_voltage
        if analysis_complete==False:
            singlePore=np.nan
            seq_pore=np.nan
        else:
            cellidx=np.where(DF.index==cellOI)[0][0]
            singlePore=DF.iloc[cellidx]['single pore']
            seq_pore=DF.iloc[cellidx]['seq pore']
            df2.loc[cellOI,'single pore']=singlePore
            df2.loc[cellOI,'sequencing pore']=seq_pore
      #  print(time.time()-t0)
        tempsubdata=0
    df2['bank id']=pd.Series(bank_id,index=df2.index)
    df2['date']=pd.Series(date,index=df2.index)
    df2['group']=pd.Series(group,index=df2.index)
    df2['station']=pd.Series(station,index=df2.index)
    df2['chip id']=pd.Series(chipid,index=df2.index)
    df2.to_pickle(runName + '_b'+bank_id+ '_getpore_Deactivations.p')
    print "Number of Cells that had littledelta of [20,200]"
    endtime = datetime.datetime.now().time()

    print "Started @: "+str(starttime)
    print "Ended @: " + str(endtime)


if __name__ == '__main__':
    scriptlocation = os.path.dirname(os.path.realpath(__file__))
    data_dir = sys.argv[1]
    P_dir = glob.glob(data_dir+"/P*")[-1]
    anno_file_path = P_dir + '/annotations.h5'
    bank_raw_dirs = glob.glob(data_dir+'/chunked_raw_data/*')
    bank_raw_dirs.sort()
    
    print '\n=== Starting Analysis ===\n'

    for raw_dir in bank_raw_dirs:
        bank_raw_path = raw_dir+'/raw/getpore_ac.h5'
        bank_id = raw_dir.split('/')[-1].split('bank')[-1]
        print bank_id
        p=multiprocessing.Process(target=execute_getpore_analysis,args=(data_dir,anno_file_path,bank_raw_path,bank_id,))
        p.start()
#            processes.append(p)
 #       for p in processes:
  #          p.start()
   #         p.join()
      #  pickle_files = glob.glob("*.p")
       # pickle_files.sort()
        
     #   df_list = [pickle_machine(pickle_file) for pickle_file in pickle_files]
    #    df = pd.concat(df_list,ignore_index=True)
   #     runDir = str(scriptlocation)
  #      fileName = runDir.split('/')[-1]
 #       df.to_csv(fileName + '_allbank_getpore_analysis.csv')
        #os.system('touch getpore_analysis.COMPLETE.flag')

