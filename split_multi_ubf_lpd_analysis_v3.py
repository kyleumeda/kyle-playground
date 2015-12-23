import h5py,os,sys,glob,time,math
import numpy as np; import pandas as pd; import matplotlib.pyplot as plt
import multiprocessing

def find_anno_file():
    P_dir=glob.glob("P*")[-1]
    anno_file=P_dir+'/annotations.h5'
    return anno_file

def get_annotation_df(anno_file):
    f=h5py.File(anno_file)
    cells=f['experiments/primary/cells'].value
    single_pores=f['experiments/primary/cell_anno/single_pore'].value
    seq_pores=f['experiments/primary/cell_anno/align_good_sequencing_pore_150909'].value
    deac_time=f['experiments/primary/cell_anno/primary_deactivation_time'].value
    life_time=f['experiments/primary/cell_anno/normalization_stop_time'].value
    pos_high=f['experiments/oc_calibration/cell_anno/oc_calibration_pos_high'].value
    df=pd.DataFrame(index=cells)
    df['single pore']=single_pores
    df['seq pore']=seq_pores
    df['deac time']=deac_time
    df['single pore end']=life_time
    df['pos high']=pos_high    
    return df

def find_zero_time(data,ts):
    zero_idx=np.where(data>0)[0][-1]
    if zero_idx != (len(ts)-1):
        zero_time=ts[zero_idx]
        deactivated=True
    else:
        zero_time = np.nan
        zero_idx = np.nan
        deactivated=False
    return zero_time,zero_idx,deactivated

def deac_type(neg_data,neg_time,ts):
    zero_idx=np.where(neg_data>0)[0][-1]
    if zero_idx != (len(ts)-1):
        neg_zero_data=neg_data[zero_idx-137:zero_idx+1].reshape(3,46)
        shortdarkdecay_all=neg_zero_data[:3,9]-neg_zero_data[:3,4]
        shortdarkdecay=np.mean(shortdarkdecay_all)
        if shortdarkdecay > 5:
            deac_cause = 'short'
        else:
            deac_cause = 'bilayer'
    else:
        deac_cause = 'None'
        shortdarkdecay = np.nan
    return deac_cause,shortdarkdecay

 
def compensate_for_dark_reduc(neg_data,pos_data,neg_time,pos_time):
   
    # deleting fake data points from dark cycles 
    dark=neg_data
    dark2=np.delete(np.array(dark), np.arange(46,np.array(dark).size,49))
    dark3=np.delete(np.array(dark2), np.arange(46,np.array(dark2).size,48))
    dark4=np.delete(np.array(dark3), np.arange(46,np.array(dark3).size,47))
    neg_data=dark4

    # deleting fake data points from dark time stamps
    start_neg_time=neg_time
    darkt2=np.delete(np.array(neg_time),np.arange(46,np.array(neg_time).size,49))
    darkt3=np.delete(np.array(darkt2),np.arange(46,np.array(darkt2).size,48))
    darkt4=np.delete(np.array(darkt3),np.arange(46,np.array(darkt3).size,47))
    neg_time=darkt4

    bwidth=30
    dwidth=46
    bright=pos_data
    brightsq=np.reshape(bright[:bright.size/bwidth*bwidth],(bright.size/bwidth,bwidth))
    brightsq_reduced=brightsq[0::4]
    bright_t=pos_time
    brightsq_t=np.reshape(bright_t[:bright_t.size/bwidth*bwidth],(bright_t.size/bwidth,bwidth))
    brightsq_t_reduced=brightsq_t[0::4]
    pos_data=brightsq_reduced.flatten()
    pos_time=brightsq_t_reduced.flatten()

    return neg_data,pos_data,neg_time,pos_time

def find_LPD_med_before_deac(pos_data,pos_zero_idx,pos_time,neg_data,neg_zero_idx,neg_time):
    if math.isnan(pos_zero_idx) == True or math.isnan(neg_zero_idx) == True:
        lpd_med_before_deac = np.nan
    else:
        pos_zero_time=pos_time[pos_zero_idx]-0.4
        pos_interval=np.where((pos_time<pos_zero_time)&(pos_time>(pos_zero_time-85)))[0]
        pos_lst_pt=np.percentile(pos_data[pos_interval],1)

        neg_zero_time=neg_time[neg_zero_idx]-0.4
        neg_interval=np.where((neg_time<neg_zero_time)&(neg_time>(neg_zero_time-85)))[0]
        neg_lst_pt=np.percentile(neg_data[neg_interval],99)

        lpd_med_before_deac=pos_lst_pt-neg_lst_pt

    return lpd_med_before_deac
    

def find_lpd_at_single_pore_end(pos_data,neg_data,pos_time,neg_time,norm_idx,singlepore,lasttimestamp):    
    
    if math.isnan(norm_idx)==False:

        pos_end=np.where(pos_time<norm_idx)[0][-1]
        neg_end=np.where(neg_time<norm_idx)[0][-1]
        pos_end_time=pos_time[pos_end]
        neg_end_time=neg_time[neg_end]

        pos_interval_before=np.where((pos_time<pos_end_time)&(pos_time>(pos_end_time-90)))[0]
        neg_interval_before=np.where((neg_time<neg_end_time)&(neg_time>(neg_end_time-90)))[0]
        pos_lst_pt_before=np.percentile(pos_data[pos_interval_before][pos_data[pos_interval_before]>0],10)
        neg_lst_pt_before=np.percentile(neg_data[neg_interval_before],90)
        lpd=pos_lst_pt_before-neg_lst_pt_before
    
        if norm_idx < lasttimestamp-95:
            # three cycles after normalization stops
            pos_interval_after=np.where((pos_time>pos_end_time)&(pos_time<(pos_end_time+90)))[0]
            neg_interval_after=np.where((neg_time>neg_end_time)&(neg_time<(neg_end_time+90)))[0]
            try:            
                pos_lst_pt_after=np.percentile(pos_data[pos_interval_after][pos_data[pos_interval_after]>0],10)
            except:
                pos_lst_pt_after = 0
            neg_lst_pt_after=np.percentile(neg_data[neg_interval_after],90)
            lpd_after=pos_lst_pt_after-neg_lst_pt_after

        else:
            lpd_after=np.nan

        special_pore=False

    elif math.isnan(norm_idx)==True and singlepore==True:
        
        pos_end_time=pos_time[lasttimestamp]
        neg_end_time=neg_time[lasttimestamp]        

        pos_interval_before=np.where((pos_time<pos_end_time)&(pos_time>(pos_end_time-90)))[0]
        neg_interval_before=np.where((neg_time<neg_end_time)&(neg_time>(neg_end_time-90)))[0]
        pos_lst_pt_before=np.percentile(pos_data[pos_interval_before],1)
        neg_lst_pt_before=np.percentile(neg_data[neg_interval_before],99)
        lpd=pos_lst_pt_before-neg_lst_pt_before
        lpd_after=np.nan      
        special_pore=True        

    else:
        lpd = np.nan
        lpd_after = np.nan
        special_pore=False

    return lpd,lpd_after,special_pore

def execute_bank_analysis(scriptlocation,anno_file_path,bank_raw_path,bank_id):
    runDir = str(scriptlocation)
    fileName = runDir.split('/')[-1]
    stationName = fileName.split('_')[-2]
    runNum = fileName.split('_')[-3]
    chipID = fileName.split('_')[-1]
    groupID = fileName.split('_')[-4]
    date = fileName.split('_')[0]
    f=h5py.File(bank_raw_path)
    F=h5py.File(anno_file_path)
    liq=f['other_traces']['liq'].value;liq=(liq/(2.0**16-1.))*(7./8.)*(2.048);
    v_pre=f['other_traces']['v_pre'].value;v_pre=(v_pre/(2.0**16-1.))*(7./8.)*(2.048);
    ts=f['other_traces/timestamp'].value;ts=ts-ts[0];ts=ts/1e6
    vapp=(v_pre-liq)*1000.;vapp=vapp.astype(np.float64)
    not_vm_zero=np.where(liq!=liq[0])[0][0]
    exp_duration=ts[-1]-ts[not_vm_zero]
    lasttimestamp=ts[-1]
    normalization_end=F['experiments/primary/cell_anno/normalization_stop_time'].value
    if np.where((f['other_traces/averaged_frame_mask'].value)>0):
        dark_cycle_reduc = True
        print 'Dark Cycle Reduction enabled'
    else:
        dark_cycle_reduc = False
        print 'Dark Cycle Reduction disabled'
    posvol=np.where(vapp>0)[0] #positive regions
    negvol=np.where(vapp<0)[0] #negative regions
    res=[]
    
    cell_list = ['b05c6560']
    cellcount=0
    cells = f['cells'].keys()
    cells.sort()
#    for cell_idx,cell in cell_list:
    for cell_idx,cell in enumerate(cells):
        print cell,cell_idx
        if cellcount>2000000:
            break
        cellcount=cellcount+1
        start_time=time.time()
        data=f['cells/'+cell].value
        print ('=== file opened in %s seconds ===' % (time.time()-start_time))
        pos_data=data[posvol]
        neg_data=data[negvol]
        pos_time=ts[posvol]
        neg_time=ts[negvol]
        cell_idx=np.where(F['cells'].value==cell)[0][0]

        if (np.where(neg_data>0)[0].size == 0) or (np.where(pos_data>0)[0].size == 0):

            print ("--- %s seconds ---" % (time.time() - start_time))
            print '+++ SKIPPED +++'
            continue
    #    elif np.where(neg_data[np.where(neg_time<neg_time[0]+200)]>0)[0].size == 0:
     #       print '****** SKIPPED ******'
      #      print ("--- %s seconds ---" % (time.time() - start_time))
       #     continue
        else:
            print ' ====== NOT SKIPPED ====== '
            singlepore=F['experiments/primary/cell_anno/single_pore'].value[cell_idx]
            if dark_cycle_reduc==True:
                neg_data,pos_data,neg_time,pos_time=compensate_for_dark_reduc(neg_data,pos_data,neg_time,pos_time)
            else:
                pass  
            pos_zero_time,pos_zero_idx,deactivated=find_zero_time(pos_data,pos_time)
            try:
                neg_zero_time,neg_zero_idx,deactivated=find_zero_time(neg_data,neg_time)
            except:
                '               WEIRD SHIT HAPPENIN HERE           '
                print cell,cell,cell,cell,cell,cell
                continue
            if deactivated == True:
                deac_cause,shortdarkdecay=deac_type(neg_data,neg_time,ts)
                lpd_med_before_deac=find_LPD_med_before_deac(pos_data,pos_zero_idx,pos_time,neg_data,neg_zero_idx,neg_time)
            else:
                deac_cause = np.nan
                shortdarkdecay = np.nan
                lpd_med_before_deac = np.nan
            sp_end_time=normalization_end[cell_idx]
            lpd_at_end_of_sp,lpd_after_end_of_sp,special_pore=find_lpd_at_single_pore_end(pos_data,neg_data,pos_time,neg_time,sp_end_time,singlepore,lasttimestamp)
            if special_pore==True:
                sp_end_time = exp_duration
            res.append({'cell':cell,'bank':bank_id,'experiment duration':exp_duration,'deac cause':deac_cause,'pos zero time':pos_zero_time, 'neg zero time':neg_zero_time,'short dark decay median':shortdarkdecay,'single pore':singlepore,'lpd median before deac':lpd_med_before_deac,'single pore end time':sp_end_time,'lpd @ end of single pore':lpd_at_end_of_sp,'lpd after single pore dies':lpd_after_end_of_sp,'date':date,'group':groupID,'station':stationName,'chip ID':chipID})
            print ("--- %s seconds ---" % (time.time() - start_time))
        
    df=pd.DataFrame(res)
    df.to_pickle(fileName+'_b'+bank_id+'_multi_ubf_analysis.p')
        

def pickle_machine(pickle_file):
    df = pd.read_pickle(pickle_file)
    return df


if __name__ == '__main__':

    scriptlocation = os.path.dirname(os.path.realpath(__file__))
    datadir = sys.argv[1]
    P_dir = glob.glob(datadir+"/P*")[-1]
    anno_file = P_dir+'/annotations.h5'
    bank_raw_dirs = glob.glob(datadir+'/chunked_raw_data/*')
    bank_raw_dirs.sort()

    print '\n=== Starting Analysis ===\n'
    for raw_dir in bank_raw_dirs:
        bank_raw_path = raw_dir+'/raw/multi_ubf.h5'
        bank_id = bank_raw_path.split('/')[1].split('bank')[-1]
        print bank_id
        p=multiprocessing.Process(target=execute_bank_analysis,args=(scriptlocation,anno_file,bank_raw_path,bank_id,))
        p.start()
     #       processes.append(p)
     #   for p in processes:
      #      p.start()
       #     p.join()

  #      pickle_files = glob.glob("*multi_ubf_analysis.p")
   #     pickle_files.sort()
        
    
#    df_list = [pickle_machine(pickle_file) for pickle_file in pickle_files]
 #       df = pd.concat(df_list,ignore_index=True)
#        runDir = str(scriptlocation)
#        fileName = runDir.split('/')[-1]
#        df.to_csv(fileName + '_allbank_multi_ubf_lpd_analysis.csv')
#        os.system('touch multi_ubf_lpd_analysis.COMPLETE.flag')







