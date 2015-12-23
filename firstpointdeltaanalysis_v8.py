import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import csv,os,sys,glob,time
import h5py
from scipy import stats
from scipy.stats import norm

# ======================================= 
# ============== Settings ===============
# =======================================

# Set to True if cell specific thinning is enabled
high_low = False
# Set to True if using Ashraf's tier thinning
tier_system = True

# ======================================= 

#defining location of files
start_time = time.ctime()
scriptlocation = os.path.dirname(os.path.realpath(__file__))
runDir = str(scriptlocation)
fileName = runDir.split('/')[-1]
stationName = fileName.split('_')[-2]
runNum = fileName.split('_')[-3]
chipID = fileName.split('_')[-1]
groupID = fileName.split('_')[-4]
date = fileName.split('_')[0]

# ======================================= 

def convert_cell_name(cell):

    if len(cell) == 8:
	    skip = True
    else:
	    skip = False
    if skip == False:
	    if len(cell.split('C')[0].split('B')[-1]) != 2:
		    bankIDprefix = str('B0')+str(cell.split('C')[0].split('B')[-1])
	    else:
		    bankIDprefix = cell.split('C')[0]
	    cellprefix = cell.split('C')[-1]
	    while len(cellprefix) != 4:
		    cellprefix = str('0') + cellprefix.split('C')[-1]
	    cellprefix = 'C'+str(cellprefix)
	    fullcellname = str(bankIDprefix)+str(cellprefix)
	    fullcellname = fullcellname.lower()
    else:
	    fullcellname = cell.lower()

    return fullcellname

def _find_max_cycles():
    FPDDatalist = []
    for file in glob.glob('reactivethinning/FPDData_*.csv'):
	    FPDDatalist.append(file)
    numcycles = len(FPDDatalist)
    return numcycles
    

def getfirstptdeltas_v2(csvfile,cell_list):
    tempdf=pd.DataFrame.from_csv(csvfile)
    oldcells = tempdf.columns.values
    renamedcells = np.asarray([convert_cell_name(cell) for cell in oldcells])
    tempdf=tempdf.T
    tempdf=tempdf.reset_index()
    del tempdf['index']
    tempdf=tempdf.set_index(renamedcells)
    tempdf=tempdf.reindex(cell_list)
    tempdf=tempdf.T
    bright=tempdf.iloc[0].values.tolist()
    dark=tempdf.iloc[1].values.tolist()
    if (tempdf.index[0] < tempdf.index[1]):	
	    deltalist=np.subtract(bright,dark,dtype='float')
    else:
	    deltalist=np.subtract(dark,bright,dtype='float')

    dfout=pd.DataFrame(deltalist).T
    dfout.columns=tempdf.keys()

    return dfout

def get_cells_from_tiers(cycle,tier_id):
    with open('reactivethinning/tier'+str(tier_id)+'_'+str(cycle)+'.cel','r') as f: 
        cells=f.read().splitlines()
    renamed_cells = [convert_cell_name(cell) for cell in cells]
    return renamed_cells

def define_tier_num(cell_id,celllist):
    for x,y in enumerate(celllist):
        if cell_id in y:
            tier_id = x+1
            break
        else:
            tier_id = np.nan
    print cell_id, '\n'        
    return tier_id

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
    df=pd.DataFrame(index=cells)
    df['single pore']=single_pores
    df['seq pore']=seq_pores
    df['deac time']=deac_time
    
    return df
    
def get_fpd_diff(cyclenum,df_list,cell_list):
    if cyclenum > 0:
        current_df = df_list[int(cyclenum)]
        previous_df = df_list[int(cyclenum)-1]
        current_FPD = current_df[current_df.columns[0]].values
        previous_FPD = previous_df[previous_df.columns[0]].values   
        FPD_diff = current_FPD - previous_FPD
    else:
        FPD_diff = 'None'
    
    df = pd.DataFrame(index=cell_list)
    df['FPD_diff_'+str(cyclenum)] = pd.Series(FPD_diff,index=df.index)
       
    return df

def cycle_df_machine(cyclenum,high_low,tier_system,full_cell_list):

    tempdf = getfirstptdeltas_v2('reactivethinning/FPDData_'+str(cyclenum)+'.csv',full_cell_list)
    df = tempdf.T
    df.columns = ['FPD_'+str(cyclenum)]  
    

    if high_low == True:
    
        hicelldf = pd.DataFrame
        hicelldf = hicelldf.from_csv('reactivethinning/hi_'+str(cyclenum)+'.cel')
        hicellList = hicelldf.index.tolist()

        lowcelldf = pd.DataFrame
        lowcelldf = lowcelldf.from_csv('reactivethinning/low_'+str(cyclenum)+'.cel')
        lowcellList = lowcelldf.index.tolist()

        hilowlist = []

        for cell in finaldf.index.values.tolist():
            if cell in hicellList:
                print cell,'high'
                hilowlist.append('high')
            elif cell in lowcellList:
                print cell,'low'
                hilowlist.append('low')
            elif cell not in hicellList or lowcellList:
                print cell,'off'			
                hilowlist.append('off')
        
        finaldf['hi-low_'+str(cyclenum)] = pd.Series(hilowlist,index=finaldf.index)

    if tier_system == True:
        tier_cell_list = [get_cells_from_tiers(cyclenum,tier_num) for tier_num in range(1,7)]
        tier_id_list = [define_tier_num(cell,tier_cell_list) for cell in full_cell_list]
        df['tier_id_'+str(cyclenum)] = tier_id_list

    return df

def define_last_values(df_list,df_diff_list):
    last_df = df_list[-1]
    last_diff_df = df_diff_list[-1]
    FPD_header = last_df.columns[0]
    tier_header = last_df.columns[-1]
    FPD_diff_header = last_diff_df.columns[0]
    last_FPD = last_df[FPD_header].values
    last_tier_id = last_df[tier_header].values
    last_FPD_diff = last_diff_df[FPD_diff_header].values

    return last_FPD, last_tier_id, last_FPD_diff
  
if __name__ == '__main__':

    numcycles = _find_max_cycles()
    f = h5py.File('raw/multi_ubf.h5','r')
    cells = f['cells'].keys()
    cells = [convert_cell_name(cell) for cell in cells]
    df = pd.DataFrame(index=cells).T

    all_dfs = [cycle_df_machine(cycle,high_low,tier_system,cells) for cycle in range(int(numcycles))]
    fpd_diffs = [get_fpd_diff(cycle,all_dfs,cells) for cycle in range(int(numcycles))]
    anno_file = find_anno_file()
    last_FPD, last_tier_id, last_FPD_diff = define_last_values(all_dfs,fpd_diffs)
    anno_df = get_annotation_df(anno_file)
    super_list = all_dfs + fpd_diffs
    super_list.append(anno_df)
    finaldf=pd.concat(super_list,axis=1)
    finaldf['last FPD'] = pd.Series(last_FPD,index=finaldf.index)
    finaldf['last tier id'] = pd.Series(last_tier_id,index=finaldf.index)
    finaldf['last FPD diff'] = pd.Series(last_FPD_diff,index=finaldf.index)
    finaldf.to_csv(fileName+'_fpd_analysis.csv')
    os.system('touch FPD_analysis.COMPLETE.flag')
    

        









