import matplotlib
matplotlib.use('Agg')
import h5py; import numpy as np; import pandas as pd; import pickle as pk; import matplotlib.pyplot as plt; import os; import seaborn as sns
scriptlocation = os.path.dirname(os.path.realpath(__file__))
runName = str(scriptlocation)
f = h5py.File(runName + '/raw/oc_calibration.h5','r')
t = f['other_traces/timestamp'].value; t=t-t[0]; t=t/1e6

oc_pores = []
oc_singlepores = []

sub = ((t[-1]-100)<t) & (t<t[-1])
print [t[0],t[-1],np.sum(sub),len(sub)]
res = []
for x in f['cells'].keys():
	y = f['cells/'+x].value
        print y
	if len(y)==0:
		continue
	z = y[sub]
	if len(z)==0:
		continue
	res.append({'cell':x,'05th':np.percentile(z,5),'15th':np.percentile(z,15),'35th':np.percentile(z,35),'55th':np.percentile(z,55),'65th':np.percentile(z,65),'85th':np.percentile(z,85),'95th':np.percentile(z,95)})
	print x
df = pd.DataFrame(res)
df['bigdelta']=df['95th']-df['05th']
df['littledelta']=df['85th']-df['15th']
df['posdroop']=df['95th']-df['65th']
df['negdroop']=df['55th']-df['15th']
sub = (0<df.littledelta) & (df.littledelta<255)
df = df[sub]
bank = [x.split('c')[0] for x in df.cell]
df['bank'] = bank
cwd = runName.split('/')[-1]
fn = runName + '/' + cwd + '_oc_calibration500' 
pk.dump(df,open(fn+'.p','wb'))
df.to_csv(fn+'.csv')
banks = list(set(df.bank))

littledeltas = []
cells = []
rowID = []
columnID = []
laneID = []


for x in banks:
	#creating negdroop/littledelta distributions for oc calibration
	sub = (df.bank==x) & (df.negdroop<7)
	sub2 = (df.bank==x)
        tmp = np.array(df.littledelta[sub])
        pik = np.sum( (25<=tmp) & (tmp<200) )
	pik2 = np.sum( (25<=tmp) & (tmp<60) )
	oc_pores.append(pik)
	oc_singlepores.append(pik2)
        plt.hist(tmp,bins=50,range=(0,250))
        plt.grid()
        plt.title(x+' '+cwd)
        plt.xlabel('25 <= '+str(pik)+' < 200, negdroop < 7')
        plt.ylabel(str(np.sum(sub))+' Total')
        plt.savefig(runName + '/' + x+'_'+cwd+'_'+str(pik)+'_'+str(np.sum(sub))+'_oc_calibration_pores.png')
	plt.xlabel('25 <= '+str(pik2)+' < 60, negdroop < 7')
	plt.ylabel(str(np.sum(sub))+' Total')
	plt.savefig(runName +'/'+ x+'_'+cwd+'_'+str(pik2)+'_'+str(np.sum(sub))+'_oc_calibration_singlepores.png')
        plt.close()
	plt.scatter(df.littledelta[sub2],df.negdroop[sub2])
	plt.grid()
	plt.title(x+' '+cwd)
	plt.xlabel('littledelta')
	plt.ylabel('negdroop')
	plt.xlim([0,200])
	plt.ylim([0,20])
        plt.savefig(runName + '/' + x + '_'+cwd+'_oc_calibration_scatter0500.png')
	plt.close()



for idx,y in enumerate(df['littledelta'].values):
	littledelta = y
	negdroop = df.iloc[idx]['negdroop']
	cell_id = df.iloc[idx]['cell']
	bank_id = df.iloc[idx]['bank']
	Quarter = int(bank_id.split('b')[-1])%4
	lane_id = int(bank_id.split('b')[-1])/4
	row_id = (int(cell_id.split('c')[-1])/128)+(int(lane_id)*74)
	column_id = (int(cell_id.split('c')[-1])%128) + (Quarter*138)
	
	if negdroop < 7:
		rowID.append(row_id)
		columnID.append(column_id)
		littledeltas.append(littledelta)
	#print cell_id, negdroop, littledelta


#frames = np.nans((286,542),dtype=int)

frames = np.empty((286,542))
frames[:] = np.NAN

d = {
	'row' : pd.Series(rowID),
	'column' : pd.Series(columnID),
	'littledelta' : pd.Series(littledeltas)
	}

DF = pd.DataFrame(d)
frames[DF['row'],DF['column']] = DF['littledelta']
heatmap = sns.heatmap(frames,cmap = plt.cm.rainbow,vmin = 0, vmax = 255)
filename = runName.split('/')[-1]
plt.axis('off')
plt.title(filename)
plt.xlabel('Columns')
plt.ylabel('Rows')
plt.savefig(filename + '_littledelta_heatmap.png')
plt.close()

os.system('touch oc_calibrationyield_heatmap_v2.py.COMPLETE.flag')

print('Finished')
