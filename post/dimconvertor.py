#!/root/Software/anaconda3/bin/python3
import os
from basic import *


src = '/gpfs01/cuigx2_work/whn/test_OFW180/data10/base/'
dst = 'mytest/'

os.system('mkdir '+dst)
os.system('mkdir '+dst+'fielddata/')
os.system('mkdir '+dst+'probedata/')
os.system('mkdir '+dst+'statdata/')
os.system('mkdir '+dst+'postdata/')


para = DataSetInfo(src)
		
logs = np.loadtxt(para.statpath+'LOG.dat', skiprows=3)
utau = np.mean([-log[6] for log in logs if int(log[0]) in para.tsteps])**.5
uave = np.mean([ log[9] for log in logs if int(log[0]) in para.tsteps])

uc = utau
# uc = uave

uc = 1.
lc = 22.6906404/1420.960*492.2115


for filename in os.listdir(para.statpath):

	if filename == 'LOG.dat':
		with open(para.statpath + filename) as fp:
			header = ''.join(fp.readlines()[:3]).strip()

		data = np.loadtxt(para.statpath+filename, skiprows=3)

		data[:,1]    /= lc/uc
		data[:,2:9]  /= uc**2
		data[:,9:12] /= uc
		data[:,12]   /= uc/lc

		np.savetxt(dst+'statdata/'+filename, data, header=header, comments='')

	elif filename == 'PROF.dat':
		with open(para.statpath + filename) as fp:
			header = ''.join(fp.readlines()[:3]).strip()

		data = np.loadtxt(para.statpath+filename, skiprows=3)

		data[:,0]     /= lc
		data[:,1:4]   /= uc
		data[:,4:11]  /= uc**2
		data[:,11:14] /= uc**3
		data[:,14]    /= uc**4
		# data[:,15]    /= uc # eddy viscosity not necessary

		np.savetxt(dst+'statdata/'+filename, data, header=header, comments='')

	else:
		with open(dst+'statdata/'+filename, 'w') as fp:
			with open(para.statpath + filename) as fp0:
				for line in fp0: fp.write(line)


for filename in os.listdir(para.fieldpath):
	print('converting ', filename)

	q = read_channel(para.fieldpath + filename)

	if   filename[:2] == 'PT':                 q /= uc**3/lc
	elif filename[:2] in ('UT','VT','WT','P'): q /= uc**2/lc
	elif filename[0]  in ('U','V','W'):        q /= uc

	write_channel(dst+'fielddata/'+filename, q)




