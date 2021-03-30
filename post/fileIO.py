import numpy as np



def write_channel(pame, q):
	ny, nz, nx = q.shape
	info = np.zeros(nx*nz, q.dtype)
	info.dtype = np.int32
	info[:3] = nx, nz, ny
	info.dtype = q.dtype
	np.concatenate((info, q), axis=None).astype(np.float64).tofile(pame)

def read_channel(pame, fmt=1):
	nx, nz, ny = np.fromfile(pame, np.int32, 3)
	if fmt==2: ny, nz = nz, ny
	q = np.fromfile(pame, np.float64, ).reshape([ny+1, nz, nx])

	if fmt==2:
		import os
		if '__trash__' not in os.listdir('.'):
			os.system('mkdir __trash__')

		os.system('mv %s __trash__/'%pame) # Linux
		# os.system('move %s __trash__/'%pame) # Windows

		write_channel(pame, q[1:]) # reorganize files of old format upon reading
	
	return q[1:]


def loadtec(filename, fmt=1, skiprows=3):
	data = np.loadtxt(filename, skiprows=skiprows)

	if skiprows == 3:
		nx, ny, nz = 0, 0, 0
	
		with open(filename) as fp:
			line = fp.readline()
			line = fp.readline()
			line = fp.readline()
			for term in line.split(','):
				if term.strip()[0] == 'i': nx = int(term.split('=')[-1])
				if term.strip()[0] == 'j': ny = int(term.split('=')[-1])
				if term.strip()[0] == 'k': nz = int(term.split('=')[-1])

		# 3D data: x, y, z, data
		if nz:
			data = data.T.reshape([-1,nz,ny,nx])
			if fmt==1: return data[0,0,0,:], data[1,0,:,0], data[2,:,0,0], data[3:]
			if fmt==2: return data[0,:,0,0], data[1,0,:,0], data[2,0,0,:], np.array([f.T for f in data[3:]])
		# 2D data: x, y, data
		elif ny:
			data = data.T.reshape([-1,ny,nx])
			if fmt==1: return data[0,0,:], data[1,:,0], data[2:]
			if fmt==2: return data[0,:,0], data[1,0,:], np.array([f.T for f in data[2:]])

	# 1D data: x, data
	return data.T[0], data.T[1:]





