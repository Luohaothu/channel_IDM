import numpy as np
from numpy.fft import fft,ifft, fft2,ifft2, hfft,ihfft
from scipy.integrate import trapz
from struct import pack, unpack
from os import listdir, system


##### common functions #####

def spec(q): return ifft(ihfft(q), axis=-2)
def phys(q): return hfft(fft(q, axis=-2)) # Nx must be even

def write_channel(pame, q):
	ny, nz, nx = q.shape
	info = np.array([nx, ny, nz]+[0]*(2*nx*nz-3), dtype=np.int32).tobytes() # 2* for 4bit -> 8bit
	np.hstack([ unpack(nx*nz*'d', info), np.ravel(q) ]).tofile(pame)

def read_channel(pame):
	with open(pame, 'rb') as fp: nx, ny, nz = unpack('3i', fp.read(12))
	# return np.fromfile(pame, np.float64, offset=(nx*nz)*8).reshape([ny, nz, nx])
	q = np.fromfile(pame, np.float64, ).reshape([ny+1, nz, nx])
	return q[1:]

###########################


class DataSetInfo:
	def __init__(self, path):
		self.datapath = path

		self.read_XIN(self.datapath)
		self.read_GRD(self.statpath)

		self.dx, self.hx = self.__get_Interval(self.x, self.xc)
		self.dy, self.hy = self.__get_Interval(self.y, self.yc)
		self.dz, self.hz = self.__get_Interval(self.z, self.zc)

		self.Nxc, self.kx = self.__get_wavenumber(self.Nx, self.Lx)
		self.Nzc, self.kz = self.__get_wavenumber(self.Nz, self.Lz)

		self.tsteps = self.__get_tsteps(self.fieldpath)

	def read_XIN(self, path):
		with open(path + "XINDAT") as fp:
			for line in fp:
				if "//" in line: line = line[:line.index("//")]

				if "fieldpath" in line: self.fieldpath = path + line.split('=')[1].split('"')[1]
				if "probepath" in line: self.probepath = path + line.split('=')[1].split('"')[1]
				if "statpath" in line: self.statpath = path + line.split('=')[1].split('"')[1]
				if "postpath" in line: self.postpath = path + line.split('=')[1].split('"')[1]

				if "Re" in line: self.Re = float( line.split('=')[1].strip().split()[0] )
				if "dt" in line: self.dt = float( line.split('=')[1].strip().split()[0] )

				if "Nx" in line: self.Nx = int( line.split('=')[1].strip().split()[0] )
				if "Ny" in line: self.Ny = int( line.split('=')[1].strip().split()[0] )
				if "Nz" in line: self.Nz = int( line.split('=')[1].strip().split()[0] )
				if "Lx" in line: self.Lx = float( line.split('=')[1].strip().split()[0] )
				if "Ly" in line: self.Ly = float( line.split('=')[1].strip().split()[0] )
				if "Lz" in line: self.Lz = float( line.split('=')[1].strip().split()[0] )

	def read_GRD(self, path):
		x, xc = [], []
		y, yc = [], []
		z, zc = [], []
		
		try:
			with open(path + 'MESH.txt') as fp:
				for line in (line.strip() for line in fp):
					if not line: continue
					elif line[0] == 'x': cor, coc = x, xc
					elif line[0] == 'y': cor, coc = y, yc
					elif line[0] == 'z': cor, coc = z, zc
					else:
						cor.append(float(line.split()[0]))
						coc.append(float(line.split()[1]))
		except:
			x = np.hstack([0, range(self.Nx)]) / (self.Nx-1) * self.Lx
			z = np.hstack([0, range(self.Nz)]) / (self.Nz-1) * self.Lz
			y = np.loadtxt(path + "CHANNEL.GRD").ravel()

			xc = np.arange(-.5,Nx) / (self.Nx-1) * self.Lx
			zc = np.arange(-.5,Nz) / (self.Nz-1) * self.Lz
			yc = np.hstack([y[1], .5 * (y[1:-1] + y[2:]), y[-1]])

		self.x, self.xc = np.array(x), np.array(xc)
		self.y, self.yc = np.array(y), np.array(yc)
		self.z, self.zc = np.array(z), np.array(zc)

	def __get_Interval(self, y, yc):
		hy = np.empty(len(y))
		dy = np.empty(len(yc))
		
		dy[1:-1] = y[2:] - y[1:-1]
		hy[1:]   = yc[1:] - yc[:-1]
		dy[0]    = 2 * (y[1] - yc[0]);
		dy[-1]   = 2 * (yc[-1]-y[-1]);
		hy[0]    = 0

		return dy, hy

	def __get_wavenumber(self, Nx, Lx):
		Nxc = (Nx-1) // 2 + 1
		kx = np.hstack([range(Nxc), range(Nxc-Nx+1, 0)]) * (2*np.pi/Lx)
		return Nxc, kx

	def __get_tsteps(self, path):
		tsteps = []
		for name in [s for s in listdir(path) if ".bin" in s]:
			tsteps.append(int("".join([c for c in name if c.isdigit()])))
		return sorted(set(tsteps)) # set() removes duplicate elements


class Field:
	def __init__(self, para):
		self.para = para

	def read(self, name, stagtyp=4):
		return self.__to_cellcenter(
			read_channel(self.para.fieldpath + name)[:,1:-1,1:-1], # for periodic channels
			stagtyp if stagtyp in (0,1,2,3) else self.__infer_stagtyp(name)	)

	def read_mean(self, name, stagtyp=4):
		return self.read_fluc_mean(name, stagtyp)[1]

	def read_fluc(self, name, stagtyp=4):
		return self.read_fluc_mean(name, stagtyp)[0]

	def read_fluc_mean(self, name, stagtyp=4):
		q = self.read(name, stagtyp)
		qm = np.mean(q, axis=(-1,-2))
		return (q.T[:]-qm).T, qm

	def __to_cellcenter(self, q, stagtyp):
		hy, dy = self.para.hy, self.para.dy
		if stagtyp == 0: return q
		if stagtyp == 1: return .5 * (q + np.roll(q,-1,axis=-1))
		if stagtyp == 3: return .5 * (q + np.roll(q,-1,axis=-2))
		if stagtyp == 2: return np.vstack([
			[ hy[1]/dy[1] * (q[1]-q[2]) + .5 * (q[1]+q[2]) ], # linear extrapolation q to yc[0]
			.5 * (q[1:-1] + q[2:]),
			[ hy[-1]/dy[-2] * (q[-1]-q[-2]) + .5 * (q[-1]+q[-2]) ] ])
		print("\nUnknown stagger type !\n")

	def __infer_stagtyp(self, name):
		return	1 if 'U' in name else \
				2 if 'V' in name else \
				3 if 'W' in name else 0



