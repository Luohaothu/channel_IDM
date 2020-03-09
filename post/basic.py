import numpy as np
from numpy.fft import fft,ifft, fft2,ifft2, hfft,ihfft
from scipy.integrate import trapz
from struct import pack, unpack
import os


##### common functions #####

def spec(q):
	return ifft(ihfft(q), axis=-2)

def phys(q):
	return hfft(fft(q, axis=-2)) # Nx must be even

def write_channel(pame, q):
	ny, nz, nx = q.shape
	info = np.array([nx, ny, nz]+[0]*(2*nx*nz-3), dtype=np.int32).tobytes() # 2* for 4bit -> 8bit
	np.hstack([ unpack(nx*nz*'d', info), np.ravel(q) ]).tofile(pame)

def read_channel(pame):
	with open(pame, 'rb') as fp: nx, ny, nz = unpack('3i', fp.read(12))
	return np.fromfile(pame, np.float64, offset=(nx*nz)*8).reshape([ny, nz, nx])

###########################


class DataSetInfo:
	def __init__(self, path):
		self.datapath = path
		self.read_XIN(self.datapath)
		self.read_GRD(self.statpath, self.Ny)

		self.Nxz = self.Nx * self.Nz
		self.Nxc = self.Nx//2+1
		self.Nzc = self.Nz//2+1
		self.kx = np.hstack([range(self.Nxc), range(self.Nxc-self.Nx, 0)]) * (2*np.pi/self.Lx)
		self.kz = np.hstack([range(self.Nzc), range(self.Nzc-self.Nz, 0)]) * (2*np.pi/self.Lz)
		self.h, self.dy = self.get_Ymesh(self.y, self.yc)

		self.tsteps = self.get_tsteps(self.fieldpath)

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

	def read_GRD(self, path, Ny):
		self.y = y = np.loadtxt(path + "CHANNEL.GRD").reshape(Ny+1)
		self.yc = np.hstack([ y[[1]], .5 * (y[1:-1] + y[2:]), y[[-1]] ])

	def get_Ymesh(self, y, yc):
		h = np.empty(len(y))
		dy= np.empty(len(yc))
		h [0]    = 0
		h [1:]   = yc[1:] - yc[:-1]
		dy[0]    = 2. * (y[1] - yc[0]);
		dy[-1]   = 2. * (yc[-1]-y[-1]);
		dy[1:-1] = y[2:] - y[1:-1]
		return h, dy

	def get_tsteps(self, path):
		tsteps = []
		for name in [s for s in os.listdir(path) if ".bin" in s]:
			tsteps.append(int("".join([c for c in name if c.isdigit()])))
		return sorted(set(tsteps)) # set() removes duplicate elements


class Field:
	def __init__(self, para):
		self.para = para

	def read(self, name, stagtyp=4):
		return self.__to_cellcenter(
			read_channel(self.para.fieldpath + name),
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
		h, dy = self.para.h, self.para.dy
		if stagtyp == 0: return q
		if stagtyp == 1: return .5 * (q + np.roll(q,-1,axis=-1))
		if stagtyp == 3: return .5 * (q + np.roll(q,-1,axis=-2))
		if stagtyp == 2: return np.vstack([
			[ h[1]/dy[1] * (q[1]-q[2]) + .5 * (q[1]+q[2]) ], # linear extrapolation q to yc[0]
			.5 * (q[1:-1] + q[2:]),
			[ h[-1]/dy[-2] * (q[-1]-q[-2]) + .5 * (q[-1]+q[-2]) ] ])
		print("\nUnknown stagger type !\n")

	def __infer_stagtyp(self, name):
		return	1 if 'U' in name else \
				2 if 'V' in name else \
				3 if 'W' in name else 0


class Statis:
	def __init__(self, feld):
		self.feld = feld
		self.para = feld.para

	def calc_umean(self, tsteps=None):
		self.Um = np.zeros(self.para.Ny+1)
		if tsteps is None: tsteps = self.para.tsteps
		for tstep in tsteps:
			print("Reading umean: tstep", tstep)
			u = self.feld.read_mean("U%08i.bin"%tstep, stagtyp=1)
			self.Um += u / len(tsteps)

	def calc_statis(self, tsteps=None):
		Ny = self.para.Ny
		Nxc= self.para.Nxc
		Nzc= self.para.Nzc
		if tsteps is None: tsteps = self.para.tsteps

		self.R11, self.R22, self.R33          = (np.zeros(Ny+1) for _ in range(3))
		self.R12, self.R23, self.R13          = (np.zeros(Ny+1) for _ in range(3))
		self.Rpu, self.Rpv, self.Rpw, self.Rpp= (np.zeros(Ny+1) for _ in range(4))
		self.Um , self.Vm , self.Wm , self.Pm = (np.zeros(Ny+1) for _ in range(4))
		self.Euu, self.Evv, self.Eww, self.Epp= (np.zeros([Ny+1, Nzc, Nxc]) for _ in range(4))
		self.Euv, self.Evw, self.Euw          = (np.zeros([Ny+1, Nzc, Nxc]) for _ in range(3))

		for tstep in tsteps:
			print("Reading statis: tstep", tstep)
			u, um = self.feld.read_fluc_mean("U%08i.bin"%tstep)
			v, vm = self.feld.read_fluc_mean("V%08i.bin"%tstep)
			w, wm = self.feld.read_fluc_mean("W%08i.bin"%tstep)
			p, pm = self.feld.read_fluc_mean("P%08i.bin"%tstep)

			self.Um += um / len(tsteps)
			self.Vm += vm / len(tsteps)
			self.Wm += wm / len(tsteps)
			self.Pm += pm / len(tsteps)

			self.R11 += np.mean(u**2,axis=(-1,-2)) / len(tsteps)
			self.R22 += np.mean(v**2,axis=(-1,-2)) / len(tsteps)
			self.R33 += np.mean(w**2,axis=(-1,-2)) / len(tsteps)
			self.R12 += np.mean(u*v, axis=(-1,-2)) / len(tsteps)
			self.R23 += np.mean(v*w, axis=(-1,-2)) / len(tsteps)
			self.R13 += np.mean(u*w, axis=(-1,-2)) / len(tsteps)

			self.Rpu += np.mean(p*u, axis=(-1,-2)) / len(tsteps)
			self.Rpv += np.mean(p*v, axis=(-1,-2)) / len(tsteps)
			self.Rpw += np.mean(p*w, axis=(-1,-2)) / len(tsteps)
			self.Rpp += np.mean(p**2,axis=(-1,-2)) / len(tsteps)

			self.Euu += self.__flipk( np.abs(spec(u))**2 ) / len(tsteps)
			self.Evv += self.__flipk( np.abs(spec(v))**2 ) / len(tsteps)
			self.Eww += self.__flipk( np.abs(spec(w))**2 ) / len(tsteps)
			self.Epp += self.__flipk( np.abs(spec(p))**2 ) / len(tsteps)
			self.Euv += self.__flipk( (np.conj(spec(u))*spec(v)).real ) / len(tsteps)
			self.Evw += self.__flipk( (np.conj(spec(v))*spec(w)).real ) / len(tsteps)
			self.Euw += self.__flipk( (np.conj(spec(u))*spec(w)).real ) / len(tsteps)

	def __flipk(self, q):
		''' fold all energy to the [:nzc,:nxc] range
		    Nx must be even, as required by hft, Nz can be even or odd  '''
		nzcu = int(np.ceil(q.shape[-2]/2))
		nzcd = q.shape[-2]//2
		q.T[1:-1] *= 2 # .T usually returns a view which can act as left value
		q.T[:,1:nzcu] += q.T[:,:nzcd:-1]
		return (q.T[:,:nzcd+1]).T

	def flipy(self):
		self.Um[:] = .5 * (self.Um + self.Um[::-1])
		self.Vm[:] = .5 * (self.Vm - self.Vm[::-1])
		self.Wm[:] = .5 * (self.Wm + self.Wm[::-1])
		self.Pm[:] = .5 * (self.Pm + self.Pm[::-1])

		self.R11[:] = .5 * (self.R11 + self.R11[::-1])
		self.R22[:] = .5 * (self.R22 + self.R22[::-1])
		self.R33[:] = .5 * (self.R33 + self.R33[::-1])
		self.R12[:] = .5 * (self.R12 - self.R12[::-1])
		self.R23[:] = .5 * (self.R23 - self.R23[::-1])
		self.R13[:] = .5 * (self.R13 + self.R13[::-1])

		self.Rpu[:] = .5 * (self.Rpu + self.Rpu[::-1])
		self.Rpv[:] = .5 * (self.Rpv - self.Rpv[::-1])
		self.Rpw[:] = .5 * (self.Rpw + self.Rpw[::-1])
		self.Rpp[:] = .5 * (self.Rpp + self.Rpp[::-1])

		self.Euu[:] = .5 * (self.Euu + self.Euu[::-1])
		self.Evv[:] = .5 * (self.Evv + self.Evv[::-1])
		self.Eww[:] = .5 * (self.Eww + self.Eww[::-1])
		self.Epp[:] = .5 * (self.Epp + self.Epp[::-1])
		self.Euv[:] = .5 * (self.Euv - self.Euv[::-1])
		self.Evw[:] = .5 * (self.Evw - self.Evw[::-1])
		self.Euw[:] = .5 * (self.Euw + self.Euw[::-1])

	def calc_wallscale(self):
		''' calculate the wall scales
		    make sure Um and R12 have been correctly assigned before
		    tauw is obtained by integrating the relation ( accurancy tested to be good O(10e-4) )
		    dU^+/dy^+ - <u'v'>^+ = 1 - y\bar '''
		Ly = self.para.Ly
		Re = self.para.Re
		yc = self.para.yc
		nudUdy = np.gradient(self.Um, yc) / Re

		# self.tauw = 0.5 * (nudUdy[0] - nudUdy[-1])
		# self.tauw = np.sum( abs(nudUdy - self.R12)[1:-1] * (self.para.y[2:] - self.para.y[1:-1]) ) / (0.25 * self.para.Ly**2)
		self.tauw = 4./Ly**2 * trapz(abs(nudUdy-self.R12), yc)

		self.utau = self.tauw**.5
		self.dnu = 1./Re / self.utau
		self.tnu = self.dnu / self.utau
		self.Ret = 1./self.dnu # channel height is taken for 2.0

	def inner_scale(self):
		self.lc = self.dnu
		self.tc = self.tnu
		self.uc = self.utau
		self.pc = self.tauw

	def outer_scale(self):
		self.lc = 1.
		self.tc = 1.
		self.uc = 1.
		self.pc = 1.


class Operator:
	def __init__(self, para):
		self.para = para

	def rotation(self, u, v, w):
		''' input on staggered grids and output on cell-centers '''
		ux,vx,wx, uy,vy,wy, uz,vz,wz = self.gradient(u,v,w)
		return wy-vz, uz-wx, vx-uy # ox, oy, oz

	def gradient(self, u, v=None, w=None):
		''' input on staggered grids and output on cell-centers '''
		return self.__gradient_scalar(u) if v is None \
		  else self.__gradient_vector(u,v,w)

	def divergence(self, u, v, w):
		''' input on staggered grids and output on cell-centers '''
		Nx = self.para.Nx
		Ny = self.para.Ny
		Nz = self.para.Nz
		dx = self.para.Lx / Nx
		dz = self.para.Lz / Nz
		dy = self.para.dy.reshape([Ny+1,1,1])

		ux = (np.roll(u,-1,axis=-1) - u) / dx
		wz = (np.roll(w,-1,axis=-2) - w) / dz
		vy = np.empty([Ny+1,Nz,Nx])
		vy[1:-1] = (v[2:] - v[1:-1]) / dy[1:-1]
		vy[0], vy[-1] = vy[1], vy[-2] # vy on layer yc[0] extrapolated from layer yc[1]

		return ux + vy + wz

	def laplace(self, q):
		Nx = self.para.Nx
		Ny = self.para.Ny
		Nz = self.para.Nz
		dx = self.para.Lx / Nx
		dz = self.para.Lz / Nz
		dy = self.para.dy.reshape([Ny+1,1,1])
		h  = self.para.h. reshape([Ny+1,1,1])
		
		qxx, qyy, qzz = (np.empty([Ny+1,Nz,Nx]) for _ in range(3))

		jc = np.arange(1,Ny)
		jm, jp = jc-1, jc+1

		qxx[:] = (np.roll(q,1,axis=-1) + np.roll(q,-1,axis=-1) - 2*q) / dx**2
		qzz[:] = (np.roll(q,1,axis=-2) + np.roll(q,-1,axis=-2) - 2*q) / dz**2
		qyy[jc]= (q[jm]/h[jc] + q[jp]/h[jp] - q[jc]*(1./h[jc]+1./h[jp])) / dy[jc]
		
		qxx[0] = qxx[-1] = 0
		qyy[0] = qyy[-1] = 0
		qzz[0] = qzz[-1] = 0

		return qxx+qyy+qzz

	def __gradient_scalar(self, q):
		''' input on cell-centers and output on staggered grids '''
		Nx = self.para.Nx
		Ny = self.para.Ny
		Nz = self.para.Nz
		dx = self.para.Lx / Nx
		dz = self.para.Lz / Nz
		h  = self.para.h.reshape([Ny+1,1,1])

		u = (q - np.roll(q,1,axis=-1)) / dx
		w = (q - np.roll(q,1,axis=-2)) / dz
		v = np.empty([Ny+1,Nz,Nx])
		v[1:] = (q[1:] - q[:-1]) / h[1:] # think of q[0] on layer yc[0]
		v[0] = 0

		return u, v, w

	def __gradient_vector(self, u, v, w):
		''' input on staggered grids and output on cell-centers '''
		Nx = self.para.Nx
		Ny = self.para.Ny
		Nz = self.para.Nz
		dx = self.para.Lx / Nx
		dz = self.para.Lz / Nz
		dy = self.para.dy.reshape([Ny+1,1,1])
		h  = self.para.h. reshape([Ny+1,1,1])

		ux, vx, wx = (np.empty([Ny+1,Nz,Nx]) for _ in range(3))
		uy, vy, wy = (np.empty([Ny+1,Nz,Nx]) for _ in range(3))
		uz, vz, wz = (np.empty([Ny+1,Nz,Nx]) for _ in range(3))

		jc = np.arange(1,Ny)
		jm, jp = jc-1, jc+1
		im, ip = np.arange(-1,Nx-1), np.arange(1,Nx+1) % Nx
		km, kp = np.arange(-1,Nz-1), np.arange(1,Nz+1) % Nz

		## compute the derivatives
		# normal derivatives on cell-centers
		ux[:] = (u[:,:,ip] - u) / dx
		wz[:] = (w[:,kp,:] - w) / dz
		vy[jc]= (v[jp] - v[jc]) / dy[jc]
		# cross derivatives on staggered grids
		uz[:] = .5/dz * (u[:,kp,:] - u[:,km,:])
		wx[:] = .5/dx * (w[:,:,ip] - w[:,:,im])
		vx[:] = .5/dx * (v[:,:,ip] - v[:,:,im])
		vz[:] = .5/dz * (v[:,kp,:] - v[:,km,:])
		uy[jc]= .5/dy[jc] * ( (u[jc]*dy[jp]+u[jp]*dy[jc]) / h[jp] - (u[jc]*dy[jm]+u[jm]*dy[jc]) / h[jc] )
		wy[jc]= .5/dy[jc] * ( (w[jc]*dy[jp]+w[jp]*dy[jc]) / h[jp] - (w[jc]*dy[jm]+w[jm]*dy[jc]) / h[jc] )

		## handle the boundaries (think of boundary on yc[0])
		# vy is extrapolated from layer yc[1] to yc[0]
		vy[0], vy[-1] = vy[1], vy[-2]
		# uy,wy are extrapolated from layer y[1] to yc[0]
		uy[0], uy[-1] = .5/h[1] * (u[1]-u[0]), .5/h[-1] * (u[-1]-u[-2])
		wy[0], wy[-1] = .5/h[1] * (w[1]-w[0]), .5/h[-1] * (w[-1]-w[-2])
		# vx,vz are determined by linear extrapolation of v to yc[0]
		vx[0] = h[1]/dy[1] * (vx[1]-vx[2]) + .5 * (vx[1]+vx[2]) # another half of the boundary is left after the CC-interpolation to avoid interference
		vz[0] = h[1]/dy[1] * (vz[1]-vz[2]) + .5 * (vz[1]+vz[2])
		# ux,uz,wx,wz are automatically on layer yc[0]
		pass

		## interpolate cross derivatives to cell-centers
		uz[:] = .5 * (uz + uz[:,:,ip])
		wx[:] = .5 * (wx + wx[:,kp,:])
		uy[:] = .5 * (uy + uy[:,:,ip])
		wy[:] = .5 * (wy + wy[:,kp,:])
		vx[jc]= .5 * (vx[jc] + vx[jp])
		vz[jc]= .5 * (vz[jc] + vz[jp])

		# the boundary left to be processed
		vx[-1]= h[-1]/dy[-2] * (vx[-1]-vx[-2]) + .5 * (vx[-1]+vx[-2])
		vz[-1]= h[-1]/dy[-2] * (vz[-1]-vz[-2]) + .5 * (vz[-1]+vz[-2])


		# # interpolate cross derivatives to cell-centers
		# uz[:] = .5 * (uz + uz[:,:,ip])
		# wx[:] = .5 * (wx + wx[:,kp,:])
		# vx[jc]= .5 * (vx[jc] + vx[jp])
		# vz[jc]= .5 * (vz[jc] + vz[jp])
		# uy[jc]= .5 * (uy + uy[:,:,ip])[jc]
		# wy[jc]= .5 * (wy + wy[:,kp,:])[jc]
		# # handle the boundaries (think of boundary on yc[0])
		# # vy is extrapolated from layer yc[1] to yc[0]
		# vy[0], vy[-1] = vy[1], vy[-2]
		# # vx,vz,uy,wy are extrapolated from layer y[1] to yc[0]
		# vx[0], vx[-1] = .5/dx * (v[1][:,ip] - v[1][:,im]), .5/dx * (v[-1][:,ip] - v[-1][:,im])	# iterable objects can only be used as slice when theres no other indeces than :
		# vz[0], vz[-1] = .5/dz * (v[1][kp,:] - v[1][km,:]), .5/dz * (v[-1][kp,:] - v[-1][km,:])
		# uy[0], uy[-1] = .5/h[1] * ((u[1]-u[0]) + (u[1]-u[0])[:,ip]), .5/h[-1] * ((u[-1]-u[-2]) + (u[-1]-u[-2])[:,ip])
		# wy[0], wy[-1] = .5/h[1] * ((w[1]-w[0]) + (w[1]-w[0])[kp,:]), .5/h[-1] * ((w[-1]-w[-2]) + (w[-1]-w[-2])[kp,:])
		# # ux,uz,wx,wz are automatically on layer yc[0]
		# pass

		return ux,vx,wx, uy,vy,wy, uz,vz,wz






# # definition of fft functions where hfft/ihfft are not avilable
# # these are tested to be valid but about 1.5 less efficient than the hfft based ones
# def spec(q):
# 	return ifft2(q).T[:q.shape[-1]//2+1].T # .T usually returns a view which is efficient
# def phys(q):
# 	qc = np.conj(np.roll(q.T[-2:0:-1][:,::-1], 1, axis=1)) # Nx must be even
# 	return fft2(np.vstack([q.T, qc]).T).real


# def stg2cc(self, u, v, w):
# 	self.__to_cellcenter(u, stagtyp=1)
# 	self.__to_cellcenter(v, stagtyp=2)
# 	self.__to_cellcenter(w, stagtyp=3)
# 	# error: return what?

# def cc2stg(self, u, v, w):
# 	h, dy = self.para.h, self.para.dy
# 	u = .5 * (u + np.roll(u,1,axis=-1))
# 	w = .5 * (w + np.roll(w,1,axis=-2))
# 	v[1:] = .5/h[1:] * (v[:-1]*dy[1:] + v[1:]*dy[:-1])
# 	v[0] = 0
# 	# error: return what?




