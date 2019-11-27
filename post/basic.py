import numpy as np
from numpy.fft import fft2 as phys
from numpy.fft import ifft2 as spec
from scipy.integrate import trapz
from struct import pack, unpack
import os




##### common functions #####

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
		self.read_GRD(self.statpath)
		self.tsteps = self.get_tsteps(self.fieldpath)

		self.Nxz = self.Nx * self.Nz
		self.Nxc = int(self.Nx / 2 + 1)
		self.Nzc = int(self.Nz / 2 + 1)
		self.kx = np.array( list(range(self.Nxc)) + list(range(self.Nxc-self.Nx, 0)) ) * (2*np.pi/self.Lx)
		self.kz = np.array( list(range(self.Nzc)) + list(range(self.Nzc-self.Nz, 0)) ) * (2*np.pi/self.Lz)

	def read_XIN(self, path):
		with open(path + "XINDAT") as fp:
			for line in fp:
				if "//" in line: line = line[:line.index("//")]

				if "fieldpath" in line: self.fieldpath = path + line.split('=')[1].split('"')[1]
				if "probepath" in line: self.probepath = path + line.split('=')[1].split('"')[1]
				if "statpath" in line: self.statpath = path + line.split('=')[1].split('"')[1]
				if "postpath" in line: self.postpath = path + line.split('=')[1].split('"')[1]

				if "Re" in line: self.Re = float( line.split('=')[1].strip().split()[0] )

				if "Nx" in line: self.Nx = int( line.split('=')[1].strip().split()[0] )
				if "Ny" in line: self.Ny = int( line.split('=')[1].strip().split()[0] )
				if "Nz" in line: self.Nz = int( line.split('=')[1].strip().split()[0] )
				if "Lx" in line: self.Lx = float( line.split('=')[1].strip().split()[0] )
				if "Ly" in line: self.Ly = float( line.split('=')[1].strip().split()[0] )
				if "Lz" in line: self.Lz = float( line.split('=')[1].strip().split()[0] )

	def read_GRD(self, path):
		Ny = self.Ny
		self.y = np.loadtxt(path + "CHANNEL.GRD").reshape(Ny+1)
		self.yc = 0.5 * (
			self.y[ [1] + list(range(1,Ny+1)) ] +
			self.y[ list(range(1,Ny+1)) + [Ny] ]	)

	def get_tsteps(self, path):
		tsteps = []
		for name in [s for s in os.listdir(path) if ".bin" in s]:
			tstep = int( "".join([c for c in name if c.isdigit()]) )
			if tstep not in tsteps: tsteps.append(tstep)
		return np.sort(tsteps)



class Field:
	def __init__(self, para):
		self.Nx, self.Ny, self.Nz, self.Nxz = para.Nx, para.Ny, para.Nz, para.Nxz
		self.fieldpath = para.fieldpath

	def read_fluc_mean(self, name, stagtyp=4):
		q = self.__to_cellcenter(
			read_channel(self.fieldpath + name),
			stagtyp if stagtyp in (0,1,2,3) else self.__infer_stagtyp(name)	)
		qm = np.mean(q, axis=(-1,-2))
		return q - qm.reshape([self.Ny+1, 1, 1]), qm

	def read_mean(self, name, stagtyp=4):
		q = np.mean(read_channel(self.fieldpath + name), axis=(-1,-2))
		return q if stagtyp != 2 else \
			0.5 * ( q[[1] + list(range(1,Ny+1))] + q[list(range(1,Ny+1)) + [Ny]] )

	def read_fluc(self, name, stagtyp=4):
		return self.read_fluc_mean(name, stagtyp)[0]

	def __to_cellcenter(self, q, stagtyp):
		Nx, Ny, Nz = self.Nx, self.Ny, self.Nz
		if stagtyp == 1: q[:] = 0.5 * ( q + q[:,:,(np.arange(Nx)+1) % Nx] )
		if stagtyp == 3: q[:] = 0.5 * ( q + q[:,(np.arange(Nz)+1) % Nz] )
		if stagtyp == 2: q[:] = 0.5 * ( q[[1] + list(range(1,Ny+1))] + q[list(range(1,Ny+1)) + [Ny]] )
		if stagtyp not in (0,1,2,3): print("\nUnknown stagger type !\n")
		return q

	def __infer_stagtyp(self, name):
		return	1 if 'U' in name else \
				2 if 'V' in name else \
				3 if 'W' in name else 0


class Statis:
	def __init__(self, para, feld):
		self.para = para
		self.feld = feld

	def calc_umean(self):
		Ny = self.para.Ny
		self.Um = np.zeros(Ny+1)
		for tstep in self.para.tsteps:
			print("Reading umean: tstep", tstep)
			u = self.feld.read_mean("U%08i.bin"%tstep, stagtyp=1)
			self.Um += u / len(self.para.tsteps)

	def calc_statis(self):
		nx, ny, nz = self.para.Nx, self.para.Ny + 1, self.para.Nz
		nxc, nzc = int(nx/2+1), int(nz/2+1)

		self.R11, self.R22, self.R33 = [np.zeros(ny) for n in range(3)]
		self.R12, self.R23, self.R13 = [np.zeros(ny) for n in range(3)]
		self.Rpu, self.Rpv, self.Rpw, self.Rpp= [np.zeros(ny) for n in range(4)]
		self.Um , self.Vm , self.Wm , self.Pm = [np.zeros(ny) for n in range(4)]
		self.Euu, self.Evv, self.Eww, self.Epp= [np.zeros([ny, nzc, nxc]) for n in range(4)]
		self.Euv, self.Evw, self.Euw = [np.zeros([ny, nzc, nxc]) for n in range(3)]

		for tstep in self.para.tsteps:
			print("Reading statis: tstep", tstep)
			u, um = self.feld.read_fluc_mean("U%08i.bin"%tstep)
			v, vm = self.feld.read_fluc_mean("V%08i.bin"%tstep)
			w, wm = self.feld.read_fluc_mean("W%08i.bin"%tstep)
			p, pm = self.feld.read_fluc_mean("P%08i.bin"%tstep)

			self.Um += um / len(self.para.tsteps)
			self.Vm += vm / len(self.para.tsteps)
			self.Wm += wm / len(self.para.tsteps)
			self.Pm += pm / len(self.para.tsteps)

			self.R11 += np.mean(u**2,axis=(-1,-2)) / len(self.para.tsteps)
			self.R22 += np.mean(v**2,axis=(-1,-2)) / len(self.para.tsteps)
			self.R33 += np.mean(w**2,axis=(-1,-2)) / len(self.para.tsteps)
			self.R12 += np.mean(u*v, axis=(-1,-2)) / len(self.para.tsteps)
			self.R23 += np.mean(v*w, axis=(-1,-2)) / len(self.para.tsteps)
			self.R13 += np.mean(u*w, axis=(-1,-2)) / len(self.para.tsteps)

			self.Rpu += np.mean(p*u, axis=(-1,-2)) / len(self.para.tsteps)
			self.Rpv += np.mean(p*v, axis=(-1,-2)) / len(self.para.tsteps)
			self.Rpw += np.mean(p*w, axis=(-1,-2)) / len(self.para.tsteps)
			self.Rpp += np.mean(p**2,axis=(-1,-2)) / len(self.para.tsteps)

			self.Euu += self.__flipk( np.abs(spec(u))**2 ) / len(self.para.tsteps)
			self.Evv += self.__flipk( np.abs(spec(v))**2 ) / len(self.para.tsteps)
			self.Eww += self.__flipk( np.abs(spec(w))**2 ) / len(self.para.tsteps)
			self.Epp += self.__flipk( np.abs(spec(p))**2 ) / len(self.para.tsteps)
			self.Euv += self.__flipk( (np.conj(spec(u))*spec(v)).real ) / len(self.para.tsteps)
			self.Evw += self.__flipk( (np.conj(spec(v))*spec(w)).real ) / len(self.para.tsteps)
			self.Euw += self.__flipk( (np.conj(spec(u))*spec(w)).real ) / len(self.para.tsteps)

	def __flipk(self, q):
		ny, nz, nx = q.shape
		nxc, nzc = int(nx/2+1), int(nz/2+1)
		q[:, :, 1:nxc] += q[:, :, :-nxc:-1]
		q[:, 1:nzc, :] += q[:, :-nzc:-1, :]
		if not (nx % 2): q[:,:,nxc-1] /= 2
		if not (nz % 2): q[:,nzc-1,:] /= 2
		return q[:,:nzc,:nxc]

	def flipy(self):
		self.Um[:] = 0.5 * (self.Um + self.Um[::-1])
		self.Vm[:] = 0.5 * (self.Vm - self.Vm[::-1])
		self.Wm[:] = 0.5 * (self.Wm + self.Wm[::-1])
		self.Pm[:] = 0.5 * (self.Pm + self.Pm[::-1])

		self.R11[:] = 0.5 * (self.R11 + self.R11[::-1])
		self.R22[:] = 0.5 * (self.R22 + self.R22[::-1])
		self.R33[:] = 0.5 * (self.R33 + self.R33[::-1])
		self.R12[:] = 0.5 * (self.R12 - self.R12[::-1])
		self.R23[:] = 0.5 * (self.R23 - self.R23[::-1])
		self.R13[:] = 0.5 * (self.R13 + self.R13[::-1])

		self.Rpu[:] = 0.5 * (self.Rpu + self.Rpu[::-1])
		self.Rpv[:] = 0.5 * (self.Rpv - self.Rpv[::-1])
		self.Rpw[:] = 0.5 * (self.Rpw + self.Rpw[::-1])
		self.Rpp[:] = 0.5 * (self.Rpp + self.Rpp[::-1])

		self.Euu[:] = 0.5 * (self.Euu + self.Euu[::-1])
		self.Evv[:] = 0.5 * (self.Evv + self.Evv[::-1])
		self.Eww[:] = 0.5 * (self.Eww + self.Eww[::-1])
		self.Epp[:] = 0.5 * (self.Epp + self.Epp[::-1])
		self.Euv[:] = 0.5 * (self.Euv - self.Euv[::-1])
		self.Evw[:] = 0.5 * (self.Evw - self.Evw[::-1])
		self.Euw[:] = 0.5 * (self.Euw + self.Euw[::-1])

	def calc_wallscale(self):
		nudUdy = np.gradient(self.Um, self.para.yc) / self.para.Re

		# self.tauw = 0.5 * (nudUdy[0] - nudUdy[-1])
		# integrate to get tauw using relation dU^+/dy^+ - <u'v'>^+ = 1 - y\bar, accurancy tested to be good ( O(10e-4) )
		# self.tauw = np.sum( abs(nudUdy - self.R12)[1:-1] * (self.para.y[2:] - self.para.y[1:-1]) ) / (0.25 * self.para.Ly**2)
		self.tauw = trapz(abs(nudUdy - self.R12), self.para.yc) / (0.25 * self.para.Ly**2)

		self.utau = self.tauw**0.5
		self.dnu = 1.0 / self.para.Re / self.utau
		self.tnu  = self.dnu / self.utau
		self.Ret = 1.0 / self.dnu # channel height is taken for 2.0

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




class Budgets:
	def __init__(self, para):
		self.para = para

	def dissipation(self):
		nx, ny, nz = self.para.Nx, self.para.Ny + 1, self.para.Nz

		dx, dz = self.para.Lx / nx, self.para.Lz / nz
		dy= np.hstack([[0], self.para.y[2:] - self.para.y[1:-1], [0]])
		h = np.hstack([[0], .5 * (dy[1:] + dy[:-1])])
		dy, h = dy.reshape([ny,1,1]), h.reshape([ny,1,1])

		jc = np.arange(1,ny-1)
		jm, jp = jc-1, jc+1
		im, ip = np.arange(-1,nx-1), np.arange(1,nx+1) % nx
		km, kp = np.arange(-1,nz-1), np.arange(1,nz+1) % nz

		self.epsl = np.zeros(ny)
		ux, uy, uz = [np.zeros([ny,nz,nx]) for n in range(3)]
		vx, vy, vz = [np.zeros([ny,nz,nx]) for n in range(3)]
		wx, wy, wz = [np.zeros([ny,nz,nx]) for n in range(3)]

		for tstep in self.para.tsteps:
			print("Reading budgets: tstep", tstep)
			u = read_channel(self.para.fieldpath + "U%08i.bin"%tstep)
			v = read_channel(self.para.fieldpath + "V%08i.bin"%tstep)
			w = read_channel(self.para.fieldpath + "W%08i.bin"%tstep)

			u -= np.mean(u, axis=(-1,-2)).reshape([ny,1,1])
			v -= np.mean(v, axis=(-1,-2)).reshape([ny,1,1])
			w -= np.mean(w, axis=(-1,-2)).reshape([ny,1,1])

			ux[:] = (u[:,:,ip] - u) / dx
			wz[:] = (w[:,kp,:] - w) / dz
			uz[:] = .5/dz * (u[:,kp,:] - u[:,km,:])
			wx[:] = .5/dx * (w[:,:,ip] - w[:,:,im])
			vx[:] = .5/dx * (v[:,:,ip] - v[:,:,im])
			vz[:] = .5/dz * (v[:,kp,:] - v[:,km,:])
			vy[jc]= (v[jp] - v[jc]) / dy[jc]
			uy[jc]= .5/dy[jc] * ( (u[jc]*dy[jp]+u[jp]*dy[jc]) / h[jp] - (u[jc]*dy[jm]+u[jm]*dy[jc]) / h[jc] )
			wy[jc]= .5/dy[jc] * ( (w[jc]*dy[jp]+w[jp]*dy[jc]) / h[jp] - (w[jc]*dy[jm]+w[jm]*dy[jc]) / h[jc] )
			
			uz[:] = .5 * (uz + uz[:,:,ip])
			wx[:] = .5 * (wx + wx[:,kp,:])
			vx[jc]= .5 * (vx[jc] + vx[jp])
			vz[jc]= .5 * (vz[jc] + vz[jp])
			uy[jc]= .5 * (uy + uy[:,:,ip])[jc]
			wy[jc]= .5 * (wy + wy[:,kp,:])[jc]

			vy[0], vy[-1] = vy[1], vy[-2]
			vx[0], vx[-1] = .5/dx * (v[1,:,ip] - v[1,:,im]), .5/dx * (v[-1,:,ip] - v[-1,:,im])
			vz[0], vz[-1] = .5/dz * (v[1,kp,:] - v[1,km,:]), .5/dz * (v[-1,kp,:] - v[-1,km,:])
			uy[0], uy[-1] = .5/h[1] * ((u[1]-u[0]) + (u[1]-u[0])[:,ip]), .5/h[-1] * ((u[-1]-u[-2]) + (u[-1]-u[-2])[:,ip])
			wy[0], wy[-1] = .5/h[1] * ((w[1]-w[0]) + (w[1]-w[0])[kp,:]), .5/h[-1] * ((w[-1]-w[-2]) + (w[-1]-w[-2])[kp,:])

			self.epsl += 1. / len(self.para.tsteps) / self.para.Re * \
				np.mean(ux**2+uy**2+uz**2+vx**2+vy**2+vz**2+wx**2+wy**2+wz**2, axis=(-1,-2))

	def flipy(self):
		self.epsl[:] = .5 * (self.epsl + self.epsl[::-1])