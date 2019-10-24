import numpy as np
from numpy.fft import fft2 as phys
from numpy.fft import ifft2 as spec
from struct import pack, unpack
import os


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
		self.para = para

	def read_mean(self, name, stagtyp=4):
		q = np.mean(self.read_channel(self.para.fieldpath + name), axis=(-1,-2))
		return q if stagtyp != 2 else \
			0.5 * ( q[[1] + list(range(1,Ny+1))] + q[list(range(1,Ny+1)) + [Ny]] )

	def read_fluc_mean(self, name, stagtyp=4):
		q = self.__to_cellcenter(
			self.read_channel(self.para.fieldpath + name),
			stagtyp if stagtyp in (0,1,2,3) else self.__infer_stagtyp(name)	)
		qm = np.mean(q, axis=(-1,-2))
		return q - qm.reshape([self.para.Ny+1, 1, 1]), qm

	def read_fluc(self, name, stagtyp=4):
		return self.read_fluc_mean(name, stagtyp)[0]

	def read_channel(self, pame):
		Nx, Ny, Nz, Nxz = self.para.Nx, self.para.Ny, self.para.Nz, self.para.Nxz

		N = (Ny+1) * Nxz;
		with open(pame, 'rb') as fp:
			fp.seek(Nxz*8) # skipping info section
			q = np.reshape( unpack(N*'d', fp.read(N*8)), [Ny+1, Nz, Nx] )
		return q

	def __to_cellcenter(self, q, stagtyp):
		Nx, Ny, Nz = self.para.Nx, self.para.Ny, self.para.Nz
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
			self.Um += u
		self.Um[:] = 0.5 * ( self.Um + self.Um[::-1] ) / len(self.para.tsteps)

	def inner_scale(self):
		dU = self.Um[1] - self.Um[0]
		dy = self.para.yc[1] - self.para.yc[0]

		self.tauw = dU / dy / self.para.Re
		self.utau = self.tauw**0.5
		self.dnu = 1.0 / self.para.Re / self.utau
		self.tnu  = self.dnu / self.utau
		self.Ret = 0.5 * self.para.Ly / self.dnu

	def calc_profs(self):
		Ny = self.para.Ny
		self.R11, self.R22, self.R33 = [np.zeros(Ny+1) for n in range(3)]
		self.R12, self.R23, self.R13 = [np.zeros(Ny+1) for n in range(3)]
		self.Rpu, self.Rpv, self.Rpw, self.Rpp= [np.zeros(Ny+1) for n in range(4)]
		self.Um , self.Vm , self.Wm , self.Pm = [np.zeros(Ny+1) for n in range(4)]

		for tstep in self.para.tsteps:
			print("Reading profs: tstep", tstep)
			u, um = self.feld.read_fluc_mean("U%08i.bin"%tstep)
			v, vm = self.feld.read_fluc_mean("V%08i.bin"%tstep)
			w, wm = self.feld.read_fluc_mean("W%08i.bin"%tstep)
			p, pm = self.feld.read_fluc_mean("P%08i.bin"%tstep)

			self.Um += um
			self.Vm += um
			self.Wm += um
			self.Pm += um

			self.R11 += np.mean(u**2, axis=(-1,-2))
			self.R22 += np.mean(v**2, axis=(-1,-2))
			self.R33 += np.mean(w**2, axis=(-1,-2))
			self.R12 += np.mean(u*v, axis=(-1,-2))
			self.R23 += np.mean(v*w, axis=(-1,-2))
			self.R13 += np.mean(u*w, axis=(-1,-2))

			self.Rpu += np.mean(p*u, axis=(-1,-2))
			self.Rpv += np.mean(p*v, axis=(-1,-2))
			self.Rpw += np.mean(p*w, axis=(-1,-2))
			self.Rpp += np.mean(p**2, axis=(-1,-2))

		self.Um[:] = 0.5 * ( self.Um + self.Um[::-1] ) / len(self.para.tsteps)
		self.Vm[:] = 0.5 * ( self.Vm - self.Vm[::-1] ) / len(self.para.tsteps)
		self.Wm[:] = 0.5 * ( self.Wm + self.Wm[::-1] ) / len(self.para.tsteps)
		self.Pm[:] = 0.5 * ( self.Pm + self.Pm[::-1] ) / len(self.para.tsteps)

		self.R11[:] = 0.5 * ( self.R11 + self.R11[::-1] ) / len(self.para.tsteps)
		self.R22[:] = 0.5 * ( self.R22 + self.R22[::-1] ) / len(self.para.tsteps)
		self.R33[:] = 0.5 * ( self.R33 + self.R33[::-1] ) / len(self.para.tsteps)
		self.R12[:] = 0.5 * ( self.R12 - self.R12[::-1] ) / len(self.para.tsteps)
		self.R23[:] = 0.5 * ( self.R23 - self.R23[::-1] ) / len(self.para.tsteps)
		self.R13[:] = 0.5 * ( self.R13 + self.R13[::-1] ) / len(self.para.tsteps)

		self.Rpu[:] = 0.5 * ( self.Rpu + self.Rpu[::-1] ) / len(self.para.tsteps)
		self.Rpv[:] = 0.5 * ( self.Rpv - self.Rpv[::-1] ) / len(self.para.tsteps)
		self.Rpw[:] = 0.5 * ( self.Rpw + self.Rpw[::-1] ) / len(self.para.tsteps)
		self.Rpp[:] = 0.5 * ( self.Rpp + self.Rpp[::-1] ) / len(self.para.tsteps)


