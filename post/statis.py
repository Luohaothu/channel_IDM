from basic import *


class Statis:
	def __init__(self, para):
		self.para = para
		self.feld = Field(para)

	def calc_umean(self, tsteps=None):
		self.Um = np.zeros(self.para.Ny+1)
		if tsteps is None: tsteps = self.para.tsteps
		for tstep in tsteps:
			print("Reading umean: tstep", tstep)
			u = self.feld.read_mean("U%08i.bin"%tstep, stagtyp=1)
			self.Um += u / len(tsteps)

	def calc_statis(self, tsteps=None):
		Ny = self.para.Ny
		Nz = self.para.Nz
		Nxc= self.para.Nxc
		if tsteps is None: tsteps = self.para.tsteps

		self.R11, self.R22, self.R33          = (np.zeros(Ny+1) for _ in range(3))
		self.R12, self.R23, self.R13          = (np.zeros(Ny+1) for _ in range(3))
		self.Rpu, self.Rpv, self.Rpw, self.Rpp= (np.zeros(Ny+1) for _ in range(4))
		self.Um , self.Vm , self.Wm , self.Pm = (np.zeros(Ny+1) for _ in range(4))
		self.Euu, self.Evv, self.Eww, self.Epp= (np.zeros([Ny+1, Nz-1, Nxc]) for _ in range(4))
		self.Euv, self.Evw, self.Euw          = (np.zeros([Ny+1, Nz-1, Nxc], dtype=complex) for _ in range(3))

		for tstep in tsteps:
			print("Reading statis: tstep", tstep)
			u, um = self.feld.read_fluc_mean("U%08i.bin"%tstep)
			v, vm = self.feld.read_fluc_mean("V%08i.bin"%tstep)
			w, wm = self.feld.read_fluc_mean("W%08i.bin"%tstep)
			p, pm = self.feld.read_fluc_mean("P%08i.bin"%tstep)

			fu = spec(u)
			fv = spec(v)
			fw = spec(w)
			fp = spec(p)

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

			self.Euu += np.abs(fu)**2 / len(tsteps)
			self.Evv += np.abs(fv)**2 / len(tsteps)
			self.Eww += np.abs(fw)**2 / len(tsteps)
			self.Epp += np.abs(fp)**2 / len(tsteps)
			self.Euv += fu.conj()*fv  / len(tsteps)
			self.Evw += fv.conj()*fw  / len(tsteps)
			self.Euw += fu.conj()*fw  / len(tsteps)

	@staticmethod
	def __flipk(q):
		''' fold all energy to the [:nzc,:nxc] range
		    Nx must be even, as required by hft, Nz can be even or odd  '''
		nzcu = int(np.ceil(q.shape[-2]/2))
		nzcd = q.shape[-2]//2
		p = np.copy((q.T[:,:nzcd+1]).T)
		p.T[:,1:nzcu] += q.T[:,:nzcd:-1]
		p.T[1:-1] *= 2
		return p

	def flipk(self):
		self.Euu = self.__flipk(self.Euu)
		self.Evv = self.__flipk(self.Evv)
		self.Eww = self.__flipk(self.Eww)
		self.Epp = self.__flipk(self.Epp)
		self.Euv = self.__flipk(self.Euv.real)
		self.Evw = self.__flipk(self.Evw.real)
		self.Euw = self.__flipk(self.Euw.real)

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

