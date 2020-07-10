from basic import *


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

		logs = np.loadtxt(self.para.statpath+'LOG.dat', skiprows=3)
		self.tauw = - np.mean([log[6] for log in logs if int(log[0]) in tsteps])

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
		    
		# Ly = self.para.Ly
		Re = self.para.Re
		# yc = self.para.yc
		# nudUdy = np.gradient(self.Um, yc) / Re
		# self.tauw = 4./Ly**2 * trapz(abs(nudUdy-self.R12), yc)

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

