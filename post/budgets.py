from basic import *
from operators import Operators


class Budgets:
	def __init__(self, para):
		self.para = para

	def dissipation(self, tsteps=None):
		fieldpath = self.para.fieldpath
		Re = self.para.Re
		Ny = self.para.Ny
		if tsteps is None: tsteps = self.para.tsteps

		opt = Operator(self.para)
		self.epsl = np.zeros(Ny+1)

		for tstep in tsteps:
			print("Reading budgets: tstep", tstep)
			u = read_channel(fieldpath + "U%08i.bin"%tstep)
			v = read_channel(fieldpath + "V%08i.bin"%tstep)
			w = read_channel(fieldpath + "W%08i.bin"%tstep)
			# compute fluctuations on their own grids
			u.T[:] -= np.mean(u, axis=(-1,-2))
			v.T[:] -= np.mean(v, axis=(-1,-2))
			w.T[:] -= np.mean(w, axis=(-1,-2))

			ux,vx,wx, uy,vy,wy, uz,vz,wz = opt.gradient(u,v,w)

			self.epsl += 1./Re * np.mean(
				ux**2 + uy**2 + uz**2 + \
				vx**2 + vy**2 + vz**2 + \
				wx**2 + wy**2 + wz**2, axis=(-1,-2)) / len(tsteps)

	def flipy(self):
		self.epsl[:] = .5 * (self.epsl + self.epsl[::-1])


	# def _dissipation(self, tsteps=None):
	# 	fieldpath = self.para.fieldpath
	# 	Re = self.para.Re
	# 	Nx = self.para.Nx
	# 	Ny = self.para.Ny
	# 	Nz = self.para.Nz
	# 	dx = self.para.Lx / Nx
	# 	dz = self.para.Lz / Nz
	# 	dy = self.para.dy.reshape([Ny+1,1,1])
	# 	h  = self.para.h.reshape([Ny+1,1,1])
	# 	if tsteps is None: tsteps = self.para.tsteps

	# 	self.epsl = np.zeros(Ny+1)
	# 	ux, vx, wx = [np.zeros([Ny+1,Nz,Nx]) for _ in range(3)]
	# 	uy, vy, wy = [np.zeros([Ny+1,Nz,Nx]) for _ in range(3)]
	# 	uz, vz, wz = [np.zeros([Ny+1,Nz,Nx]) for _ in range(3)]

	# 	jc = np.arange(1,Ny)
	# 	jm, jp = jc-1, jc+1
	# 	im, ip = np.arange(-1,Nx-1), np.arange(1,Nx+1) % Nx
	# 	km, kp = np.arange(-1,Nz-1), np.arange(1,Nz+1) % Nz

	# 	for tstep in tsteps:
	# 		print("Reading budgets: tstep", tstep)
	# 		u = read_channel(fieldpath + "U%08i.bin"%tstep)
	# 		v = read_channel(fieldpath + "V%08i.bin"%tstep)
	# 		w = read_channel(fieldpath + "W%08i.bin"%tstep)
	# 		# compute fluctuations on their own grids
	# 		u.T[:] -= np.mean(u, axis=(-1,-2))
	# 		v.T[:] -= np.mean(v, axis=(-1,-2))
	# 		w.T[:] -= np.mean(w, axis=(-1,-2))
	# 		# compute normal derivatives on cell-centers
	# 		ux[:] = (u[:,:,ip] - u) / dx
	# 		wz[:] = (w[:,kp,:] - w) / dz
	# 		vy[jc]= (v[jp] - v[jc]) / dy[jc]
	# 		# compute cross derivatives on their own grids
	# 		uz[:] = .5/dz * (u[:,kp,:] - u[:,km,:])
	# 		wx[:] = .5/dx * (w[:,:,ip] - w[:,:,im])
	# 		vx[:] = .5/dx * (v[:,:,ip] - v[:,:,im])
	# 		vz[:] = .5/dz * (v[:,kp,:] - v[:,km,:])
	# 		uy[jc]= .5/dy[jc] * ( (u[jc]*dy[jp]+u[jp]*dy[jc]) / h[jp] - (u[jc]*dy[jm]+u[jm]*dy[jc]) / h[jc] )
	# 		wy[jc]= .5/dy[jc] * ( (w[jc]*dy[jp]+w[jp]*dy[jc]) / h[jp] - (w[jc]*dy[jm]+w[jm]*dy[jc]) / h[jc] )
	# 		# interpolate cross derivatives to cell-centers
	# 		uz[:] = .5 * (uz + uz[:,:,ip])
	# 		wx[:] = .5 * (wx + wx[:,kp,:])
	# 		vx[jc]= .5 * (vx[jc] + vx[jp])
	# 		vz[jc]= .5 * (vz[jc] + vz[jp])
	# 		uy[jc]= .5 * (uy + uy[:,:,ip])[jc]
	# 		wy[jc]= .5 * (wy + wy[:,kp,:])[jc]
	# 		# handle the boundaries
	# 		vy[0], vy[-1] = vy[1], vy[-2]
	# 		vx[0], vx[-1] = .5/dx * (v[1][:,ip] - v[1][:,im]), .5/dx * (v[-1][:,ip] - v[-1][:,im])	# iterable objects can only be used as slice when theres no other indeces than :
	# 		vz[0], vz[-1] = .5/dz * (v[1][kp,:] - v[1][km,:]), .5/dz * (v[-1][kp,:] - v[-1][km,:])
	# 		uy[0], uy[-1] = .5/h[1] * ((u[1]-u[0]) + (u[1]-u[0])[:,ip]), .5/h[-1] * ((u[-1]-u[-2]) + (u[-1]-u[-2])[:,ip])
	# 		wy[0], wy[-1] = .5/h[1] * ((w[1]-w[0]) + (w[1]-w[0])[kp,:]), .5/h[-1] * ((w[-1]-w[-2]) + (w[-1]-w[-2])[kp,:])

	# 		self.epsl += 1./Re * np.mean(
	# 			ux**2 + uy**2 + uz**2 + \
	# 			vx**2 + vy**2 + vz**2 + \
	# 			wx**2 + wy**2 + wz**2, axis=(-1,-2)) / len(tsteps)





	