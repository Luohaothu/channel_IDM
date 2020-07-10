from basic import *


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

