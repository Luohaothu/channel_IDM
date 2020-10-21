from basic import *


class Operators:
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
		dx = self.para.dx
		dz = self.para.dz.reshape([-1,1])
		dy = self.para.dy.reshape([-1,1,1])

		ux, wz, vy = (np.empty([Ny+1, Nz+1, Nx+1]) for _ in range(3))

		ux[:,:,1:-1] = (u[:,:,2:] - u[:,:,1:-1]) / dx[1:-1]
		wz[:,1:-1]   = (w[:,2:]   - w[:,1:-1]  ) / dz[1:-1]
		vy[1:-1]     = (v[2:]     - v[1:-1]    ) / dy[1:-1]
		
		ux[:,:,[0,-1]] = ux[:,:,[-2,1]] # periodic boundary
		wz[:,[0,-1]]   = wz[:,[-2,1]]   # periodic boudnary
		vy[[0,-1]]     = vy[[1,-2]]     # linear extrapolated v

		return ux + vy + wz

	def laplace(self, q):
		''' compute Laplacian of cell-centered q to cell centers
			x, z boundaries are periodic, y boundary linear extrapolation for q '''
		Nx = self.para.Nx
		Ny = self.para.Ny
		Nz = self.para.Nz
		dx = self.para.dx
		dz = self.para.dz.reshape([-1,1])
		dy = self.para.dy.reshape([-1,1,1])
		hx = self.para.hx
		hz = self.para.hz.reshape([-1,1])
		hy = self.para.hy.reshape([-1,1,1])
		
		qxx, qyy, qzz = (np.empty([Ny+1, Nz+1, Nx+1]) for _ in range(3))

		im, ic, ip = np.arange(Nx-1), np.arange(1,Nx), np.arange(2,Nx+1)
		km, kc, kp = np.arange(Nz-1), np.arange(1,Nz), np.arange(2,Nz+1)
		jm, jc, jp = np.arange(Ny-1), np.arange(1,Ny), np.arange(2,Ny+1)

		qxx[:,:,ic] = ((q[:,:,ip] - q[:,:,ic]) / hx[ip] - (q[:,:,ic] - q[:,:,im]) / hx[ic]) / dx[ic]
		qzz[:,kc]   = ((q[:,kp]   - q[:,kc]  ) / hz[kp] - (q[:,kc]   - q[:,km]  ) / hz[kc]) / dz[kc]
		qyy[jc]     = ((q[jp]     - q[jc]    ) / hy[jp] - (q[jc]     - q[jm]    ) / hy[jc]) / dy[jc]
		
		qxx[:,:,[0,-1]] = qxx[:,:,[-2,1]] # periodic boundary
		qzz[:,[0,-1]]   = qzz[:,[-2,1]]   # periodic boudnary
		qyy[[0,-1]]     = 0               # linear extrapolated q

		return qxx + qyy + qzz

	def __gradient_scalar(self, q):
		''' input on cell-centers and output on staggered grids '''
		Nx = self.para.Nx
		Ny = self.para.Ny
		Nz = self.para.Nz
		hx = self.para.hx
		hz = self.para.hz.reshape([-1,1])
		hy = self.para.hy.reshape([-1,1,1])

		qx, qy, qz = (np.zeros([Ny+1, Nz+1, Nx+1]) for _ in range(3))

		qx[:,:,1:] = (q[:,:,1:] - q[:,:,:-1]) / hx[1:]
		qz[:,1:]   = (q[:,1:]   - q[:,:-1]  ) / hz[1:]
		qy[1:]     = (q[1:]     - q[:-1]    ) / hy[1:]

		return qx, qy, qz

	def __gradient_vector(self, u, v, w):
		''' input on staggered grids and output on cell-centers '''
		Nx = self.para.Nx
		Ny = self.para.Ny
		Nz = self.para.Nz
		dx = self.para.dx
		dz = self.para.dz.reshape([-1,1])
		dy = self.para.dy.reshape([-1,1,1])
		hx = self.para.hx
		hz = self.para.hz.reshape([-1,1])
		hy = self.para.hy.reshape([-1,1,1])

		ux, vx, wx = (np.empty([Ny+1, Nz+1, Nx+1]) for _ in range(3))
		uy, vy, wy = (np.empty([Ny+1, Nz+1, Nx+1]) for _ in range(3))
		uz, vz, wz = (np.empty([Ny+1, Nz+1, Nx+1]) for _ in range(3))

		im, ic, ip = np.arange(Nx-1), np.arange(1,Nx), np.arange(2,Nx+1)
		km, kc, kp = np.arange(Nz-1), np.arange(1,Nz), np.arange(2,Nz+1)
		jm, jc, jp = np.arange(Ny-1), np.arange(1,Ny), np.arange(2,Ny+1)

		## compute the derivatives
		# normal derivatives on cell-centers
		ux[:,:,ic] = (u[:,:,ip] - u[:,:,ic]) / dx[ic]
		wz[:,kc]   = (w[:,kp]   - w[:,kc]  ) / dz[kc]
		vy[jc]     = (v[jp]     - v[jc]    ) / dy[jc]
		# cross derivatives on staggered grids
		vx[:,:,ic] = .5/dx[ic] * ((v[:,:,ic]*dx[ip] + v[:,:,ip]*dx[ic]) / hx[ip] - (v[:,:,ic]*dx[im] + v[:,:,im]*dx[ic]) / hx[ic])
		wx[:,:,ic] = .5/dx[ic] * ((w[:,:,ic]*dx[ip] + w[:,:,ip]*dx[ic]) / hx[ip] - (w[:,:,ic]*dx[im] + w[:,:,im]*dx[ic]) / hx[ic])
		uz[:,kc]   = .5/dz[kc] * ((u[:,kc]  *dz[kp] + u[:,kp]  *dz[kc]) / hz[kp] - (u[:,kc]  *dz[km] + u[:,km]  *dz[kc]) / hz[kc])
		vz[:,kc]   = .5/dz[kc] * ((v[:,kc]  *dz[kp] + v[:,kp]  *dz[kc]) / hz[kp] - (v[:,kc]  *dz[km] + v[:,km]  *dz[kc]) / hz[kc])
		uy[jc]     = .5/dy[jc] * ((u[jc]    *dy[jp] + u[jp]    *dy[jc]) / hy[jp] - (u[jc]    *dy[jm] + u[jm]    *dy[jc]) / hy[jc])
		wy[jc]     = .5/dy[jc] * ((w[jc]    *dy[jp] + w[jp]    *dy[jc]) / hy[jp] - (w[jc]    *dy[jm] + w[jm]    *dy[jc]) / hy[jc])

		## y boundary process that need to be done before cell-center interpolation (think of boundary on yc[0])
		# uy & wy on yc[0] taken from y[1]
		uy[[0,-1]] = .5/hy[1] * (u[1]-u[0]), .5/hy[-1] * (u[-1]-u[-2])
		wy[[0,-1]] = .5/hy[1] * (w[1]-w[0]), .5/hy[-1] * (w[-1]-w[-2])
		# linear extrapolation for boundary v
		bvx = hy[1]/dy[1] * (vx[1]-vx[2]) + .5 * (vx[1]+vx[2]), hy[-1]/dy[-2] * (vx[-1]-vx[-2]) + .5 * (vx[-1]+vx[-2])
		bvz = hy[1]/dy[1] * (vz[1]-vz[2]) + .5 * (vz[1]+vz[2]), hy[-1]/dy[-2] * (vz[-1]-vz[-2]) + .5 * (vz[-1]+vz[-2])

		## interpolate cross derivatives to cell-centers
		uz[:,:,ic] = .5 * (uz[:,:,ic] + uz[:,:,ip])
		uy[:,:,ic] = .5 * (uy[:,:,ic] + uy[:,:,ip])
		wx[:,kc]   = .5 * (wx[:,kc]   + wx[:,kp])
		wy[:,kc]   = .5 * (wy[:,kc]   + wy[:,kp])
		vx[jc]     = .5 * (vx[jc]     + vx[jp])
		vz[jc]     = .5 * (vz[jc]     + vz[jp])

		## handle the boundaries
		vx[[0,-1]] = bvx
		vz[[0,-1]] = bvz
		vy[[0,-1]] = vy[[1,-2]]

		ux[:,:,[0,-1]] = ux[:,:,[-2,1]]
		vx[:,:,[0,-1]] = vx[:,:,[-2,1]]
		wx[:,:,[0,-1]] = wx[:,:,[-2,1]]
		uy[:,:,[0,-1]] = uy[:,:,[-2,1]]
		uz[:,:,[0,-1]] = uz[:,:,[-2,1]]

		uz[:,[0,-1]] = uz[:,[-2,1]]
		vz[:,[0,-1]] = vz[:,[-2,1]]
		wz[:,[0,-1]] = wz[:,[-2,1]]
		wx[:,[0,-1]] = wx[:,[-2,1]]
		wy[:,[0,-1]] = wy[:,[-2,1]]

		return ux,vx,wx, uy,vy,wy, uz,vz,wz

