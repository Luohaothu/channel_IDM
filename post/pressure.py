#!/root/Software/anaconda3/bin/python3
from basic import *
from operators import Operators
from tools import Tools


class Pressure:
	def __init__(self, para):
		self.para = para

	def get_Q(self, u, v, w):
		''' compute Q = -1/2 * \nabla \dot N where N(u) = \nabla \dot (uu) '''
		Nx = self.para.Nx
		Ny = self.para.Ny
		Nz = self.para.Nz
		dx = self.para.dx
		dz = self.para.dz.reshape([-1,1])
		dy = self.para.dy.reshape([-1,1,1])
		hx = self.para.hx
		hz = self.para.hz.reshape([-1,1])
		hy = self.para.hy.reshape([-1,1,1])

		uu, vv, ww = (np.empty([Ny+1, Nz+1, Nx+1]) for _ in range(3))
		uv, vw, wu = (np.empty([Ny+1, Nz+1, Nx+1]) for _ in range(3))
		N1, N2, N3 = (np.zeros([Ny+1, Nz+1, Nx+1]) for _ in range(3))

		im, ic, ip = np.arange(Nx-1), np.arange(1,Nx), np.arange(2,Nx+1)
		km, kc, kp = np.arange(Nz-1), np.arange(1,Nz), np.arange(2,Nz+1)
		jm, jc, jp = np.arange(Ny-1), np.arange(1,Ny), np.arange(2,Ny+1)

		## compute 2nd order terms

		# normal terms on cell centers
		uu[:,:,ic]  = .25 * (u[:,:,ic] + u[:,:,ip])**2 # on cell centers
		ww[:,kc]    = .25 * (w[:,kc]   + w[:,kp]  )**2 # on cell centers
		vv[jc]      = .25 * (v[jc]     + v[jp]    )**2 # on cell centers
		# boundaries of normal terms
		uu[:,:,[0,-1]] = uu[:,:,[-2,1]] # periodic u boundary
		ww[:,[0,-1]]   = ww[:,[-2,1]]   # periodic w boundary
		vv[[0,-1]]     = (hy[1]/dy[1] * (v[1]-v[2]) + .5 * (v[1]+v[2]))**2, \
					 (hy[-1]/dy[-2] * (v[-1]-v[-2]) + .5 * (v[-1]+v[-2]))**2 # linear extrapolation of v
		# cross terms on edges
		uv[1:]      = .5/hy[1:] * (u[1:]    *dy[:-1] + u[:-1]    *dy[1:]) # on z-edges
		uv[:,:,1:] *= .5/hx[1:] * (v[:,:,1:]*dx[:-1] + v[:,:,:-1]*dx[1:])
		vw[:,1:]    = .5/hz[1:] * (v[:,1:]  *dz[:-1] + v[:,:-1]  *dz[1:]) # on x-edges
		vw[1:]     *= .5/hy[1:] * (w[1:]    *dy[:-1] + w[:-1]    *dy[1:])
		wu[:,:,1:]  = .5/hx[1:] * (w[:,:,1:]*dx[:-1] + w[:,:,:-1]*dx[1:]) # on y-edges
		wu[:,1:]   *= .5/hz[1:] * (u[:,1:]  *dz[:-1] + u[:,:-1]  *dz[1:])

		## compute the convection operator N(u) = \nabla \dot (uu)

		# 3 components of N on staggered grids
		N1[:,:,1:] += (uu[:,:,1:] - uu[:,:,:-1]) / hx[1:]
		N1[:,kc]   += (wu[:,kp]   - wu[:,kc]   ) / dz[kc]
		N1[jc]     += (uv[jp]     - uv[jc]     ) / dy[jc]

		N2[:,:,ic] += (uv[:,:,ip] - uv[:,:,ic]) / dx[ic]
		N2[:,kc]   += (vw[:,kp]   - vw[:,kc]  ) / dz[kc]
		N2[1:]     += (vv[1:]     - vv[:-1]   ) / hy[1:]

		N3[:,:,ic] += (wu[:,:,ip] - wu[:,:,ic]) / dx[ic]
		N3[:,1:]   += (ww[:,1:]   - ww[:,:-1] ) / hz[1:]
		N3[jc]     += (vw[jp]     - vw[jc]    ) / dy[jc]
		# boundaries of N
		# N1, N3 on yc[0] are automatically zero
		N1[:,[0,-1]]   = N1[:,[-2,1]]
		N2[:,[0,-1]]   = N2[:,[-2,1]]
		N2[:,:,[0,-1]] = N2[:,:,[-2,1]]
		N3[:,:,[0,-1]] = N3[:,:,[-2,1]]

		## Q = -1/2 * \nabla \dot N, the returned result's boundary not handled
		return -.5 * Operators(self.para).divergence(N1,N2,N3)

		##### version1: basic (suffer to unacceptable error)
		# ux,vx,wx, uy,vy,wy, uz,vz,wz = Operators(self.para).gradient(u,v,w)
		# return -.5 * (ux**2+vy**2+wz**2 + 2*(uy*vx+uz*wx+wy*vz))

	def rhs_rapid_slow(self, u, v, w):
		''' note: rapid term is more recommended to solve by get_Q(u)-get_Q(u') '''
		Nx = self.para.Nx
		Ny = self.para.Ny
		Nz = self.para.Nz
		dx = self.para.dx
		dz = self.para.dz.reshape([-1,1])
		dy = self.para.dy.reshape([-1,1,1])
		hx = self.para.hx
		hz = self.para.hz.reshape([-1,1])
		hy = self.para.hy.reshape([-1,1,1])

		im, ic, ip = np.arange(Nx-1), np.arange(1,Nx), np.arange(2,Nx+1)
		km, kc, kp = np.arange(Nz-1), np.arange(1,Nz), np.arange(2,Nz+1)
		jm, jc, jp = np.arange(Ny-1), np.arange(1,Ny), np.arange(2,Ny+1)

		um = np.mean(u[:,1:-1,1:-1], axis=(-1,-2)).reshape([-1,1,1])
		vm = np.mean(v[:,1:-1,1:-1], axis=(-1,-2)).reshape([-1,1,1])
		wm = np.mean(w[:,1:-1,1:-1], axis=(-1,-2)).reshape([-1,1,1])

		uf = u - um
		vf = v - vm
		wf = w - wm

		# dUm/dy on cell centers
		Uy = np.empty([Ny+1,1,1])
		Uy[jc] = .5/dy[jc] * (
			(um[jc]*dy[jp] + um[jp]*dy[jc]) / hy[jp] - \
			(um[jc]*dy[jm] + um[jm]*dy[jc]) / hy[jc] )
		# use value on y[1] to approximate yc[0]
		Uy[[0,-1]] = .5/hy[1] * (um[1]-um[0]), .5/hy[-1]* (um[-1]-um[-2])

		# dv'/dx on cell centers
		vx = np.zeros([Ny+1, Nz+1, Nx+1])
		vx[:,:,ic] = .5/dx[ic] * ((vf[:,:,ic]*dx[ip] + vf[:,:,ip]*dx[ic]) / hx[ip] - (vf[:,:,ic]*dx[im] + vf[:,:,im]*dx[ic]) / hx[ic])
		bvx = hy[1]/dy[1] * (vx[1]-vx[2]) + .5 * (vx[1]+vx[2]), hy[-1]/dy[-2] * (vx[-1]-vx[-2]) + .5 * (vx[-1]+vx[-2]) # linear extrapolation of v
		vx[jc]= .5 * (vx[jc] + vx[jp])
		vx[[0,-1]] = bvx

		return -2.*Uy*vx, 2.*self.get_Q(uf,vf,wf)

	def poisson(self, rhs, bcoef=2):
		''' solve the poisson equation at cell-centers with specified BC
		    the boundaries of rhs are BC and the solution includes boundaries '''
		Nx = self.para.Nx
		Ny = self.para.Ny
		Nz = self.para.Nz
		Nxc= self.para.Nxc
		dx = self.para.dx
		dz = self.para.dz
		dy = self.para.dy
		hy = self.para.hy

		a = np.empty(Ny+1)
		b = np.empty(Ny+1)
		c = np.empty(Ny+1)
		
		jc = np.arange(1,Ny)
		jp = np.arange(2,Ny+1)

		# differencing coefficients
		a[jc] =  1./dy[jc] / hy[jc]
		b[jc] = -1./dy[jc] * (1./hy[jc] + 1./hy[jp])
		c[jc] =  1./dy[jc] / hy[jp]
		# boundary coefficients (a[0] = c[-1] = 0 are not relevant)
		if   bcoef is 2: a[-1], b[[0,-1]], c[0] = -1./hy[-1], (-1./hy[1], 1./hy[-1]), 1./hy[1]
		elif bcoef is 1: a[-1], b[[0,-1]], c[0] = 0, 1, 0
		else:            a[-1], b[[0,-1]], c[0] = bcoef

		# # assemble coefficients into matrix
		# mat          = np.diag(b)
		# mat[1:,:-1] += np.diag(a[1:])
		# mat[:-1,1:] += np.diag(c[:-1])

		# solve Poisson equation in spectral space
		phi = spec(rhs[:,1:-1,1:-1])
		kxs = 2./dx[1]**2 * (1 - np.cos(self.para.kx*dx[1])) # kx**2
		kzs = 2./dz[1]**2 * (1 - np.cos(self.para.kz*dz[1])) # kz**2

		for k, kz2 in enumerate(kzs):
			if not k%10: print('%i%%'%(100.*k/Nz))
			for i, kx2 in enumerate(kxs[:Nxc]):

				a_, b_, c_ = a, b.copy(), c
				b_[jc] -= kx2 + kz2

				if bcoef is 2 and not (i or k): b_[1] = 1e20

				phi[:,k,i] = Tools.chasing3(a_, b_, c_, phi[:,k,i])

				# mat[jc,jc] = b[jc] - (kx2 + kz2)
				# if bcoef is 2 and not (i or k): mat[1,1] = 1e20
				# phi[:,k,i] = np.linalg.solve(mat, phi[:,k,i]) # boundaries are also obtained

		p = np.empty([Ny+1, Nz+1, Nx+1])
		p[:,1:-1,1:-1] = phys(phi)
		p[:,:,[0,-1]]  = p[:,:,[-2,1]]
		p[:,[0,-1]]    = p[:,[-2,1]]

		return p

	def helmholtz(self, u, v, w):
		''' break a vector into solenoidal & dilatational parts '''
		Nx = self.para.Nx
		Ny = self.para.Ny
		Nz = self.para.Nz
		dx = self.para.dx
		dz = self.para.dz.reshape([-1,1])
		dy = self.para.dy.reshape([-1,1,1])
		hx = self.para.hx
		hz = self.para.hz.reshape([-1,1])
		hy = self.para.hy.reshape([-1,1,1])
		dy = self.para.dy.reshape([-1,1,1])
		hy = self.para.hy.reshape([-1,1,1])

		opt = Operators(self.para)
		
		# compute divergence and rotation of the vector field
		phi = opt.divergence(u,v,w)
		xi1, xi2, xi3 = opt.rotation(u,v,w)
		# implement homogeneous BC
		phi[[0,-1]] = 0
		xi1[[0,-1]] = 0
		xi2[[0,-1]] = 0
		xi3[[0,-1]] = 0

		# solve poisson equations with Dirichlet BC for xi2 and Neumann BC for others
		phi = self.poisson(phi)
		xi1 = self.poisson(-xi1)
		xi2 = self.poisson(-xi2, bcoef=1)
		xi3 = self.poisson(-xi3)

		# interpolate xi to staggered grids
		xi1[:,:,1:] = .5/hx[1:] * (xi1[:,:,1:]*dx[:-1] + xi1[:,:,:-1]*dx[1:])
		xi3[:,1:]   = .5/hz[1:] * (xi3[:,1:]  *dz[:-1] + xi3[:,:-1]  *dz[1:])
		xi2[1:]     = .5/hy[1:] * (xi2[1:]    *dy[:-1] + xi2[:-1]    *dy[1:])
		
		xi1[:,:,0] = 0
		xi3[:,0]   = 0
		xi2[0]     = 0

		# return Us (centered), Ud (staggered), xi (staggered), phi (centered)
		return opt.rotation(xi1,xi2,xi3), opt.gradient(phi), (xi1,xi2,xi3), phi







