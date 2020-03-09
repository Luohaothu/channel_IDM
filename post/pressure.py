#!/root/Software/anaconda3/bin/python3
from basic import *


class Pressure:
	def __init__(self, para):
		self.para = para

	def get_Q(self, u, v, w):
		Nx = self.para.Nx
		Ny = self.para.Ny
		Nz = self.para.Nz
		dx = self.para.Lx / Nx
		dz = self.para.Lz / Nz
		dy = self.para.dy.reshape([Ny+1,1,1])
		h  = self.para.h. reshape([Ny+1,1,1])

		us, vc, ws = (np.empty([Ny+1,Nz,Nx]) for _ in range(3))
		N1, N2, N3 = (np.empty([Ny+1,Nz,Nx]) for _ in range(3))

		jc = np.arange(1,Ny)
		jm, jp = jc-1, jc+1
		im, ip = np.arange(-1,Nx-1), np.arange(1,Nx+1) % Nx
		km, kp = np.arange(-1,Nz-1), np.arange(1,Nz+1) % Nz

		## compute the convection operator N(u) = \nabla(u \dot u)
		# y-interpolation of velocities (with boundary handled)
		us[1:] = .5/h[1:] * (u[:-1]*dy[1:] + u[1:]*dy[:-1]) # ZE
		ws[1:] = .5/h[1:] * (w[:-1]*dy[1:] + w[1:]*dy[:-1]) # XE
		vc[jc] = .5 * (v[jc] + v[jp]) # CC
		us[0] = 0
		ws[0] = 0
		vc[0] = h[ 1]/dy[ 1] * (v[ 1]-v[ 2]) + .5 * (v[ 1]+v[ 2]) # linear extrapolation of v
		vc[-1]= h[-1]/dy[-2] * (v[-1]-v[-2]) + .5 * (v[-1]+v[-2])
		# cross products on edges or self products on cell-centers
		vv = vc**2 # CC (0~Ny)
		uu = .25 * (u + u[:,:,ip])**2 # CC (0~Ny)
		ww = .25 * (w + w[:,kp,:])**2 # CC (0~Ny)
		uw = .25 * (u + u[:,km,:]) * (w + w[:,:,im]) # YE (0~Ny)
		uv = .5 * (v + v[:,:,im]) * us # ZE (1~Ny)
		vw = .5 * (v + v[:,km,:]) * ws # XE (1~Ny)
		# get the 3 components of N
		N1[jc] = (uu-uu[:,:,im])[jc] / dx + (uv[jp]-uv[jc]) / dy[jc] + (uw[:,kp]-uw)[jc] / dz
		N2[1:] = (uv[:,:,ip]-uv)[1:] / dx + (vv[1:]-vv[:-1]) / h[1:] + (vw[:,kp]-vw)[1:] / dz
		N3[jc] = (uw[:,:,ip]-uw)[jc] / dx + (vw[jp]-vw[jc]) / dy[jc] + (ww-ww[:,km])[jc] / dz
		N1[0] = N1[-1] = N3[0] = N3[-1] = 0

		# Q = -1/2 * \nabla(N), the returned result's boundary not handled
		return -.5 * Operator(self.para).divergence(N1,N2,N3)

		##### version1: basic (suffer to unacceptable error)
		# ux,vx,wx, uy,vy,wy, uz,vz,wz = Operator(self.para).gradient(u,v,w)
		# return -.5 * (ux**2+vy**2+wz**2 + 2*(uy*vx+uz*wx+wy*vz))

	def rhs_rs(self, u, v, w):
		Nx = self.para.Nx
		Ny = self.para.Ny
		dx = self.para.Lx / Nx
		dy = self.para.dy
		h  = self.para.h

		jc = np.arange(1,Ny)
		jm, jp = jc-1, jc+1
		im, ip = np.arange(-1,Nx-1), np.arange(1,Nx+1) % Nx

		um = np.mean(u, axis=(-1,-2))
		vm = np.mean(v, axis=(-1,-2)) # vm is aligned to y, not yc
		wm = np.mean(w, axis=(-1,-2))
		uf = (u.T - um).T
		vf = (v.T - vm).T
		wf = (w.T - wm).T

		# dUm/dy
		uy = np.empty(Ny+1)
		uy[jc]= .5/dy[jc] * (
			(um[jc]*dy[jp]+um[jp]*dy[jc]) / h[jp] - \
			(um[jc]*dy[jm]+um[jm]*dy[jc]) / h[jc] )
		uy[0] = .5/h[1] * (um[1] -um[0]) # use value on y[1] to approximate yc[0]
		uy[-1]= .5/h[-1]* (um[-1]-um[-2])
		# dv'/dx
		vx = .5/dx * (vf[:,:,ip] - vf[:,:,im])
		vx[jc]= .5 * (vx[jc] + vx[jp])
		vx[0] = h[ 1]/dy[ 1] * (vx[ 1]-vx[ 2]) + .5 * (vx[ 1]+vx[ 2]) # linear extrapolation of v
		vx[-1]= h[-1]/dy[-2] * (vx[-1]-vx[-2]) + .5 * (vx[-1]+vx[-2])

		return -2.*(uy*vx.T).T, 2.*self.get_Q(uf,vf,wf)

	def poisson(self, rhs, bcoef=2):
		''' solve the poisson equation at cell-centers with specified BC
		    the boundaries of rhs are BC and the solution includes boundaries '''
		Nx = self.para.Nx
		Ny = self.para.Ny
		Nz = self.para.Nz
		Nxc= self.para.Nxc
		dx = self.para.Lx / Nx
		dz = self.para.Lz / Nz
		dy = self.para.dy
		h  = self.para.h
		kx = self.para.kx
		kz = self.para.kz

		a = np.empty(Ny+1)
		b = np.empty(Ny+1)
		c = np.empty(Ny+1)
		diag = range(Ny+1)[1:-1]

		# boundary coefficients (a[0] = c[-1] = 0 are not relevant)
		if   bcoef is 2: b[0], c[0], a[-1], b[-1] = -1./h[1], 1./h[1], -1./h[-1],1./h[-1]
		elif bcoef is 1: b[0], c[0], a[-1], b[-1] = 1, 0, 0, 1
		else:            b[0], c[0], a[-1], b[-1] = bcoef
		# differencing coefficients
		a[1:-1] = 1./dy[1:-1] / h[1:-1]
		b[1:-1] =-1./dy[1:-1] * (1./h[1:-1] + 1./h[2:])
		c[1:-1] = 1./dy[1:-1] / h[2:]
		# assemble coefficients into matrix
		mat = np.diag(b)
		mat[1:,:-1] += np.diag(a[1:])
		mat[:-1,1:] += np.diag(c[:-1])
		# solve Poisson equation in spectral space
		phi = spec(rhs)
		for k in range(Nz):
			for i in range(Nxc):
				mat[diag,diag] = b[diag] - ( # kx[i]**2 + kz[k]**2)
					2./dx**2 * (1 - np.cos(kx[i]*dx)) + \
					2./dz**2 * (1 - np.cos(kz[k]*dz)) )
				if bcoef is 2 and not (i or k): mat[1,1] = 1e20
				phi[:,k,i] = np.linalg.solve(mat, phi[:,k,i]) # boundaries are also obtained

		return phys(phi)

	def helmholtz(self, u, v, w):
		''' break a vector into solenoidal & dilatational parts '''
		opt = Operator(self.para)
		dy= self.para.dy.reshape([self.para.Ny+1,1,1])
		h = self.para.h. reshape([self.para.Ny+1,1,1])

		# compute divergence and rotation of the vector field
		phi = opt.divergence(u,v,w)
		xi1, xi2, xi3 = opt.rotation(u,v,w)
		# implement homogeneous BC
		phi[0] = phi[-1] = 0
		xi1[0] = xi1[-1] = 0
		xi2[0] = xi2[-1] = 0
		xi3[0] = xi3[-1] = 0
		# solve poisson equations with Dirichlet BC for xi2 and Riemann BC for others
		phi = self.poisson(phi)
		xi1 = self.poisson(-xi1)
		xi2 = self.poisson(-xi2, bcoef=1)
		xi3 = self.poisson(-xi3)
		# interpolate xi to staggered grids
		xi1 = .5 * (xi1 + np.roll(xi1,1,axis=-1))
		xi3 = .5 * (xi3 + np.roll(xi3,1,axis=-2))
		xi2[1:] = .5/h[1:] * (xi2[:-1]*dy[1:] + xi2[1:]*dy[:-1]) # xi2[0] = 0 has been handled by BC
		# return Us (centered), Ud (staggered), xi (staggered), phi (centered)
		return opt.rotation(xi1,xi2,xi3), opt.gradient(phi), (xi1,xi2,xi3), phi




if __name__ == '__main__':

	# para = DataSetInfo("../test_DAOFW_fake100/base/")
	para = DataSetInfo("../refdata32/")

	tstep = para.tsteps[0]

	u = Field(para).read('U%08i.bin'%tstep, stagtyp=0)
	v = Field(para).read('V%08i.bin'%tstep, stagtyp=0)
	w = Field(para).read('W%08i.bin'%tstep, stagtyp=0)
	p = Field(para).read('P%08i.bin'%tstep, stagtyp=0)

	rhs = 2.*Pressure(para).get_Q(u,v,w)
	rhs[0] = rhs[-1] = 0

	q = Pressure(para).poisson(rhs)


	#################### plot ####################
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt
	from matplotlib.animation import FuncAnimation

	fig, ax = plt.subplots(tight_layout=True, squeeze=True)
	ax.plot(para.yc, np.mean(p,axis=(-1,-2)), marker='o')
	ax.plot(para.yc, np.mean(q,axis=(-1,-2)))
	ax.plot(para.yc, np.mean((p.T-np.mean(p,axis=(-1,-2)))**2, axis=(0,1))**.5, marker='o')
	ax.plot(para.yc, np.mean((q.T-np.mean(q,axis=(-1,-2)))**2, axis=(0,1))**.5)
	ax.legend(['p', 'q'])
	fig.savefig('../PHI.png', dpi=300)
	plt.close()



# if __name__ == '__main__':

# 	para = DataSetInfo("../test_DAOFW_fake100/test/")

# 	if para.Ly < 2.:
# 		para0 = DataSetInfo("../test_DAOFW_fake100/base/")

# 		jb = (para0.Ny - para.Ny) // 2;
# 		para.yc[0] = para0.yc[jb];
# 		para.yc[-1]= para0.yc[para0.Ny-jb];

# 		para.h, para.dy = para.get_Ymesh(para.y, para.yc)

# 	tstep = para.tsteps[500]

# 	u = Field(para).read('U%08i.bin'%tstep, stagtyp=0)
# 	v = Field(para).read('V%08i.bin'%tstep, stagtyp=0)
# 	w = Field(para).read('W%08i.bin'%tstep, stagtyp=0)
# 	p = Field(para).read('P%08i.bin'%tstep, stagtyp=0)

# 	us, ud, xi, q = Pressure(para).helmholtz(u,v,w)


# 	import matplotlib
# 	matplotlib.use('Agg')
# 	import matplotlib.pyplot as plt
# 	from matplotlib.animation import FuncAnimation

# 	fig, ax = plt.subplots(tight_layout=True, squeeze=True)

# 	u = .5 * (u + np.roll(u,-1,axis=-1))
# 	ud[0][:] = .5 * (ud[0] + np.roll(ud[0],-1,axis=-1))

# 	c = ax.contourf(np.arange(para.Nx)*para.Lx/para.Nx, para.yc, (us[0]+ud[0])[:,0,:])
# 	c = ax.contour (np.arange(para.Nx)*para.Lx/para.Nx, para.yc, u            [:,0,:], levels=c.levels, colors='black')

# 	fig.savefig('PHI.png', dpi=300)
# 	plt.close()











# # Gaussian elimination of lower boundary
# a[1] -= b[0] * (a[1]/b[0]) # a[1] = 0
# b[1] -= c[0] * (a[1]/b[0])
# bcl = -rhs[0] * (a[1]/b[0]) # layer j==0
# # Gaussian elimination of upper boundary
# b[-2] -= a[-1] * (c[-2]/b[-1])
# c[-2] -= b[-1] * (c[-2]/b[-1]) # c[-2] = 0
# bcu = -rhs[-1] * (c[-2]/b[-1]) # layer j==Ny


# # chasing method for tri-diag equation
# def chasing(a,b,c,d):
# 	N = len(a)
# 	l = np.zeros(N)
# 	r = np.zeros(N)
# 	m = np.zeros(N)
# 	q = np.zeros(N)

# 	r[0] = b[0]
# 	m[0] = d[0]
# 	for j in range(1,N):
# 		l[j] = a[j] / r[j-1]
# 		r[j] = b[j] - l[j]*c[j-1]
# 		m[j] = d[j] - l[j]*m[j-1]
# 	q[-1] = m[-1] / r[-1]
# 	for j in range(N-1)[::-1]:
# 		q[j] = (m[j] - c[j]*q[j+1]) / r[j]

# 	return q






