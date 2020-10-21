#!/root/Software/anaconda3/bin/python3
import numpy as np
import torch as tc
import numpy.fft as ft
from scipy.signal import hilbert

'''
Pack up some often used no-state functions.
Use class for inheritance when some functions need modification
'''

class Tools:

	@staticmethod
	def spec(q): return ft.ifft(ft.ihfft(q), axis=-2)
	@staticmethod
	def phys(q): return ft.hfft(ft.fft(q, axis=-2)) # Nx must be even

	@staticmethod
	def envelop(q): return np.abs(hilbert(q))
	@staticmethod
	def envelup(q): return np.abs(hilbert(q * (q>0))) * 2
	@staticmethod
	def envelow(q): return np.abs(hilbert(q * (q<0))) * 2

	@staticmethod
	def normalize(q): return q / np.max(np.abs(q))
	
	@staticmethod
	def roll2(q, i, k): return np.roll(np.roll(q, i, axis=-1), k, axis=-2)

	@staticmethod
	def corr(u, v):
		''' <u'v'> / sqrt(<u'^2><v'^2>) '''
		return (np.mean(u*v) - (np.mean(u)*np.mean(v))) / (np.std(u) * np.std(v))

	@staticmethod
	def corr2p(u, v):
		''' <u'(x,z)v'(x+dx,z+dz)> / sqrt(<u'^2><v'^2>) '''
		fu = spec(u)
		fv = spec(v)
		fu.T[0,0] -= np.mean(fu.T[0,0])
		fv.T[0,0] -= np.mean(fv.T[0,0])
		return phys(np.conj(fu) * fv) / (np.std(u) * np.std(v))

	@staticmethod
	def newton(pfun, gma0, maxiter=20):
		''' iteratively optimize the residual function pfun(gma) towards 0 '''
		gma, dgma = gma0, .1 * gma0

		F0 = pfun(gma-dgma)

		for i in range(maxiter):

			F = pfun(gma)

			grad = (F-F0) / dgma
			if not np.all(grad): break
			dgma = - F / grad
			if not np.all(dgma): break
			gma += dgma

			F0 = F

		else:
			return gma
		
		print('Iteration terminated at ', i)
		return gma

	@staticmethod
	def chasing3(a, b, c, d):
		''' alternate direction method for tri-diag equation,
			dimensions other than 0 are batch dimensions for d '''

		# a |b c      |      l |1        |     |u c      |
		#   |a b c    |        |l 1      |     |  u c    |
		#   |  a b c  |   ->   |  l 1    |  *  |    u c  |
		#   |    a b c|        |    l 1  |     |      u c|
		#   |      a b| c      |      l 1|     |        u| c

		N = len(a)
		if not (len(b)==len(c)==len(d)==N): print('dimension incompatible')

		r = np.empty_like(a)
		y = np.empty_like(d)

		r[0] = b[0]
		y[0] = d[0]
		for j in range(1,N):
			l    = a[j] / r[j-1]
			r[j] = b[j] - l*c[j-1]
			y[j] = d[j] - l*y[j-1]
			
		y[-1] /= r[-1]
		for j in range(N-1)[::-1]:
			y[j] = (y[j] - c[j]*y[j+1]) / r[j]

		return y

	@staticmethod
	def chasing5(a, b, c, d, e, f):
		''' alternate direction method for penta-diag equation,
			dimensions other than 0 are batch dimensions for f '''

		# a b |c d e    |      zt gm |af            |     |1  bt qt      |
		#   a |b c d e  |         zt |gm af         |     |   1  bt qt   |
		#     |a b c d e|    ->      |zt gm af      |  *  |      1  bt qt|
		#     |  a b c d| e          |   zt gm af   |     |         1  bt| qt
		#     |    a b c| d e        |      zt gm af|     |            1 | bt qt
		
		N  = len(a)
		if not (len(b)==len(c)==len(d)==len(e)==len(f)==N): print('dimension incompatible')

		bt = np.empty_like(a)
		qt = np.empty_like(a)
		y  = np.empty_like(f)

		af0   = c[0]
		bt[0] = d[0] / af0
		gm    = b[1]
		af1   = c[1] - gm * bt[0]

		y[0] =  f[0] / af0
		y[1] = (f[1] - gm * y[0]) / af1
		
		# direction 1
		for j in range(2,N):
			qt[j-2] = e[j-2] / af0
			bt[j-1] = (d[j-1] - gm * qt[j-2]) / af1
			zt = a[j]
			gm = b[j] - zt * bt[j-2]
			af = c[j] - zt * qt[j-2] - gm * bt[j-1]

			y[j] = (f[j] - zt * y[j-2] - gm * y[j-1]) / af

			af0, af1 = af1, af

		# direction 2
		y[-2] -= bt[-2] * y[-1]
		for j in range(N-2)[::-1]:
			y[j] = y[j] - qt[j] * y[j+2] - bt[j] * y[j+1]

		return y


class Tools_cuda(Tools):

	@staticmethod
	def corr(u, v):
		''' <u'v'> / sqrt(<u'^2><v'^2>) '''
		cor = tc.mean(u*v) - tc.mean(u) * tc.mean(v)
		norm = tc.std(u, unbiased=False) * tc.std(v, unbiased=False)
		return (cor/norm).item()

	@staticmethod
	def corr2p(u, v):
		''' <u'(x,z)v'(x+dx,z+dz)> / sqrt(<u'^2><v'^2>) '''
		fur, fui = tc.rfft(u, signal_ndim=2).T
		fvr, fvi = tc.rfft(v, signal_ndim=2).T
		fur[0,0] -= tc.mean(fur[0,0])
		fui[0,0] -= tc.mean(fui[0,0])
		fvr[0,0] -= tc.mean(fvr[0,0])
		fvi[0,0] -= tc.mean(fvi[0,0])
		n2, n1 = u.shape[-2:]
		e = tc.stack([fur*fvr+fui*fvi, fur*fvi-fui*fvr]).T / (n2*n1)
		norm = tc.std(u, unbiased=False) * tc.std(v, unbiased=False)
		return tc.irfft(e, signal_ndim=2, signal_sizes=(n2,n1)) / norm

	@staticmethod
	def envelop(q):
		''' mannually implement scipy.signal.hilbert() by fft on tc tensors '''
		fq = tc.rfft(q, signal_ndim=1)
		fq.T[:,1:-1] *= 2
		qa = tc.cat([fq, tc.zeros(fq.T[:,1:-1].T.shape).cuda()], dim=-2)
		qa = tc.ifft(qa, signal_ndim=1)
		return (qa.T[0]**2 + qa.T[1]**2).T**.5
	@staticmethod
	def envelup(q):
		return Tools_cuda.envelop(q * (q>0)) * 2		
	@staticmethod
	def envelow(q):
		return Tools_cuda.envelop(q * (q<0)) * 2

	@staticmethod
	def spec(q):
		''' do the same spec on GPU '''
		p = tc.tensor(q.astype(np.float32), device='cuda')
		pr, pi = tc.rfft(p, signal_ndim=2).T
		pr /= tc.prod(tc.tensor(p.shape[-2:]))
		pi /= tc.prod(tc.tensor(p.shape[-2:])) * (-1)
		return pr.T.cpu().numpy() + 1j * pi.T.cpu().numpy()

	@staticmethod
	def phys(q):
		''' do the same phys on GPU '''
		nz, nxc = q.shape[-2:]
		pr = tc.tensor(q.real.astype(np.float32), device='cuda')
		pi = tc.tensor(q.imag.astype(np.float32), device='cuda') * (-1)
		p  = tc.irfft(tc.stack([pr.T, pi.T]).T, signal_ndim=2, signal_sizes=[nz, (nxc-1)*2]) # Nx must be even
		p *= tc.prod(tc.tensor(p.shape[-2:]))
		return p.cpu().numpy()

	@staticmethod
	def chasing3(a_, b_, c_, d_):
		''' alternate direction method for tri-diag equation on GPU,
			dimensions other than 0 are batch dimensions for d '''
		a = tc.tensor(a_.astype(np.float32), device='cuda')
		b = tc.tensor(b_.astype(np.float32), device='cuda')
		c = tc.tensor(c_.astype(np.float32), device='cuda')
		d = tc.stack([ # treat real & imag parts of d as a trailing batch dimension
			tc.tensor(d_.real.astype(np.float32), device='cuda').T,
			tc.tensor(d_.imag.astype(np.float32), device='cuda').T]).T

		N = len(a)
		if not (len(b)==len(c)==len(d)==N): print('dimension incompatible')

		r = tc.empty_like(a)
		y = tc.empty_like(d)

		r[0] = b[0]
		y[0] = d[0]
		for j in range(1,N):
			l    = a[j] / r[j-1]
			r[j] = b[j] - l*c[j-1]
			y[j] = d[j] - l*y[j-1]
			
		y[-1] /= r[-1]
		for j in range(N-1)[::-1]:
			y[j] = (y[j] - c[j]*y[j+1]) / r[j]

		return y.T[0].T.cpu().numpy() + 1j * y.T[1].T.cpu().numpy()

	@staticmethod
	def chasing5(a_, b_, c_, d_, e_, f_):
		''' alternate direction method for penta-diag equation on GPU,
			dimensions other than 0 are batch dimensions for f,
			a~e can be 1D or same size as f '''
		a = tc.tensor(a_.astype(np.float32), device='cuda')
		b = tc.tensor(b_.astype(np.float32), device='cuda')
		c = tc.tensor(c_.astype(np.float32), device='cuda')
		d = tc.tensor(d_.astype(np.float32), device='cuda')
		e = tc.tensor(e_.astype(np.float32), device='cuda')
		f = tc.stack([ # treat real & imag parts of d as a trailing batch dimension
			tc.tensor(f_.real.astype(np.float32), device='cuda').T,
			tc.tensor(f_.imag.astype(np.float32), device='cuda').T]).T
		
		N  = len(a)
		if not (len(b)==len(c)==len(d)==len(e)==len(f)==N): print('dimension incompatible')

		bt = tc.empty_like(a)
		qt = tc.empty_like(a)
		y  = tc.empty_like(f)

		af0   = c[0]
		bt[0] = d[0] / af0
		gm    = b[1]
		af1   = c[1] - gm * bt[0]

		y[0].T[:] =  f[0].T / af0.T
		y[1].T[:] = (f[1].T - gm.T * y[0].T) / af1.T
		
		# direction 1
		for j in range(2,N):
			qt[j-2] = e[j-2] / af0
			bt[j-1] = (d[j-1] - gm * qt[j-2]) / af1
			zt = a[j]
			gm = b[j] - zt * bt[j-2]
			af = c[j] - zt * qt[j-2] - gm * bt[j-1]

			y[j].T[:] = (f[j].T - zt.T * y[j-2].T - gm.T* y[j-1].T) / af.T

			af0, af1 = af1, af

		# direction 2
		y[-2].T[:] -= bt[-2].T * y[-1].T
		for j in range(N-2)[::-1]:
			y[j].T[:] = y[j].T - qt[j].T * y[j+2].T - bt[j].T * y[j+1].T

		return y.T[0].T.cpu().numpy() + 1j * y.T[1].T.cpu().numpy()



if __name__ == '__main__':

	# test codes for chasing
	N = 30
	M = 5

	a = np.random.rand(N)
	b = np.random.rand(N)
	c = np.random.rand(N)
	d = np.random.rand(N)
	e = np.random.rand(N)
	f = np.random.rand(N,M) + 1j * np.random.rand(N,M)

	mat = np.diag(b)
	mat[1:,:-1] += np.diag(a[1:])
	mat[:-1,1:] += np.diag(c[:-1])
	print(np.mean(np.abs( Tools.chasing3(a,b,c,f) - np.linalg.solve(mat,f) )))

	mat = np.diag(c)
	mat[2:,:-2] += np.diag(a[2:])
	mat[1:,:-1] += np.diag(b[1:])
	mat[:-1,1:] += np.diag(d[:-1])
	mat[:-2,2:] += np.diag(e[:-2])
	print(np.mean(np.abs( Tools.chasing5(a,b,c,d,e,f) - np.linalg.solve(mat,f) )))


	