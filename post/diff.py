import numpy as np
from numpy.polynomial import chebyshev as cb


class Schemes:
	def __init__(self):
		pass

	@staticmethod
	def cheb(ys):
		''' derivative matrix for discrete function deployed on full cheb points, f'=D1 f, f''=D2 f '''
		n = len(ys)	
			
		# evalute base functions at discrete points
		T = np.array([np.cos(i*np.arccos(ys)) for i in range(n)])

		M = Schemes.__chebrec_der(n)
		D1 = np.linalg.solve(T, M @ T).T
		D2 = np.linalg.solve(T, M @ M @ T).T

		I1 = np.diag(Schemes.__cheb_quad(n))

		return D1, D2, I1

	@staticmethod
	def compact6(ys):

		CL, CR = Schemes.__compact6_d1(ys)
		D1 = np.linalg.solve(CL, CR)

		CL, CR = Schemes.__compact6_d2(ys)
		D2 = np.linalg.solve(CL, CR)

		# Integration matrix, with homogeneous Dirichlet BC on both ends
		CR[0,0] = CR[-1,-1] = 1e20
		I2 = np.linalg.solve(CR, CL)

		return D1, D2, I2

	@staticmethod
	def compact6_coef(ys, order):
		alf, bet, a, b, c, d, e = (np.zeros(len(ys)) for _ in range(7))

		if   order==1: CL, CR = Schemes.__compact6_d1(ys)
		elif order==2: CL, CR = Schemes.__compact6_d2(ys)

		alf[1:]  = np.diag(CL[1:,:-1])
		bet[:-1] = np.diag(CL[:-1,1:])
		a[2:]    = np.diag(CR[2:,:-2])
		b[1:]    = np.diag(CR[1:,:-1])
		c[:]     = np.diag(CR)
		d[:-1]   = np.diag(CR[:-1,1:])
		e[:-2]   = np.diag(CR[:-2,2:])

		# # Neumann BC coefficients for 2nd order integration
		# c2[0] =  (ys[2]-ys[0])**2                   / ((ys[1]-ys[0]) * (ys[2]-ys[0]) * (ys[2]-ys[1]))
		# d2[0] = -(ys[1]-ys[0])**2                   / ((ys[1]-ys[0]) * (ys[2]-ys[0]) * (ys[2]-ys[1]))
		# e2[0] = ((ys[1]-ys[0])**2-(ys[2]-ys[0])**2) / ((ys[1]-ys[0]) * (ys[2]-ys[0]) * (ys[2]-ys[1]))

		# c2[-1] =  (ys[-3]-ys[-1])**2                     / ((ys[-2]-ys[-1]) * (ys[-3]-ys[-1]) * (ys[-3]-ys[-2]))
		# b2[-1] = -(ys[-2]-ys[-1])**2                     / ((ys[-2]-ys[-1]) * (ys[-3]-ys[-1]) * (ys[-3]-ys[-2]))
		# a2[-1] = ((ys[-2]-ys[-1])**2-(ys[-3]-ys[-1])**2) / ((ys[-2]-ys[-1]) * (ys[-3]-ys[-1]) * (ys[-3]-ys[-2]))

		return alf, bet, a, b, c, d, e

	@staticmethod
	def center2(ys):
		n = len(ys)

		D1 = np.zeros([n,n])
		D2 = np.zeros([n,n])

		for j, y in enumerate(ys):

			if       j ==  0: j1, j2, j3 = 0,   1, 2
			elif 0 < j < n-1: j1, j2, j3 = j-1, j, j+1
			else:             j1, j2, j3 = -3, -2, -1

			y1, y2, y3 = ys[j1], ys[j2], ys[j3]

			mat = np.linalg.inv([
				[1, y1-y, .5*(y1-y)**2],
				[1, y2-y, .5*(y2-y)**2],
				[1, y3-y, .5*(y3-y)**2]])

			D1[j,[j1,j2,j3]] = mat[1]
			D2[j,[j1,j2,j3]] = mat[2]

		I2 = D2.copy()
		I2[0,0] = I2[-1,-1] = 1e20
		I2 = np.linalg.inv(I2)

		return D1, D2, I2


	@staticmethod
	def __chebrec_der(n):
		''' recurrence relationship of the derivative of Chebyshev polynomials
			2 T_n = 1/(n+1) T_{n+1}^' - 1/(n-1) T_{n-1}^' '''
		M = np.zeros([n,n])
		for j in range(n):
			k = j//2
			if j%2:
				for i in range(k+1):
					M[j,2*i] = 4 * k + 2
				M[j,0] /= 2
			else:
				for i in range(1,k+1):
					M[j,2*i-1] = 4 * k
		return M

	@staticmethod
	def __chebrec_int(n):
		''' recurrence relationship of the integration of Chebyshev polynomials
			int_0^x T_n = 1/2 ( 1/(n+1) T_{n+1} - 1/(n-1) T_{n-1} ) '''
		M = np.zeros([n,n])
		for j in range(n):
			if   j == 0:   M[j, 1       ] = 1
			elif j == 1:   M[j,(0,2)    ] = .25
			elif j == n-1: M[j, j-1     ] = -.5/(j-1) # higher-order are considered to have zero weight
			else:          M[j,(j-1,j+1)] = -.5/(j-1), .5/(j+1)
		return M

	@staticmethod
	def __cheb_quad(n):
		''' weights for Chebyshev–Gauss quadrature,
			https://keisan.casio.com/exec/system/1282551763 '''
		W = np.zeros(n)
		a = np.pi / (n-1)
		for i,y in enumerate(cb.chebpts2(n)):
			if 0 < i < n-1:
				W[i] = a * np.sin(a*i)**2 / (1 - y**2)**.5
		return W / 2 # /2 for normalizing the integration on range [-1,1]

	@staticmethod
	def __compact6_d1(ys):
		''' CL f' = CR f, with CL tri-diagonal and CR penta-diagonal,
			alf f'_2 + f'_3 + bet f'_4 = a f_1 + b f_2 + c f_3 + d f_4 + e f_5 '''
		n = len(ys)
		CL = np.eye(n)
		CR = np.zeros([n, n])

		for j, y in enumerate(ys):

			j1, j2, j3, j4, j5 = j-2, j-1, j, j+1, j+2

			hm2 = ys[j1] - y
			hm1 = ys[j2] - y
			hp1 = ys[j4%n] - y
			hp2 = ys[j5%n] - y

			if 1 < j < n-2:
				# Taylor series coefficients
				mat = np.array([
					[1,      1,      1, 1,      1,       0,         0       ],
					[hm2,    hm1,    0, hp1,    hp2,    -1,        -1       ],
					[hm2**2, hm1**2, 0, hp1**2, hp2**2, -2*hm1,    -2*hp1   ],
					[hm2**3, hm1**3, 0, hp1**3, hp2**3, -3*hm1**2, -3*hp1**2],
					[hm2**4, hm1**4, 0, hp1**4, hp2**4, -4*hm1**3, -4*hp1**3],
					[hm2**5, hm1**5, 0, hp1**5, hp2**5, -5*hm1**4, -5*hp1**4],
					[hm2**6, hm1**6, 0, hp1**6, hp2**6, -6*hm1**5, -6*hp1**5]])
				# Taylor series combination
				rhs = np.array([0, 1, 0, 0, 0, 0, 0])
				# stencil point weights
				a, b, c, d, e, alf, bet = np.linalg.solve(mat, rhs)
				# put weights into the whole matrix
				CL[j, [j2, j4]] = alf, bet
				CR[j, [j1, j2, j3, j4, j5]] = a, b, c, d, e

			elif j == 1 or j == n-2:

				mat = np.array([
					[1,      1, 1,       0,         0       ],
					[hm1,    0, hp1,    -1,        -1       ],
					[hm1**2, 0, hp1**2, -2*hm1,    -2*hp1   ],
					[hm1**3, 0, hp1**3, -3*hm1**2, -3*hp1**2],
					[hm1**4, 0, hp1**4, -4*hm1**3, -4*hp1**3]])

				rhs = np.array([0, 1, 0, 0, 0])

				b, c, d, alf, bet = np.linalg.solve(mat, rhs)

				CL[j, [j2, j4]] = alf, bet
				CR[j, [j2, j3, j4]] = b, c, d

			elif j == 0:

				mat = np.array([
					[1, 1,      1,       0       ],
					[0, hp1,    hp2,    -1       ],
					[0, hp1**2, hp2**2, -2*hp1   ],
					[0, hp1**3, hp2**3, -3*hp1**2]])

				rhs = np.array([0, 1, 0, 0])

				c, d, e, bet = np.linalg.solve(mat, rhs)

				CL[j, j4] = bet
				CR[j, [j3, j4, j5]] = c, d, e

			elif j == n-1:

				mat = np.array([
					[1, 1,      1,       0       ],
					[0, hm1,    hm2,    -1       ],
					[0, hm1**2, hm2**2, -2*hm1   ],
					[0, hm1**3, hm2**3, -3*hm1**2]])

				rhs = np.array([0, 1, 0, 0])

				c, b, a, alf = np.linalg.solve(mat, rhs)

				CL[j, j2] = alf
				CR[j, [j1, j2, j3]] = a, b, c

		return CL, CR

	@staticmethod
	def __compact6_d2(ys):
		''' CL f'' = CR f, with CL tri-diagonal and CR penta-diagonal,
			alf f''_2 + f''_3 + bet f''_4 = a f_1 + b f_2 + c f_3 + d f_4 + e f_5 '''
		n = len(ys)
		CL = np.eye(n)
		CR = np.zeros([n, n])

		for j, y in enumerate(ys):

			j1, j2, j3, j4, j5 = j-2, j-1, j, j+1, j+2

			hm2 = ys[j1] - y
			hm1 = ys[j2] - y
			hp1 = ys[j4%n] - y
			hp2 = ys[j5%n] - y

			if 1 < j < n-2:
				# Taylor series coefficients
				mat = np.array([
					[1,      1,      1, 1,      1,       0,          0        ],
					[hm2,    hm1,    0, hp1,    hp2,     0,          0        ],
					[hm2**2, hm1**2, 0, hp1**2, hp2**2, -2,         -2        ],
					[hm2**3, hm1**3, 0, hp1**3, hp2**3, -6 *hm1,    -6 *hp1   ],
					[hm2**4, hm1**4, 0, hp1**4, hp2**4, -12*hm1**2, -12*hp1**2],
					[hm2**5, hm1**5, 0, hp1**5, hp2**5, -20*hm1**3, -20*hp1**3],
					[hm2**6, hm1**6, 0, hp1**6, hp2**6, -30*hm1**4, -30*hp1**4]])
				# Taylor series combination
				rhs = np.array([0, 0, 2, 0, 0, 0, 0])
				# stencil point weights
				a, b, c, d, e, alf, bet = np.linalg.solve(mat, rhs)
				# put weights into the whole matrix
				CL[j, [j2, j4]] = alf, bet
				CR[j, [j1, j2, j3, j4, j5]] = a, b, c, d, e

			elif j == 1:

				mat = np.array([
					[1,      1, 1,      1,       0        ],
					[hm1,    0, hp1,    hp2,     0        ],
					[hm1**2, 0, hp1**2, hp2**2, -2        ],
					[hm1**3, 0, hp1**3, hp2**3, -6 *hp1   ],
					[hm1**4, 0, hp1**4, hp2**4, -12*hp1**2]])

				rhs = np.array([0, 0, 2, 0, 0])

				b, c, d, e, bet = np.linalg.solve(mat, rhs)

				CL[j, j4] = bet
				CR[j, [j2, j3, j4, j5]] = b, c, d, e

			elif j == n-2:

				mat = np.array([
					[1,      1, 1,      1,       0        ],
					[hp1,    0, hm1,    hm2,     0        ],
					[hp1**2, 0, hm1**2, hm2**2, -2        ],
					[hp1**3, 0, hm1**3, hm2**3, -6 *hm1   ],
					[hp1**4, 0, hm1**4, hm2**4, -12*hm1**2]])

				rhs = np.array([0, 0, 2, 0, 0])

				d, c, b, a, alf = np.linalg.solve(mat, rhs)

				CL[j, j2] = alf
				CR[j, [j1, j2, j3, j4]] = a, b, c, d

			# # singularity occurs when grid is nearly uniform
			# elif j == 0:
			# 	mat = np.array([
			# 		[1, 1,      1,       0    ],
			# 		[0, hp1,    hp2,     0    ],
			# 		[0, hp1**2, hp2**2, -2    ],
			# 		[0, hp1**3, hp2**3, -6*hp1]])
			# 	rhs = np.array([0, 0, 2, 0])
			# 	c, d, e, bet = np.linalg.solve(mat, rhs)
			# 	CL[j, j4] = bet
			# 	CR[j, [j3, j4, j5]] = c, d, e
			# elif j == n-1:
			# 	mat = np.array([
			# 		[1, 1,      1,       0    ],
			# 		[0, hm1,    hm2,     0    ],
			# 		[0, hm1**2, hm2**2, -2    ],
			# 		[0, hm1**3, hm2**3, -6*hm1]])
			# 	rhs = np.array([0, 0, 2, 0])
			# 	c, b, a, alf = np.linalg.solve(mat, rhs)
			# 	CL[j, j2] = alf
			# 	CR[j, [j1, j2, j3]] = a, b, c

			elif j == 0:   CR[j, [j3, j4, j5]] = 2/hp1/hp2 * np.array([1, -hp2/(hp2-hp1), hp1/(hp2-hp1)])
			elif j == n-1: CR[j, [j3, j2, j1]] = 2/hm1/hm2 * np.array([1, -hm2/(hm2-hm1), hm1/(hm2-hm1)])

		return CL, CR


	@staticmethod
	def check(ys, us, D1, D2, I2=None):
		from matplotlib import pyplot as plt

		du1 = np.gradient(us, ys, edge_order=2)
		du2 = np.gradient(np.gradient(us, ys, edge_order=2), ys, edge_order=2)

		scl1 = np.std(us) / np.std(du1) if np.std(du1) > 1e-10 else 1
		scl2 = np.std(us) / np.std(du2) if np.std(du2) > 1e-10 else 1

		plt.plot(ys, us, 'k')
		plt.plot(ys, du1*scl1, 'k')
		plt.plot(ys, du2*scl2, 'k')

		du1 = D1 @ us
		du2 = D2 @ us
		if I2 is not None: idu = I2 @ du2

		if I2 is not None: plt.plot(ys, idu, '.', label='0th')
		plt.plot(ys, du1*scl1, '.', label='1th')
		plt.plot(ys, du2*scl2, '.', label='2th')

		plt.legend()
		plt.show()

		plt.close()



if __name__ == '__main__':

	Ny = 201
	Ly = 1.

	# ys = Ly * np.linspace(0, 1, Ny)
	# ys = Ly * (1 - np.cos(.5*np.pi*np.linspace(0,1,Ny)))
	ys = cb.chebpts2(Ny)

	# D1, D2, I2 = Schemes.compact6(ys)
	# D1, D2, I2 = Schemes.center2(ys)
	D1, D2, I1 = Schemes.cheb(ys)

	# us = np.sin(2*np.pi*ys) * np.cos(np.pi*ys**2)
	us = (ys - ys[0]) * (ys[-1] - ys)

	print(cb.chebval(ys, cb.chebint(cb.chebfit(ys, us, len(ys)-1), lbnd=ys[0]))[-1])
	print(np.diag(I1) @ us * (ys[-1]-ys[0]))

	Schemes.check(ys, us, D1, D2)







