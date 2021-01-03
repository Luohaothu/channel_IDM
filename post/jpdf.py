import numpy as np
from scipy.interpolate import interp2d
import basic


def calc_jpdf(uset_, vset_, n1=100, n2=100):

	uset = np.ravel(uset_)
	vset = np.ravel(vset_)

	if len(uset) != len(vset):
		print('Samples do not align!')
		exit()

	u0, v0 = np.min(uset), np.min(vset)
	u1, v1 = np.max(uset), np.max(vset)

	us = np.linspace(u0, u1, n1)
	vs = np.linspace(v0, v1, n2)

	du = us[1] - us[0]
	dv = vs[1] - vs[0]

	pdf = np.zeros([len(vs), len(us)])

	idx = np.vstack([(uset-u0)/du+.5, (vset-v0)/dv+.5]).astype(int)
	idx, cnt = np.unique(idx, return_counts=True, axis=-1)

	pdf[idx[1], idx[0]] += cnt / len(uset) / (du * dv)

	return us, vs, pdf


def get_samples(para, yp, tskip=1):

	ny = para.Ny
	yps = para.yc / para.lc
	tsteps = para.tsteps[::tskip]
	get_fluc_plus = lambda q: (q - np.mean(q)) / para.uc

	j = np.argmin(np.abs(para.yc - yp*para.lc))

	uset, vset, wset = [], [], []

	for t in tsteps:
		print('reading step %i'%t)

		uset.append( get_fluc_plus(basic.Field(para).read_fluc('U%08i.bin'%t)[[j,ny-j]]) )
		vset.append( get_fluc_plus(basic.Field(para).read_fluc('V%08i.bin'%t)[[j,ny-j]] * np.reshape([1,-1], [-1,1,1])) )
		wset.append( get_fluc_plus(basic.Field(para).read_fluc('W%08i.bin'%t)[[j,ny-j]]) )

	return [np.ravel(_) for _ in (uset, vset, wset)]


def get_jpdf(para, yp, tskip=1):

	yps = para.yc / para.lc

	# find the 2 neighbouring layers
	j0 = np.argmin(np.abs(yps - yp))
	j1 = j0+1 if yps[j0] < yp+1e-6 else 1 if j0==0 else j0-1

	yp0 = yps[j0]
	yp1 = yps[j1]

	# compute JPDF of one layer
	uset, vset, wset = get_samples(para, yp0, tskip=tskip)

	us0, vs0, pdf0_uv = calc_jpdf(uset, vset)
	us0, ws0, pdf0_uw = calc_jpdf(uset, wset)

	# compute JPDF of the other layer
	uset, vset, wset = get_samples(para, yp1, tskip=tskip)

	us1, vs1, pdf1_uv = calc_jpdf(uset, vset)
	us1, ws1, pdf1_uw = calc_jpdf(uset, wset)

	# interpolate layer1 to layer0 and linearly combine them
	us = us0
	vs = vs0
	ws = ws0

	pdf_uv = (yp1-yp)/(yp1-yp0) * pdf0_uv + (yp0-yp)/(yp0-yp1) * interp2d(us1, vs1, pdf1_uv)(us, vs)
	pdf_uw = (yp1-yp)/(yp1-yp0) * pdf0_uw + (yp0-yp)/(yp0-yp1) * interp2d(us1, ws1, pdf1_uw)(us, ws)

	return us, vs, ws, pdf_uv, pdf_uw


def write_jpdf(pame, us, vs, pdf, casename='jpdf'):
	header = \
		'Title = "Joint PDF of u and v"\n' + \
		'variables = "%s", "%s", "%s"\n' \
		% (	'u<sup>+</sup>',
			'v<sup>+</sup>',
			'pdf' ) + \
		'zone t = "%s", i = %i, j = %i' %(casename, len(us), len(vs))

	data = np.empty([3, len(vs), len(us)])

	for j, v in enumerate(vs):
		for i, u in enumerate(us):
			data[:,j,i] = u, v, pdf[j,i]

	data = np.transpose([col.ravel() for col in data])

	np.savetxt(pame, data, header=header, comments='')

	return pame



if __name__ == '__main__':

	para = basic.DataSetInfo('/mnt/TurbNAS/whn/cases_OFW/group2/data8_DSM2000_dx100/test/')
	para.scale_inner()

	us, vs, ws, pdf_uv, pdf_uw = get_jpdf(para, yp=1200, tskip=1)

	write_jpdf('jpdf_uv.dat', us, vs, pdf_uv, casename='DSM2000')
	write_jpdf('jpdf_uw.dat', us, ws, pdf_uw, casename='DSM2000')





