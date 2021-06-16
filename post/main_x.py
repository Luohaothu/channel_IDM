import numpy as np

from basic import DataSetInfo
from statis import Statis_x
from write_statis import Write


if __name__ == '__main__':

	para = DataSetInfo("/mnt/disk2/whn/etbl/TBL_1420/test/")
	# para.scale_velo(1./22.6906404, 1./492.2115)
	para.scale_none()

	stas = Statis_x(para)
	stas.calc_statis()
	stas.calc_develops(para)

	writer = Write(
		para.postpath,
		casename = para.datapath.split('/')[-2],
		irange = range(0, para.Nx-1),
		jrange = range(0, para.Ny+1) # sum(divmod(para.Ny,2)))
		)

	writer.profiles_x(para, stas)
	writer.develops(para, stas)


