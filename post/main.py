import numpy as np
from struct import pack, unpack
from numpy.fft import fft2 as phys
from numpy.fft import ifft2 as spec
from matplotlib import pyplot as plt

# data = np.loadtxt("PROF.dat")
# data[:, 1:] += data[::-1, 1:]
# data[:, 1:] /= 2.0

# np.savetxt("PROF.dat", data)

Nx, Ny, Nz = 192, 193, 36
Lx, Ly, Lz = 3.1416, 2.0, 0.31416
Nxz = Nx * Nz
Nxc = int (Nx/2+1)
Nzc = int (Nz/2+1)


kx = np.hstack( (np.arange(Nxc), np.arange(Nxc-Nx, 0)) ) * (2*np.pi/Lx)
kz = np.hstack( (np.arange(Nzc), np.arange(Nzc-Nz, 0)) ) * (2*np.pi/Lz)

def read_channel(file_path_name):
	recl = (Ny+1) * Nxz * 8;
	with open(file_path_name, 'rb') as fp:
		fp.seek(Nxz*8) # skipping info section
		q = np.reshape( unpack(recl/8*'d', fp.read(recl)), [Ny+1, Nz, Nx] )
	return q - np.reshape( np.mean(q, axis=(-1,-2)), [Ny+1, 1, 1] )

q = read_channel("V00004000.bin")
Evv = abs(spec(q))**2
Evv = np.sum(Evv+Evv[::-1], axis=-1) / 2
plt.plot(np.sort(kz), Evv[15, np.argsort(kz)], '.-')

	
q = read_channel("V00063000.bin")
Evv = abs(spec(q))**2
Evv = np.sum(Evv+Evv[::-1], axis=-1) / 2
plt.plot(np.sort(kz), Evv[15, np.argsort(kz)], '.-')

plt.show()