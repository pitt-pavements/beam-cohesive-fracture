import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14})

file_list = ['symm_Tt-63Tm+70Tb-63_0.8.inp_coh-stress-form.txt']
xi = [0.8]
theta = [0.9]
numfiles = len(file_list)
h = 150.
ft = -3.4

for i in range(numfiles):
	filename = file_list[i]
	print(filename)
	data = np.loadtxt(filename,unpack=True).T
	data = data[-1,1:]/ft
	if filename == file_list[0]:
		nelem = np.shape(data)[0]
		dh = h/float(nelem)
		z = dh*np.ones(nelem)
		z[0] *= 0.5
		z = np.cumsum(z)/h
	
		plt.plot(data,z,'k',label=r'$\theta = $'+str(theta[i]))

plt.ylim((1.,0.))
plt.xlim((-1.5,3.0))
plt.legend()
plt.xlabel(r'$t_n/f_t$')
plt.ylabel(r'$z/h$')
#plt.savefig('plot.jpg')
plt.show()

