import numpy as np
import matplotlib.pyplot as plt

t,gp,Tn = np.loadtxt('fort.80',unpack=True)

nsteps = t[len(t)-1]
nsteps = int(nsteps)

nx = len(t)/nsteps
nx = int(nx) #No of gauss points
dn = np.zeros((nsteps,nx+1)) #extra 1 for time step

for i in range(nsteps):
	dn[i,0] = i
	dn[i,1:] = Tn[i*nx:i*nx+(nx)]
	
np.savetxt('Tncoh.txt',dn)