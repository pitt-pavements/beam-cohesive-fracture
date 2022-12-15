import numpy as np
import sys

input_file = sys.argv[1]
coh_stress_file = input_file + '_coh-stress.txt'
output_file = input_file + '_coh-stress-form.txt'

t,gp,Tn = np.loadtxt(coh_stress_file,unpack=True)

nsteps = t[len(t)-1]
nsteps = int(nsteps)

nx = len(t)/nsteps
nx = int(nx) #No of gauss points, always even
nelem = nx//2 #number of elements
dn_all = np.zeros((nsteps,nx+1)) #extra 1 for time step
dn_avg = np.zeros((nsteps,nelem+1)) #extra 1 for time step

for i in range(nsteps):
	dn_all[i,0] = i+1
	dn_all[i,1:] = Tn[i*nx:i*nx+(nx)]
	
dn_avg[:,0] = dn_all[:,0]
dn_avg[:,1:] = 0.5*(dn_all[:,1::2] + dn_all[:,2::2]) 
	
np.savetxt(output_file,dn_avg)