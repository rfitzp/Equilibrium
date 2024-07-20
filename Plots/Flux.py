import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn = 'Equilibrium.nc'
ds = nc.Dataset(fn)
R  = ds['R']
Z  = ds['Z']
Rc = ds['R_c']
Zc = ds['Z_c']

RR = np.asarray(R)
nf = RR.shape[0]

fig = plt.figure(figsize=(8.0, 8.0))
plt.rc('xtick', labelsize=17) 
plt.rc('ytick', labelsize=17) 

plt.subplot(1, 1, 1)

plt.xlim(0.45, 1.55)
plt.ylim(-0.55, 0.55)

for n in range (nf-1):
    plt.plot(R[n],  Z[n],  color = 'blue', linewidth = 0.5, linestyle = 'solid')
plt.plot(R[nf-1],  Z[nf-1],  color = 'red', linewidth = 1.0, linestyle = 'solid')    

plt.plot([1.], [0.], marker='o', markersize=2, color="red")

plt.plot(Rc, Zc, marker='o', markersize=3, color="green", linewidth = 0)

plt.xlabel(r'$R/R_0$', fontsize="20")
plt.ylabel(r'$Z/Z_0$',  fontsize="20")

plt.tight_layout()

plt.show()    
