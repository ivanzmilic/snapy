
import numpy as np
import sys

pathsource = sys.argv[1]

path3D = pathsource + '3D/'
path2D = pathsource + '2D/'

#You get the information from the cubes:
header = np.loadtxt(path3D+'Header.024000')

nx = np.int64(header[0])
ny = np.int64(header[1])
nz = np.int64(header[2])
dx = np.float64(header[3])
dy = np.float64(header[4])
dz = np.float64(header[5])
time= np.float64(header[6])
delta_t = np.float64(header[7])
maxva= np.float64(header[8])

print ('nx,ny,nz,dx,dy,dz,time,delta_t,maxva =',nx,ny,nz,dx,dy,dz,time,delta_t,maxva)

#Then you read the actual physical quantities from the cubes like this:
#rho = np.memmap('3D/result_prim_0.024000',dtype=np.float32,mode='r',shape=(nz,ny,nx))

'''result_prim_0 = density (g/cm^3)
result_prim_1 = vx (vertical) (cm/s)
result_prim_2 = vy (horziontal) (cm/s)
result_prim_3 = vz (horizontal) (cm/s)
result_prim_4 = energy (erg/cm^3)
result_prim_5 = Bx (vertical) 1/sqrt(4pi) G
result_prim_6 = By (horziontal) 1/sqrt(4pi) G
result_prim_7 = Bz (horizontal) 1/sqrt(4pi) G
eosT = Temperature (K)
eosP = pressure
tau = continuum optical depth, 10A band @ 500nm
Qres = resistive heating (erg/cm^3/s)
Qvis = viscous heating (erg/cm^3/s)'''

#Then the cubes are all float32 binaries:
# Here comes the eternal question on the velocities direction, and the ordering of the axes:
#rho = np.memmap(path3D+'result_prim_0.'+str(sys.argv[2]),dtype=np.float32,mode='r',shape=(nz,ny,nx))
#vx = np.memmap(path3D+'result_prim_1.'+str(sys.argv[2]),dtype=np.float32,mode='r',shape=(nz,ny,nx))
#vy = np.memmap(path3D+'result_prim_2.'+str(sys.argv[2]),dtype=np.float32,mode='r',shape=(nz,ny,nx))
#vz = np.memmap(path3D+'result_prim_3.'+str(sys.argv[2]),dtype=np.float32,mode='r',shape=(nz,ny,nx))
#energy = np.memmap(path3D+'result_prim_4.'+str(sys.argv[2]),dtype=np.float32,mode='r',shape=(nz,ny,nx))
#Bx = np.memmap(path3D+'result_prim_5.'+str(sys.argv[2]),dtype=np.float32,mode='r',shape=(nz,ny,nx))
#By = np.memmap(path3D+'result_prim_6.'+str(sys.argv[2]),dtype=np.float32,mode='r',shape=(nz,ny,nx))
#Bz = np.memmap(path3D+'result_prim_7.'+str(sys.argv[2]),dtype=np.float32,mode='r',shape=(nz,ny,nx))



#Slices are done in the following way
# The dimensions and physical extent are below, you always have TIME X ONEDIMENSION X OTHERDIMENSION
nx = 1536 
ny = 1536
nz = 1056
Lx = 10e8
Ly = 10e8
Lz =  6.875e8

# To read the information on the time instances themselves, you read the time file:
time = np.fromfile(path2D+'Pore_10Mm_6x6km_res_Bz400G_time_012000_to_044450.dat',dtype=np.float64)
# the file contains the iteration number and simulation solar time [s]
nsnaps = time.size//2
time = time.reshape(2,nsnaps)
print (nsnaps)
print (time[1,:])

# Then we can read the timeseries of specific 2D slices, for instance, the horizontal velocity v_x at the photosphere (tau=1) is in the file:
vx_tau1 = np.memmap(path2D+'Pore_10Mm_6x6km_res_Bz400G_tau_1.000_vy_012000_to_044450.dat',dtype=np.float32,mode='r',shape=(nsnaps,nx,ny))
# Like in the most of the runs, both the arrangement on indices and the vector notation is swapped, so:
# What MURAM calls v_y is actually v_x
# What MURAM calls v_x is actually v_z (vertical)
# What MURAM calls v_z is actually v_y
# Here I assume what you call x is your variable that is contained in the outer loop i.e. slower index, so:
# vx[i,j] -> i corresponds to x and j to y
# For this to make sense, you have to imshow the array with origin='lower', and TRANSPOSED. This is because imshow
# plots images like matrices (first index is which row, second index is which column), and starts from the top.
# So we have:

print(vx_tau1.shape)
# Plotting here just to make sure:
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
plt.figure(figsize= (12.0,10))
plt.imshow(vx_tau1[10,:,:].T/1E5,cmap='bwr',vmin=-5,vmax=5,origin='lower')
plt.tight_layout()
plt.colorbar(label='v_x [km/s]')  
plt.tight_layout()  
plt.savefig('vx_tau1.png',dpi=150, bbox_inches='tight')

plt.close()

# When you look at this v_x should be aligned along the x axis, which given granular velocities makes sense.