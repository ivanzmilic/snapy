import numpy as np
import pyana 
from astropy.io import fits 
import sys
from scipy.interpolate import RegularGridInterpolator
from tqdm import tqdm

# The goal is to rotate, accurately, a 3D atmosphere. x,y do not matter so we do it y,z plane at a time.

# load the atmosphere
file_in = sys.argv[1]

cubein = 0

if (file_in[-3:] == '.f0'):

	cubein = pyana.fzread(file_in)["data"]

elif (file_in[-3:] == 'its'):

	cubein = fits.open(file_in)[0].data

else:

	print("error::the file format is not supported. exiting...")
	exit();

print ("info::atmos shape is: ", cubein.shape)

cubein = cubein[:,:,:,::-1]

NP, NX, NY, NZ = cubein.shape

# The scheme is go x at the time: 

# Specify spacings:
dx = 32.0
dy = 32.0
dz = 16.0

x = np.arange(NX) * dx 
y = np.arange(NY) * dy
z = np.arange(NZ) * dz


# Specify angles:

theta = float(sys.argv[2])
phi = float(sys.argv[3])

# Now, we want to preserve dz, then our steps backwards are:

theta = np.radians(theta)
phi = np.radians(phi)

dzt = dz * np.cos(theta)
dxt = dz * np.sin(theta) * np.cos(phi)
dyt = dz * np.sin(theta) * np.sin(phi)

zt_max = z[-1] / np.cos(theta) # longest possible

NZ_new = int(zt_max // dz) # number of Z in the new cube

# allocate the new cube

cube_tilted = np.zeros([NP,NX,NY,NZ_new])

# And do the log of pressures in the old cube:
cubein[3,:,:,:] = np.log10(cubein[3,:,:,:])
cubein[4,:,:,:] = np.log10(cubein[4,:,:,:])

for p in range(2,NP): # these are the only ones you are interpolating

	interpolator = RegularGridInterpolator((x,y,z), cubein[p])

	print ("info::now interpolating quantity :",p)
	for i in tqdm(range(0,NX)):
		for j in range(0,NY):

			# Create coordinates on the appropriate ray:
			zt = z[-1] - np.arange(NZ_new) * dzt
			xt = x[i] - np.arange(NZ_new) * dxt
			yt = y[j] - np.arange(NZ_new) * dyt

			x_t_min = np.min(xt)
			x_t_max = np.max(xt)
			N_fold_min = (x[0] - x_t_min) // x[-1] + 1 
			N_fold_max = (x_t_max - x[-1]) // x[-1] + 1 
			N_fold = int(np.max([N_fold_min, N_fold_max]))

			for ii in range(0,N_fold):
				x_t_smaller = np.where(xt < x[0])
				x_t_larger = np.where(xt > x[-1])
				xt[x_t_smaller] += x[-1]
				xt[x_t_larger] -= x[-1]

			y_t_min = np.min(yt)
			y_t_max = np.max(yt)
			N_fold_min = (y[0] - y_t_min) // y[-1] + 1 
			N_fold_max = (y_t_max - y[-1]) // y[-1] + 1 
			N_fold = int(np.max([N_fold_min, N_fold_max]))

			for ii in range(0,N_fold):
				y_t_smaller = np.where(yt < y[0])
				y_t_larger = np.where(yt > y[-1])
				yt[y_t_smaller] += y[-1]
				yt[y_t_larger] -= y[-1]

			grid = np.concatenate((xt[:,np.newaxis],yt[:,np.newaxis],zt[:,np.newaxis]),axis=1)
			cube_tilted[p,i,j,:] = interpolator(grid)

cube_tilted[3,:,:,:] = 10.0 ** cube_tilted[3,:,:,:]
cube_tilted[4,:,:,:] = 10.0 ** cube_tilted[4,:,:,:]

pyana.fzwrite("debug.f0", cube_tilted, 0, 'bla')




