# This is the routine made for david where we make an iso-tau cube from a MURAM cube
# Basically, just a copy-paste from the notebook, but with some extra stuff to make it work for any cube

import numpy as np
import matplotlib.pyplot as plt
import sys
import muram as mio

pathsource = sys.argv[1]

path3D = pathsource + '3D/'
path2D = pathsource + '2D/'

snapno = int(sys.argv[2])
cube = mio.MuramSubSnap(pathsource, snapno)

# Little debug line to see the tau:
tau = cube.tau
tau_mean = np.mean(tau.transpose(1,2,0), axis=(0,1))
tau0_index = np.argmin(np.abs(tau_mean - 1), axis=0)  # Find the index where tau is closest to 1
print("info::the tau = 0 index is at: ", tau0_index)


# Copy relevant range of quantities to new arrays, and transpose them to have the shape (x,y,z) instead of (z,x,y)
T_photosphere = np.copy(cube.Temp[:,:,:])
T_photosphere = T_photosphere.transpose(1,2,0)  # Now shape is (x,y,z)
p_photosphere = np.copy(cube.Pres[:,:,:])
p_photosphere = p_photosphere.transpose(1,2,0)  # Now shape is (x,y,z)

# Same for tau 

tau_photosphere = np.copy(cube.tau[:,:,:])
tau_photosphere = tau_photosphere.transpose(1,2,0)  # Now shape is (x,y,z)

Bz_photosphere = np.copy(cube.Bx[:,:,:] * np.sqrt(4.0 * np.pi))
Bz_photosphere = Bz_photosphere.transpose(1,2,0)  # Now shape is (x,y,z)
Bx_photosphere = np.copy(cube.By[:,:,:] * np.sqrt(4.0 * np.pi))
Bx_photosphere = Bx_photosphere.transpose(1,2,0)  # Now shape is (x,y,z)
By_photosphere = np.copy(cube.Bz[:,:,:] * np.sqrt(4.0 * np.pi))
By_photosphere = By_photosphere.transpose(1,2,0)  # Now shape is (x,y,z)

vz_photosphere = np.copy(cube.vx[:,:,:])
vz_photosphere = vz_photosphere.transpose(1,2,0)  # Now shape is (x,y,z)
vx_photosphere = np.copy(cube.vy[:,:,:])
vx_photosphere = vx_photosphere.transpose(1,2,0)  # Now shape is (x,y,z)
vy_photosphere = np.copy(cube.vz[:,:,:])
vy_photosphere = vy_photosphere.transpose(1,2,0)  # Now shape is (x,y,z)

# Ok, finally, let's interpolate our cubes onto a fixed iso-tau grid. 
from tqdm import tqdm
z_og = np.arange(T_photosphere.shape[2]) * 16e5 # height in cm, assuming 20 km grid spacing

tau_values = np.array([1.0,0.1,0.01,0.001])

# Now interpolate:
from scipy.interpolate import interp1d
T_iso_tau = np.zeros((T_photosphere.shape[0], T_photosphere.shape[1], len(tau_values)))
p_iso_tau = np.zeros((T_photosphere.shape[0], p_photosphere.shape[1], len(tau_values)))
vx_iso_tau = np.zeros((T_photosphere.shape[0], vx_photosphere.shape[1], len(tau_values)))
vy_iso_tau = np.zeros((T_photosphere.shape[0], vy_photosphere.shape[1], len(tau_values)))
vz_iso_tau = np.zeros((T_photosphere.shape[0], vz_photosphere.shape[1], len(tau_values)))
Bx_iso_tau = np.zeros((T_photosphere.shape[0], Bx_photosphere.shape[1], len(tau_values)))
By_iso_tau = np.zeros((T_photosphere.shape[0], By_photosphere.shape[1], len(tau_values)))
Bz_iso_tau = np.zeros((T_photosphere.shape[0], Bz_photosphere.shape[1], len(tau_values)))
z_iso_tau = np.zeros((T_photosphere.shape[0], Bz_photosphere.shape[1], len(tau_values)))

from tqdm import tqdm

for i in tqdm(range(T_photosphere.shape[0])):  # over x
    for j in range(T_photosphere.shape[1]):  # over y
        f = interp1d(tau_photosphere[i,j,:], T_photosphere[i,j,:], bounds_error=False, fill_value="extrapolate")
        T_iso_tau[i,j,:] = f(tau_values)
        f = interp1d(tau_photosphere[i,j,:], p_photosphere[i,j,:], bounds_error=False, fill_value="extrapolate")
        p_iso_tau[i,j,:] = f(tau_values)
        f = interp1d(tau_photosphere[i,j,:], vx_photosphere[i,j,:], bounds_error=False, fill_value="extrapolate")
        vx_iso_tau[i,j,:] = f(tau_values)
        f = interp1d(tau_photosphere[i,j,:], vy_photosphere[i,j,:], bounds_error=False, fill_value="extrapolate")
        vy_iso_tau[i,j,:] = f(tau_values)
        f = interp1d(tau_photosphere[i,j,:], vz_photosphere[i,j,:], bounds_error=False, fill_value="extrapolate")
        vz_iso_tau[i,j,:] = f(tau_values)
        f = interp1d(tau_photosphere[i,j,:], Bx_photosphere[i,j,:], bounds_error=False, fill_value="extrapolate")
        Bx_iso_tau[i,j,:] = f(tau_values)
        f = interp1d(tau_photosphere[i,j,:], By_photosphere[i,j,:], bounds_error=False, fill_value="extrapolate")
        By_iso_tau[i,j,:] = f(tau_values)
        f = interp1d(tau_photosphere[i,j,:], Bz_photosphere[i,j,:], bounds_error=False, fill_value="extrapolate")
        Bz_iso_tau[i,j,:] = f(tau_values)
        f = interp1d(tau_photosphere[i,j,:], z_og, bounds_error=False, fill_value="extrapolate")
        z_iso_tau[i,j,:] = f(tau_values)

z_iso_tau -= np.mean(z_iso_tau[:,:,0]) # set z=0 at tau=1 on average

# Looks good, let's pack them in a cube and save to fits:
from astropy.io import fits
hdu = fits.PrimaryHDU(T_iso_tau)
hdu1 = fits.ImageHDU(p_iso_tau)
hdu2 = fits.ImageHDU(vx_iso_tau)
hdu3 = fits.ImageHDU(vy_iso_tau)
hdu4 = fits.ImageHDU(vz_iso_tau)
hdu5 = fits.ImageHDU(Bx_iso_tau)
hdu6 = fits.ImageHDU(By_iso_tau)
hdu7 = fits.ImageHDU(Bz_iso_tau)
hdu8 = fits.ImageHDU(z_iso_tau)
hdul = fits.HDUList([hdu, hdu1, hdu2, hdu3, hdu4, hdu5, hdu6, hdu7, hdu8])
hdul.writeto(pathsource+'muram_iso_tau_cube_'+str(snapno)+'.fits', overwrite=True)



