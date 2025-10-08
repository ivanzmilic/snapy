import numpy as np
import matplotlib.pyplot as plt 
import firtez_dz as frz
import sys
import pyana


# Transforms a single atmospheric cube from firtez-dz format in to snapi format
# mostly used to calculate departure coefficients for JM

ffile = sys.argv[1]

atmin = frz.read_model(ffile)

NX,NY,NZ = atmin.tau.shape

print ("info :: dimensions are: ", NX, NY, NZ)

atmout = np.zeros([12, NX, NY, NZ])

atmout[0,:,:,:] = atmin.tau

if (atmout[0,0,0,0] < 1E-8): # seems like tau is empty:
	print("info::log tau seems to be not provided, initializing my own:")
	atmout[0,:,:,:] = np.linspace(-8,2,NZ)[None,None,:]
atmout[1,:,:,:] = atmin.z * 1E5
atmout[2,:,:,:] = atmin.tem
atmout[3,:,:,:] = atmin.pg
atmout[4,:,:,:] = atmin.pel #atmout[3,:,:,:] * 0.05
if (atmout[4,0,0,0] < 1E-5 * atmout[3,0,0,0]): # electron pressure actually zero:
	print("info::electron pressure seems to be not provided, initializing my own:")
	atmout[4,:,:,:] = atmout[3,:,:,:] * 0.05
atmout[8,:,:,:] = np.sqrt(atmin.vmic)
atmout[9,:,:,:] = atmin.vz * 0.0

# Polish temperatures

#smallT = np.where(atmout[2] < 3200.0)
#atmout[2,smallT] = 3200.0

B_mag = np.sqrt(atmin.bz**2.0 + atmin.bx**2.0 + atmin.by**2.0)
theta = np.arccos(atmin.bz/(B_mag+0.1)) # make sure it's not dividing by zero
phi   = np.arctan2(atmin.by, atmin.bx)

atmout[7,:,:,:] = B_mag
atmout[10,:,:,:] = theta
atmout[11,:,:,:] = phi

# Gonna smooth out some stuff: 
smooth = int(sys.argv[3])
if (smooth == 1):
	print("info::i am going to smooth the atmosphere in all three dimensions")
	from scipy.ndimage import gaussian_filter

	atmout[3] = np.log10(atmout[3])
	for i in range(2,12):
		atmout[i] = gaussian_filter(atmout[i],(1,1,1))

	atmout[9,:,:,:] = 0.0
	atmout[3] = 10.0 ** atmout[3]

if (smooth == 2):
	print("info::i am going to smooth the atmosphere in x and y")
	from scipy.ndimage import gaussian_filter

	atmout[3] = np.log10(atmout[3])
	for i in range(2,12):
		atmout[i] = gaussian_filter(atmout[i],(1,1,0))

	atmout[9,:,:,:] = 0.0
	atmout[3] = 10.0 ** atmout[3]

if (smooth == 3):
	print("info::i am going to smooth the atmosphere in z")
	from scipy.ndimage import gaussian_filter

	atmout[3] = np.log10(atmout[3])
	for i in range(2,12):
		atmout[i] = gaussian_filter(atmout[i],(0,0,1))

	atmout[9,:,:,:] = 0.0
	atmout[3] = 10.0 ** atmout[3]

print("info::smoothing finished. now going to write...")

skip_xy = int(sys.argv[4])

if (NX==1 and NY==1):
	np.savetxt(sys.argv[2], atmout[:,0,0,::-1].T, header=str(NZ)+" WHATWHAT", comments='', fmt="%1.6e")
else:
	pyana.fzwrite(sys.argv[2],atmout[:,::skip_xy,::skip_xy,::-1],0,'temp')
