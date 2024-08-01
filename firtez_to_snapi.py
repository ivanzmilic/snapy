import numpy as np
import matplotlib.pyplot as plt 
import firtez_dz as frz
import sys
import pyana


# Transforms a single atmospheric cube from cobold h5py format into f0 (pyana) snapi format

ffile = sys.argv[1]

atmin = frz.read_model(ffile)

NX,NY,NZ = atmin.tau.shape

print ("info :: dimensions are: ", NX, NY, NZ)

atmout = np.zeros([12, NX, NY, NZ])

atmout[0,:,:,:] = atmin.tau
atmout[1,:,:,:] = atmin.z * 1E5
atmout[2,:,:,:] = atmin.tem
atmout[3,:,:,:] = atmin.pg
atmout[9,:,:,:] = atmin.vz

B_mag = np.sqrt(atmin.bz**2.0 + atmin.bx**2.0 + atmin.by[:,:,:]**2.0)
theta = np.arccos(atmin.vz/(B_mag+0.1)) # make sure it's not dividing by zero
phi   = np.arctan(atmin.by / atmin.bx)

atmout[7,:,:,:] = B_mag
atmout[10,:,:,:] = theta
atmout[11,:,:,:] = phi

pyana.fzwrite(sys.argv[2],atmout[:,:,:,::-1],0,'temp')
