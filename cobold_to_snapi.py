import numpy as np
import matplotlib.pyplot as plt 
import pyana
import h5py 
import sys

# Transforms a single atmospheric cube from cobold h5py format into f0 (pyana) snapi format

cobfile = sys.argv[1]

atmin = h5py.File(cobfile,'r')

print ("info:: file keys are: ", atmin.keys())

NX = len(atmin['x'])
NY = len(atmin['y'])
NZ = len(atmin['z'])

print ("info :: original dimensions are: ", NX, NY, NZ)

atmout = np.zeros([12, NX, NY, NZ])

atmout[0,:,:,:] = np.linspace(-6,2,NZ)[None,None,:]
z = atmin['z'][:]
atmout[1,:,:,:] = z[None,None,:]
atmout[2,:,:,:] = atmin['T'][:,:,:].transpose(2,1,0)
atmout[3,:,:,:] = atmin['p'][:,:,:].transpose(2,1,0)
atmout[9,:,:,:] = atmin['vz'][:,:,:].transpose(2,1,0)

B_mag = np.sqrt(atmin['bz'][:,:,:]**2.0 + atmin['bx'][:,:,:]**2.0 + atmin['by'][:,:,:]**2.0)
theta = np.arccos(atmin['bz'][:,:,:]/(B_mag[:,:,:]+0.1)) # make sure it's not dividing by zero
phi   = np.arctan2(atmin['by'][:,:,:], atmin['bx'][:,:,:])

atmout[7,:,:,:] = B_mag.transpose(2,1,0)
atmout[10,:,:,:] = theta.transpose(2,1,0)
atmout[11,:,:,:] = phi.transpose(2,1,0)

# Hardcode cut in z, fix later:

atmout = atmout[:,::2,::2,105:]

print ("info :: final dimensions are: ", atmout.shape)

pyana.fzwrite(sys.argv[2],atmout[:,:,:,::-1],0,'temp')
