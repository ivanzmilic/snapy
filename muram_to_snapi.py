import numpy as np
import matplotlib.pyplot as plt 
import pyana
import h5py 
import sys
import muram as mio

# Reads a single muram snapshot and packs it into a .f0 snapi atmosphere: 

path = sys.argv[1]
it_no = int(sys.argv[2])

snap=muram.MuramSubSnap(path,it_no)

T = snap.Temp.transpose(1,2,0)
	
	print ("info::muram_binary_loader::the original dimensions are: ", T.shape)

NX, NY, NZ = T.shape
print ("info :: dimensions are: ", NX, NY, NZ)

atmout = np.zeros([12, NX, NY, NZ])

#atmout[0,:,:,:] = np.linspace(-6,2,NZ)[None,None,:]
#z = atmin['z'][:]
#atmout[1,:,:,:] = z[None,None,:]
#atmout[2,:,:,:] = atmin['T'][:,:,:].transpose(2,1,0)
#atmout[3,:,:,:] = atmin['p'][:,:,:].transpose(2,1,0)
#atmout[9,:,:,:] = atmin['vz'][:,:,:].transpose(2,1,0)

#B_mag = np.sqrt(atmin['bz'][:,:,:]**2.0 + atmin['bx'][:,:,:]**2.0 + atmin['by'][:,:,:]**2.0)
#theta = np.arccos(atmin['bz'][:,:,:]/(B_mag[:,:,:]+0.1)) # make sure it's not dividing by zero
#phi   = np.arctan(atmin['by'][:,:,:] / atmin['bx'][:,:,:])

#atmout[7,:,:,:] = B_mag.transpose(2,1,0)
#atmout[10,:,:,:] = theta.transpose(2,1,0)
#atmout[11,:,:,:] = phi.transpose(2,1,0)

pyana.fzwrite(sys.argv[2],atmout[:,:,:,::-1],0,'temp')
