import numpy as np
import matplotlib.pyplot as plt 
from astropy.io import fits
import pyana
import sys

# Transforms a single atmospheric cube from sir fits format into f0 (pyana) snapi format

sirfile = sys.argv[1]

atmin = fits.open(cobfile)[0].data[:,:,::2,::2]

NP, NZ, NX, NY = atmin.shape

print ("info :: dimensions are: ", NX, NY, NZ)

atmin = atmin.transpose(0,2,3,1)

atmout = np.zeros([12, NX, NY, NZ])

atmout[0,:,:,:] = atmin[0,:,:,:]
atmout[1,:,:,:] = atmin[9,:,:,:]
atmout[2,:,:,:] = atmin[1,:,:,:]
atmout[3,:,:,:] = atmin[2,:,:,:]
atmout[9,:,:,:] = atmin[5,:,:,:]*(-1.0)
atmout[7,:,:,:] = atmin[4,:,:,:]
atmout[10,:,:,:] = np.radians(atmin[6,:,:,:])
atmout[11,:,:,:] = np.radians(atmin[7,:,:,:])

pyana.fzwrite(sys.argv[2],atmout[:,:,:,::-1],0,'temp')
