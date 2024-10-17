import numpy as np
import matplotlib.pyplot as plt 
import firtez_dz as frz
import sys
import pyana
import os


# Transforms a single atmospheric cube coming from Shah's muram simulations to an adequate SNAPI cube
# magnetic field is basically ignored
# we only have pressure, temperature and los velocity :-) 

ffile = sys.argv[1]

atmin = np.load(ffile)

NT,NZ,NX = atmin[:,:,:,0].shape
atmin = atmin.transpose(0,2,1,3)
print ("info :: dimensions are: ", atmin.shape)


atmout = np.zeros([12, NT, NX, NZ])

atmout[0,:,:,:] = np.linspace(-6,2,NZ)
atmout[1,:,:,:] = np.arange(NZ)*16.0*1E5 # to cm
atmout[2,:,:,:] = atmin[:,:,:,0]
atmout[3,:,:,:] = atmin[:,:,:,2]
atmout[4,:,:,:] = atmout[3,:,:,:] * 0.05
atmout[9,:,:,:] = atmin[:,:,:,1]

# Polish temperatures
print("info::debug, min T = ", np.min(atmout[2]))
print("info::debug, max T =  ", np.max(atmout[2]))

T = np.copy(atmout[2])
smallT = np.where(T<3200)
T[smallT] = 3200.0
atmout[2] = np.copy(T)

print("info::debug, min T = ", np.min(atmout[2]))
print("info::debug, max T =  ", np.max(atmout[2]))


atmout[7,:,:,:] = 0.0
atmout[10,:,:,:] = 0.1
atmout[11,:,:,:] = 0.1

pyana.fzwrite(sys.argv[2],atmout[:,:,:,::-1],0,'temp')
