import numpy as np
import matplotlib.pyplot as plt 
from astropy.io import fits 
import pyana 

import sys

# Take a MiHi cube in one temporal instance and prepare stuff for snapi inversion.
# hardcode some stuff because we are in a rush 

stokes = fits.open(sys.argv[1])[0].data

stokes = stokes.astype("double")

ll = fits.open(sys.argv[1])[1].data

qs_mean = np.mean(stokes[:,:,0,220:260])

print("mean qs = ", qs_mean)

stokes /= qs_mean

pyana.fzwrite(sys.argv[2],stokes,0,'bla')

np.savetxt("na_mihi_lgrid.dat", ll.T, header=str(len(ll)), comments='')

mask = np.ones(len(ll))

mask[:20] = 0.0

mask[-20:] = 0.0

np.savetxt("na_mihi_mask.dat", mask.T, header=str(len(ll)), comments='')

