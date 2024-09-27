import numpy as np 
import matplotlib.pyplot as plt 
import pyana
from astropy.io import fits # not necessary but 
import firtez_dz as frz

import sys

def svd_compress(A,limit):
	
	U,s,V = np.linalg.svd(A,full_matrices=True)

	dims = A.shape

	small = np.where(s<max(s)*limit)
	s[small] = 0.0
	s_m = np.zeros(dims)
	#print (small)
	s_m[:dims[0],:dims[0]] = np.diag(s)

	A_sparse = np.dot(U,np.dot(s_m,V))

	return A_sparse

''' This function takes two 3D cubes from snapi and packs them into several individual 3D cubes of 
departure coefficients that can be read by firtez. It is modeled after the make_beta.pro by JMB\
'''

f_lte = sys.argv[1]
f_nlte = sys.argv[2]

atmext = sys.argv[3] 

# In principle one can devise the atmname from other names, but better give the guys a freedom

popslte = pyana.fzread(f_lte)["data"]
popsnlte = pyana.fzread(f_nlte)["data"]

print("info::the dimensions of the population cubes are:", popslte.shape)

dc = popsnlte / popslte

dc = dc[:,:,::-1,:]

tags=['mg5173','ca8542','na5896']
indexlo=[8,16,26]
indexuo=[11,18,27]

indexln=[14,8,26]
indexun=[17,10,27]

nlines = len(tags)



### Now some visualization:
print(dc.shape)
dcmean = np.mean(dc[:10,:10], axis=(0,1))

dcfalc = np.loadtxt(sys.argv[4], skiprows = 1)

print ("info: the dimensions of the referent model atmosphere are:", dcfalc.shape)

for i in range(0, nlines):

	plt.figure(figsize=[7,14])
	plt.subplot(311)
	plt.semilogy(dcfalc[:, indexlo[i]], label='falc')
	plt.semilogy(np.mean(dc[:2,:2], axis=(0,1))[:,indexln[i]], label='2x2')
	plt.semilogy(np.mean(dc[:10,:10], axis=(0,1))[:,indexln[i]], label='10x10')
	plt.semilogy(dcmean[:, indexln[i]], label = 'mean')
	plt.legend()

	plt.subplot(312)
	plt.semilogy(dcfalc[:, indexuo[i]], label='falc')
	plt.semilogy(np.mean(dc[:2,:2], axis=(0,1))[:,indexun[i]], label='2x2')
	plt.semilogy(np.mean(dc[:10,:10], axis=(0,1))[:,indexun[i]], label='10x10')
	plt.semilogy(dcmean[:, indexun[i]], label = 'mean')
	plt.legend()

	plt.subplot(313)
	plt.plot(dcfalc[:, indexuo[i]]/dcfalc[:, indexlo[i]], label='falc')
	plt.plot(dcmean[:, indexun[i]]/dcmean[:, indexln[i]], label='mean')
	plt.legend()

	plt.tight_layout()
	plt.savefig(tags[i]+'_plot.png', bbox_inches='tight')
	plt.close("all")

# --------------------------------------------------------------------------------------------------------------
# define betas and come up with filters to smoothen things nicely:

import numpy as np
import pywt
import sys
import skimage

from skimage.restoration import (denoise_wavelet, estimate_sigma)

NX,NY,NZ = dc[:,:,:,0].shape

betas = np.zeros([nlines, 2, NX, NY, NZ])

for l in range(0, nlines):

	betas[l,0,:,:,:] = dc[:,:,:, indexln[l]]
	betas[l,1,:,:,:] = dc[:,:,:, indexun[l]]


for l in range(0, 0):

	for i in range(0,1):
		for d in range(0,NZ):

			#reconstruct
			#kek = svd_compress(betas[l,i,:,:,d].T, 2E-3)

			kek = betas[l,i,:,:,d]

			coeffs = pywt.dwt2(kek,'db8')

			#assign it to these mighty arrays
			cA,(cH,cV,cD) = coeffs
			#print cA.shape
			#regularize
			cH[:,:] = 0.0
			cV[:,:] = 0.0
			cD[:,:] = 0.0
			#cA[47:57,:] = 0.0

			#reconstruct
			kek = pywt.idwt2(coeffs,'db8')

			if (i == 0 and d == 96):

				print("info::debug line about the dimensions", betas[l,i,:,:,d].shape)
				print("info::debug line about the dimensions", kek.shape)

				plt.figure(figsize=[12,8])
				plt.subplot(221)
				plt.imshow(betas[l,i,:30,:30,d].T, origin='lower', cmap='cividis')
				plt.colorbar()
				plt.subplot(222)
				plt.imshow(kek[:30,:30].T, origin='lower', cmap='cividis')
				plt.colorbar()

				plt.subplot(223)
				plt.imshow(betas[l,i,:,:,d].T, origin='lower', cmap='cividis')
				plt.colorbar()
				plt.subplot(224)
				plt.imshow(kek[:,:].T, origin='lower', cmap='cividis')
				plt.colorbar()
				
				plt.savefig('debug'+tags[l]+'.png', bbox_inches='tight')

			betas[l,i,:,:,d] = kek[:-1,:]

print("info::test line : ", betas.shape)

# OK THIS DOES NOT WORK PROPERLY
#small = np.where(betas < 1E-3)
#betas[small] = 1E-3

#big = np.where(betas > 1E3)
#betas[big] = 1E3

for l in range(0,nlines):

	betafname = 'beta_' + tags[l] + '_'+ atmext + '.bin'

	frz.write_beta_atom(betafname, betas[l,0], betas[l,1])

	testbetafname = 'beta_' + tags[l] + atmext + '_1col.bin'

	frz.write_beta_atom(testbetafname, betas[l,0,0,0].reshape(1,1,-1), betas[l,1,0,0].reshape(1,1,-1))

print ("info::everything seems to be fine, but please check the 1D file")


