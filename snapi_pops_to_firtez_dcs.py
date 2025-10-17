import numpy as np 
import matplotlib.pyplot as plt 
import pyana
from astropy.io import fits # not necessary but 
import firtez_dz as frz

import sys

def svd_compress(A,limit): # WIP
	
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

atmext = sys.argv[3] # how to call the output finally

# In principle one can devise the atmname from other names, but better give the guys a freedom

popslte = pyana.fzread(f_lte)["data"]
popsnlte = pyana.fzread(f_nlte)["data"]

print("info::the dimensions of the population cubes are:", popslte.shape)

dc = popsnlte / popslte

dc = dc[:,:,::-1,:]
del(popslte)
del(popsnlte)

## Below is to be changed according to what you did (which atoms, what transitions, etc...)
## ----------------------------------------------------------------------------------------

# Both sodium lines:
#tags=['na5890','na5896']
#indexln=[14,14]
#indexun=[15,16]

# 4-line setup for SST/CRISP (5172 + 5890 + 6302 + 8542)
tags = ['mg5172', 'na5896','ca8542']
indexln=[6,12,26]
indexun=[9,13,28]


#indexl0=[14,14]
#indexu0=[15,16]

## ----------------------------------------------------------------------------------------

nlines = len(tags)

print("info::total number of lines to model: ", nlines)

### Now some visualization:
print(dc.shape)
NX = dc.shape[0]
NY = dc.shape[1]
dcmean = np.mean(dc, axis=(0,1))


#dcfalc = np.loadtxt(sys.argv[4], skiprows = 1)

#print ("info: the dimensions of the referent model atmosphere are:", dcfalc.shape)

for i in range(0, nlines):

	plt.figure(figsize=[7,14])
	plt.subplot(311)
	#plt.semilogy(dcfalc[:, indexlo[i]], label='falc')
	if (NX > 1 and NY > 1):
		plt.semilogy(np.mean(dc[:2,:2], axis=(0,1))[:,indexln[i]], label='2x2')
	if (NX > 9 and NY > 9):
		plt.semilogy(np.mean(dc[:10,:10], axis=(0,1))[:,indexln[i]], label='10x10')
	plt.semilogy(dcmean[:, indexln[i]], label = 'mean')
	plt.legend()

	plt.subplot(312)
	#plt.semilogy(dcfalc[:, indexuo[i]], label='falc')
	if (NX > 1 and NY > 1):
		plt.semilogy(np.mean(dc[:2,:2], axis=(0,1))[:,indexun[i]], label='2x2')
	if (NX > 9 and NY > 9):
		plt.semilogy(np.mean(dc[:10,:10], axis=(0,1))[:,indexun[i]], label='10x10')
	plt.semilogy(dcmean[:, indexun[i]], label = 'mean')
	plt.legend()

	plt.subplot(313)
	#plt.plot(dcfalc[:, indexuo[i]]/dcfalc[:, indexlo[i]], label='falc')
	plt.plot(dcmean[:, indexun[i]]/dcmean[:, indexln[i]], label='mean')
	plt.legend()

	plt.tight_layout()
	plt.savefig(tags[i]+'_plot.png', bbox_inches='tight')
	plt.close("all")

	# 2D plots of departure coefficients

	d_to_plot = (12,36,60,84)
	NDP = len(d_to_plot)

	plt.figure(figsize=[10.0,17])
	pltindex = 1

	for d in d_to_plot:

		plt.subplot(NDP,2,pltindex)

		plt.imshow(np.log10(dc[:,:,d,indexln[i]].T), origin='lower', cmap='cividis')
		plt.colorbar()

		pltindex+=1

		plt.subplot(NDP,2,pltindex)

		plt.imshow(np.log10(dc[:,:,d,indexun[i]].T), origin='lower', cmap='cividis')
		plt.colorbar()

		pltindex+=1

	plt.savefig(tags[i]+'_maps_plot.png', bbox_inches='tight')	
		
# --------------------------------------------------------------------------------------------------------------
# define betas and come up with filters to smoothen things nicely:

expansion = int(sys.argv[4])

NX,NY,NZ = dc[:,:,:,0].shape

betas = np.zeros([nlines, 2, NX, NY, NZ])

for l in range(0, nlines):

	betas[l,0,:,:,:] = dc[:,:,:, indexln[l]]
	betas[l,1,:,:,:] = dc[:,:,:, indexun[l]]


for l in range(0,nlines):

	betafname = 'beta_' + tags[l] + '_'+ atmext + '.bin'

	frz.write_beta_atom(betafname, betas[l,0], betas[l,1])

	#testbetafname = 'beta_' + tags[l] + atmext + '_1col.bin'

	#frz.write_beta_atom(testbetafname, betas[l,0,0,0].reshape(1,1,-1), betas[l,1,0,0].reshape(1,1,-1))

print ("info::we output the original betas for the specified lines")
print ("info::everything seems to be fine, but please check the 1D file")

# smoothing:
tosmooth = int(sys.argv[5])

if (not tosmooth):
	exit();

import scipy.signal as sgl
# Median filter the data:

betasm = np.copy(betas)

for l in range(0, nlines):
	for i in range(0,2):
		for d in range(0,NZ):
			betasm[l,i,:,:,d] = sgl.medfilt2d(betas[l,i,:,:,d], kernel_size=tosmooth)

for l in range(0,nlines):

	betafname = 'beta_smoothed_' + tags[l] + '_'+ atmext + '.bin'

	frz.write_beta_atom(betafname, betasm[l,0], betasm[l,1])

	#testbetafname = 'beta_smoothed_' + tags[l] + atmext + '_1col.bin'

	#frz.write_beta_atom(testbetafname, betasm[l,0,0,0].reshape(1,1,-1), betasm[l,1,0,0].reshape(1,1,-1))

print ("info::we output the original betas for the specified lines")
print ("info::everything seems to be fine, but please check the 1D file")


del(betas)
del(dc)

print("info::the data was median filtered ", nlines)


for l in range(0, nlines):
	
	d_to_plot = (12,36,60,84)
	NDP = len(d_to_plot)

	plt.figure(figsize=[10.0,17])
	pltindex = 1
	for d in d_to_plot:

		plt.subplot(NDP,2,pltindex)

		plt.imshow(np.log10(betasm[l,0,:,:,d].T), origin='lower', cmap='cividis')
		plt.colorbar()

		pltindex+=1

		plt.subplot(NDP,2,pltindex)

		plt.imshow(np.log10(betasm[l,1,:,:,d].T), origin='lower', cmap='cividis')
		plt.colorbar()

		pltindex+=1


	plt.savefig(tags[l]+'_maps_filtered_plot.png', bbox_inches='tight')	
	plt.close('all')

from scipy.interpolate import RegularGridInterpolator

if (not expansion):
	exit()

# Extrapolate the cube on the desired size
# Also median filter to remove outliers

NX_full = 1358
NY_full = 1321

x_old = np.arange(NX) * 4.0
y_old = np.arange(NY) * 4.0
z_old = np.arange(NZ)

x_new = np.arange(NX_full)
y_new = np.arange(NY_full)
z_new = z_old

X,Y = np.meshgrid(x_new, y_new)

betasl = np.zeros([nlines, 2, NX_full, NY_full, NZ])

for l in range (0,nlines):

	for d in range(0,NZ):

		f = RegularGridInterpolator((x_old,y_old),betasm[l,0,:,:,d],bounds_error=False, fill_value=None)
		betasl[l,0,:,:,d] = f((X,Y)).T
		
		f = RegularGridInterpolator((x_old,y_old),betasm[l,1,:,:,d],bounds_error=False, fill_value=None)
		betasl[l,1,:,:,d] = f((X,Y)).T

		print ('Interpolated for depth = ', d)


		
	d_to_plot = (12,36,60,84)
	NDP = len(d_to_plot)

	plt.figure(figsize=[10.0,17])
	pltindex = 1
	for d in d_to_plot:

		plt.subplot(NDP,2,pltindex)

		plt.imshow(np.log10(betasl[l,0,:,:,d].T), origin='lower', cmap='cividis')
		plt.colorbar()

		pltindex+=1

		plt.subplot(NDP,2,pltindex)

		plt.imshow(np.log10(betasl[l,1,:,:,d].T), origin='lower', cmap='cividis')
		plt.colorbar()

		pltindex+=1


	plt.savefig(tags[l]+'_maps_filtered_expanded_plot.png', bbox_inches='tight')	
	plt.close('all')



for l in range(0,nlines):

	betafname = 'beta_full_' + tags[l] + '_'+ atmext + '.bin'

	frz.write_beta_atom(betafname, betasl[l,0], betasl[l,1])

	#testbetafname = 'beta_full_' + tags[l] + atmext + '_1col.bin'

	#frz.write_beta_atom(testbetafname, betasl[l,0,0,0].reshape(1,1,-1), betasl[l,1,0,0].reshape(1,1,-1))


print ("info::everything seems to be fine, but please check the 1D file")


