import numpy as np
import matplotlib.pyplot as plt
import pyana
import h5py
import sys

# Transforms a single atmospheric cube from cobold h5py format into f0 (pyana) snapi format

# Switch: when True, plot x,y horizontal maps of the packed atmosphere for a quick sanity check
PLOT_MAPS = bool(int(sys.argv[3]))

# Undersampling factor in x and y (integer >= 1). skip=1 keeps full resolution.
# Note: undersampling in z is not advised, so it is deliberately not applied there.
SKIP = int(sys.argv[4])
if SKIP < 1:
    sys.exit("error :: skip must be an integer >= 1")

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

# Undersample in x,y by SKIP; hardcode cut in z, fix later:

z_lower = 100

atmout = atmout[:,::SKIP,::SKIP,z_lower:]

print ("info :: final dimensions are: ", atmout.shape)

if PLOT_MAPS:
    # x,y horizontal maps of the packed atmosphere at a chosen height index
    #iz = atmout.shape[3] // 2  # mid-height slice; change as needed
    # iz should be where the mean temperature in x,y plane is closest to 6200 K
    iz = np.argmin(np.abs(atmout[2,:,:,:].mean(axis=(0,1)) - 6200))

    quantities = [
        (2,  'T [K]'),
        (3,  'p [dyn/cm^2]'),
        (9,  'vz [cm/s]'),
        (7,  '|B| [G]'),
        (10, 'theta [rad]'),
        (11, 'phi [rad]'),
    ]

    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    for ax, (idx, label) in zip(axes.ravel(), quantities):
        im = ax.imshow(atmout[idx, :, :, iz].T, origin='lower', aspect='auto')
        ax.set_title(label)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        fig.colorbar(im, ax=ax)
    fig.suptitle('Packed atmosphere x,y maps at iz = %d' % iz)
    fig.tight_layout()

    figname = sys.argv[2] + '_maps.png'
    fig.savefig(figname, dpi=120)
    print ("info :: saved x,y maps to: ", figname)
    plt.show()

pyana.fzwrite(sys.argv[2],atmout[:,:,:,::-1],0,'temp')
