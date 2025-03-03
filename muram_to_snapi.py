import numpy as np
import matplotlib.pyplot as plt 
import pyana
import sys
import muram as mio

# Reads a single muram snapshot and packs it into a .f0 snapi atmosphere: 

path = sys.argv[1]
n_iter = int(sys.argv[2])
type = sys.argv[10]

xmin = int(sys.argv[3])
xmax = int(sys.argv[4])
ymin = int(sys.argv[5])
ymax = int(sys.argv[6])
zmin = int(sys.argv[7])
zmax = int(sys.argv[8])
skip = int(sys.argv[9])

if (type=='muramsub'):

	snap=mio.MuramSubSnap(path,n_iter)

	T = snap.Temp.transpose(1,2,0)
	
	print ("info::muram_binary_loader::the original dimensions are: ", T.shape)

	NX, NY, NZ = T.shape
	
	print ("info :: original dimensions are: ", NX, NY, NZ)

	T = snap.Temp[zmin:zmax, xmin:xmax:skip, ymin:ymax:skip].transpose(1,2,0)

	NX, NY, NZ = T.shape
	print ("info :: output dimensions are: ", NX, NY, NZ)

	atmout = np.zeros([12,NX,NY,NZ])
	Tc = np.copy(T)
	Tc[np.where(T<3200.0)] = 3200.0
	#Tc[np.where(T>50000.0)] = 50000.0
	z = np.arange(zmax-zmin) * 16E5
	tau = np.linspace(-6,2,NZ)
	p = snap.Pres[zmin:zmax, xmin:xmax:skip, ymin:ymax:skip].transpose(1,2,0)
	vz = snap.vx[zmin:zmax, xmin:xmax:skip, ymin:ymax:skip].transpose(1,2,0)

	atmout[0,:,:,:] = tau[None, None, :]
	atmout[1,:,:,:] = z[None, None, :]
	atmout[2,:,:,:] = Tc
	atmout[3,:,:,:] = p
	atmout[4,:,:,:] = p * 0.05
	atmout[9,:,:,:] = vz

	Bz = snap.Bx[zmin:zmax, xmin:xmax:skip, ymin:ymax:skip].transpose(1,2,0) * np.sqrt(4.0*np.pi)
	Bx = snap.By[zmin:zmax, xmin:xmax:skip, ymin:ymax:skip].transpose(1,2,0) * np.sqrt(4.0*np.pi)
	By = snap.Bz[zmin:zmax, xmin:xmax:skip, ymin:ymax:skip].transpose(1,2,0) * np.sqrt(4.0*np.pi)
	
	B = np.sqrt(Bx**2.0 + By**2.0 + Bz**2.0)
	theta = np.arccos(Bz/(B+0.001))
	phi = np.arctan2(By,Bx)

	atmout[7,:,:,:] = B
	atmout[10,:,:,:] = theta
	atmout[11,:,:,:] = phi

	outputname = sys.argv[11] + '_' + sys.argv[2]+ '.f0'

	pyana.fzwrite(outputname, atmout[:,:,:,::-1],0,'bla')

elif (type=='muram'):

	snap=mio.MuramSnap(path,n_iter)

	T = snap.Temp.transpose(1,2,0)
	
	print ("info::muram_binary_loader::the original dimensions are: ", T.shape)

	NX, NY, NZ = T.shape
	
	print ("info :: original dimensions are: ", NX, NY, NZ)

	T = snap.Temp[zmin:zmax, xmin:xmax:skip, ymin:ymax:skip].transpose(1,2,0)

	NX, NY, NZ = T.shape
	print ("info :: output dimensions are: ", NX, NY, NZ)

	atmout = np.zeros([12,NX,NY,NZ])
	Tc = np.copy(T)
	Tc[np.where(T<3200.0)] = 3200.0
	#Tc[np.where(T>50000.0)] = 50000.0
	z = np.arange(zmax-zmin) * 16E5
	tau = np.linspace(-6,2,NZ)
	p = snap.Pres[zmin:zmax, xmin:xmax:skip, ymin:ymax:skip].transpose(1,2,0)
	vz = snap.vx[zmin:zmax, xmin:xmax:skip, ymin:ymax:skip].transpose(1,2,0)

	atmout[0,:,:,:] = tau[None, None, :]
	atmout[1,:,:,:] = z[None, None, :]
	atmout[2,:,:,:] = Tc
	atmout[3,:,:,:] = p
	atmout[4,:,:,:] = p * 0.05
	atmout[9,:,:,:] = vz

	Bz = snap.Bx[zmin:zmax, xmin:xmax:skip, ymin:ymax:skip].transpose(1,2,0) * np.sqrt(4.0*np.pi)
	Bx = snap.By[zmin:zmax, xmin:xmax:skip, ymin:ymax:skip].transpose(1,2,0) * np.sqrt(4.0*np.pi)
	By = snap.Bz[zmin:zmax, xmin:xmax:skip, ymin:ymax:skip].transpose(1,2,0) * np.sqrt(4.0*np.pi)
	
	B = np.sqrt(Bx**2.0 + By**2.0 + Bz**2.0)
	theta = np.arccos(Bz/(B+0.001))
	phi = np.arctan2(By,Bx)

	atmout[7,:,:,:] = B
	atmout[10,:,:,:] = theta
	atmout[11,:,:,:] = phi

	outputname = sys.argv[11] + '_' + sys.argv[2]+ '.f0'

	pyana.fzwrite(outputname, atmout[:,:,:,::-1],0,'bla')	


elif (type=='muramt'):

	# Transposed in a different way to "default one" NX, NZ, NY

	snap=mio.MuramSnap(path,n_iter)

	T = snap.Temp.transpose(0,2,1)
	
	print ("info::muram_binary_loader::the original dimensions are: ", T.shape)

	NX, NY, NZ = T.shape
	
	print ("info :: original dimensions are: ", NX, NY, NZ)

	T = snap.Temp[xmin:xmax:skip,zmin:zmax, ymin:ymax:skip].transpose(0,2,1)

	NX, NY, NZ = T.shape
	print ("info :: output dimensions are: ", NX, NY, NZ)

	atmout = np.zeros([12,NX,NY,NZ])
	Tc = np.copy(T)
	Tc[np.where(T<3200.0)] = 3200.0
	#Tc[np.where(T>50000.0)] = 50000.0
	z = np.arange(zmax-zmin) * 16E5
	tau = np.linspace(-6,2,NZ)
	p = snap.Pres[xmin:xmax:skip, zmin:zmax, ymin:ymax:skip].transpose(0,2,1)
	vz = snap.vy[xmin:xmax:skip, zmin:zmax, ymin:ymax:skip].transpose(0,2,1)

	atmout[0,:,:,:] = tau[None, None, :]
	atmout[1,:,:,:] = z[None, None, :]
	atmout[2,:,:,:] = Tc
	atmout[3,:,:,:] = p
	atmout[4,:,:,:] = p * 0.05
	atmout[9,:,:,:] = vz

	Bz = snap.By[xmin:xmax:skip, zmin:zmax, ymin:ymax:skip].transpose(0,2,1) * np.sqrt(4.0*np.pi)
	Bx = snap.Bx[xmin:xmax:skip, zmin:zmax, ymin:ymax:skip].transpose(0,2,1) * np.sqrt(4.0*np.pi)
	By = snap.Bz[xmin:xmax:skip, zmin:zmax, ymin:ymax:skip].transpose(0,2,1) * np.sqrt(4.0*np.pi)
	
	B = np.sqrt(Bx**2.0 + By**2.0 + Bz**2.0)
	theta = np.arccos(Bz/(B+0.001))
	phi = np.arctan2(By,Bx)

	atmout[7,:,:,:] = B
	atmout[10,:,:,:] = theta
	atmout[11,:,:,:] = phi

	outputname = sys.argv[11] + '_' + sys.argv[2]+ '.f0'

	pyana.fzwrite(outputname, atmout[:,:,:,::-1],0,'bla')	

else:
        print ("Uknown file type. Exiting...")

