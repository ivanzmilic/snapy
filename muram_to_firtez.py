import numpy as np
import matplotlib.pyplot as plt 
import pyana
import sys
import muram as mio
import firtez_dz as frz

# Reads a single muram snapshot and packs it into a .bin firtez atmosphere 

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

	atmout = frz.atm_model3D(NX, NY, NZ)
	
	Tc = np.copy(T)
	Tc[np.where(T<3200.0)] = 3200.0
	#Tc[np.where(T>50000.0)] = 50000.0
	z = np.arange(zmax-zmin) * 16.0
	z3d = np.zeros([NX, NY, NZ])
	z3d[:,:,:] = z[None,None,:]
	tau = np.linspace(-6,2,NZ) # Should not matter
	p = snap.Pres[zmin:zmax, xmin:xmax:skip, ymin:ymax:skip].transpose(1,2,0)
	vz = snap.vx[zmin:zmax, xmin:xmax:skip, ymin:ymax:skip].transpose(1,2,0)

	atmout.set_z(z3d)
	atmout.set_pg(p)
	atmout.set_tem(Tc)
	atmout.set_vz(-vz)

	Bz = snap.Bx[zmin:zmax, xmin:xmax:skip, ymin:ymax:skip].transpose(1,2,0) * np.sqrt(4.0*np.pi)
	Bx = snap.By[zmin:zmax, xmin:xmax:skip, ymin:ymax:skip].transpose(1,2,0) * np.sqrt(4.0*np.pi)
	By = snap.Bz[zmin:zmax, xmin:xmax:skip, ymin:ymax:skip].transpose(1,2,0) * np.sqrt(4.0*np.pi)
	
	atmout.set_bx(Bx)
	atmout.set_by(By)
	atmout.set_bz(Bz)
	
	outputname = sys.argv[11] + '_' + sys.argv[2]+ '.bin'

	atmout.write_model(outputname)

elif (type=='muram'):

	snap=mio.MuramSnap(path,n_iter)

	T = snap.Temp.transpose(1,2,0)
	
	print ("info::muram_binary_loader::the original dimensions are: ", T.shape)

	NX, NY, NZ = T.shape
	
	print ("info :: original dimensions are: ", NX, NY, NZ)

	T = snap.Temp[zmin:zmax, xmin:xmax:skip, ymin:ymax:skip].transpose(1,2,0)

	NX, NY, NZ = T.shape
	print ("info :: output dimensions are: ", NX, NY, NZ)

	atmout = frz.atm_model3D(NX, NY, NZ)
	
	Tc = np.copy(T)
	Tc[np.where(T<3200.0)] = 3200.0
	#Tc[np.where(T>50000.0)] = 50000.0
	z = np.arange(zmax-zmin) * 16.0
	z3d = np.zeros([NX, NY, NZ])
	z3d[:,:,:] = z[None,None,:]
	tau = np.linspace(-6,2,NZ) # Should not matter
	p = snap.Pres[zmin:zmax, xmin:xmax:skip, ymin:ymax:skip].transpose(1,2,0)
	vz = snap.vx[zmin:zmax, xmin:xmax:skip, ymin:ymax:skip].transpose(1,2,0)

	atmout.set_z(z3d)
	atmout.set_pg(p)
	atmout.set_tem(Tc)
	atmout.set_vz(-vz)

	Bz = snap.Bx[zmin:zmax, xmin:xmax:skip, ymin:ymax:skip].transpose(1,2,0) * np.sqrt(4.0*np.pi)
	Bx = snap.By[zmin:zmax, xmin:xmax:skip, ymin:ymax:skip].transpose(1,2,0) * np.sqrt(4.0*np.pi)
	By = snap.Bz[zmin:zmax, xmin:xmax:skip, ymin:ymax:skip].transpose(1,2,0) * np.sqrt(4.0*np.pi)
	
	atmout.set_bx(Bx)
	atmout.set_by(By)
	atmout.set_bz(Bz)
	
	outputname = sys.argv[11] + '_' + sys.argv[2]+ '.bin'

	atmout.write_model(outputname)	


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

	atmout = frz.atm_model3D(NX, NY, NZ)
	
	Tc = np.copy(T)
	Tc[np.where(T<3200.0)] = 3200.0
	#Tc[np.where(T>50000.0)] = 50000.0
	z = np.arange(zmax-zmin) * 16.0
	z3d = np.zeros([NX, NY, NZ])
	z3d[:,:,:] = z[None,None,:]
	tau = np.linspace(-6,2,NZ) # Should not matter
	p = snap.Pres[zmin:zmax, xmin:xmax:skip, ymin:ymax:skip].transpose(1,2,0)
	vz = snap.vx[zmin:zmax, xmin:xmax:skip, ymin:ymax:skip].transpose(1,2,0)

	atmout.set_z(z3d)
	atmout.set_pg(p)
	atmout.set_tem(Tc)
	atmout.set_vz(-vz)

	Bz = snap.Bx[zmin:zmax, xmin:xmax:skip, ymin:ymax:skip].transpose(0,2,1) * np.sqrt(4.0*np.pi)
	Bx = snap.By[zmin:zmax, xmin:xmax:skip, ymin:ymax:skip].transpose(0,2,1) * np.sqrt(4.0*np.pi)
	By = snap.Bz[zmin:zmax, xmin:xmax:skip, ymin:ymax:skip].transpose(0,2,1) * np.sqrt(4.0*np.pi)
	
	atmout.set_bx(Bx)
	atmout.set_by(By)
	atmout.set_bz(Bz)
	
	outputname = sys.argv[11] + '_' + sys.argv[2]+ '.bin'

	atmout.write_model(outputname)

else:
        print ("Uknown file type. Exiting...")

