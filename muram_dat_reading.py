

#Yeah, you get the information from the cubes:
header = np.loadtxt('3D/Header.024000')
nx = np.int64(header[0])
ny = np.int64(header[1])
nz = np.int64(header[2])
dx = np.float64(header[3])
dy = np.float64(header[4])
dz = np.float64(header[5])
time= np.float64(header[6])
delta_t = np.float64(header[7])
maxva= np.float64(header[8])

rho = np.memmap('3D/result_prim_0.024000',dtype=np.float32,mode='r',shape=(nz,ny,nx))

#Then the cubes are all float32 binaries:

#They are all float32 binaries:
result_prim_0 = density (g/cm^3)
result_prim_1 = vx (vertical) (cm/s)
result_prim_2 = vy (horziontal) (cm/s)
result_prim_3 = vz (horizontal) (cm/s)
result_prim_4 = energy (erg/cm^3)
result_prim_5 = Bx (vertical) 1/sqrt(4pi) G
result_prim_6 = By (horziontal) 1/sqrt(4pi) G
result_prim_7 = Bz (horizontal) 1/sqrt(4pi) G
eosT = Temperature (K)
eosP = pressure
tau = continuum optical depth, 10A band @ 500nm
Qres = resistive heating (erg/cm^3/s)
Qvis = viscous heating (erg/cm^3/s)

Slices are:
nx = 1536 
ny = 1536
nz = 1056
Lx = 10e8
Ly = 10e8
Lz =  6.875e8

time = np.fromfile('2D/Pore_10Mm_6x6km_res_Bz400G_time_012000_to_044450.dat',dtype=np.float64)

# the file contains the iteration number and simulation solar time [s]
nsnaps = time.size//2
time = time.reshape(2,nsnaps)

Iout = np.memmap('2D/Pore_10Mm_6x6km_res_Bz400G_Iout_bref_012000_to_044450.dat',dtype=np.float32,mode='r',shape=(nsnaps,nx,ny))
