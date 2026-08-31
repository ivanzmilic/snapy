import numpy as np
import pyana
import sys
import os

# muram_npy_to_snapi.py
# -----------------------------------------------------------------------------
# Convert a MURaM-derived atmospheric .npy cube (as produced by muram_to_cube.py
# or downstream processing) into a SNAPI .f0 atmosphere.
#
# The .npy cubes carry only the physical MHD variables -- they are MISSING the
# optical-depth axis and the geometric-height axis that SNAPI expects. This
# script adds both:
#   * tau : a GENERIC placeholder, linspace(-6, 2, NZ). SNAPI recomputes the
#           real optical depth from the physical structure, so the values here
#           are only a monotonic filler (matches muram_to_snapi.py exactly).
#   * z   : the geometric height, arange(NZ) * delta_z, where delta_z (the
#           vertical grid spacing, in km) is supplied on the command line.
#
# Two input layouts are auto-detected (both seen in the downloaded cubes):
#   Format A  [9, NX, NY, NZ] = [T, p, rho, vx, vy, vz, Bx, By, Bz]  (quantity axis FIRST)
#   Format B  [NX, NY, NZ, 8] = [T, p,      vx, vy, vz, Bx, By, Bz]  (quantity axis LAST)
# In both layouts z-index 0 is the DEEP/bottom layer and the last z-index is the
# top of the atmosphere (verified from the pressure stratification), and the
# magnetic field is already in Gauss (the sqrt(4*pi) factor is baked in upstream).
#
# The SNAPI .f0 cube written out is [12, NX, NY, NZ]:
#    0: log(tau)      generic placeholder
#    1: z   [cm]      geometric height (from delta_z)
#    2: T   [K]       floored at 3200 K
#    3: p   [dyn/cm2] gas pressure (floored at a tiny positive value)
#    4: pe  [dyn/cm2] electron pressure ~ 0.05 * p  (crude guess, refined by SNAPI)
#    5,6: unused (0)
#    7: |B| [G]
#    8: unused (0)
#    9: v_LOS [cm/s]  vertical velocity
#   10: theta [rad]   field inclination from the vertical (Bz)
#   11: phi   [rad]   field azimuth
# It is written with the z-axis reversed (index 0 -> top of atmosphere), which is
# the SNAPI convention used by muram_to_snapi.py and by the existing *.f0 files.
#
# Usage:
#   python muram_npy_to_snapi.py <input.npy> <delta_z_km> [output.f0]
# Example:
#   python muram_npy_to_snapi.py cut_phch_50.npy 16.0
#   python muram_npy_to_snapi.py 10_sloja_no_phys.npy 12.0 sloja_no_phys.f0
# -----------------------------------------------------------------------------

if len(sys.argv) < 3:
    print("usage: python muram_npy_to_snapi.py <input.npy> <delta_z_km> [output.f0]")
    sys.exit(1)

infile = sys.argv[1]
delta_z_km = float(sys.argv[2])                       # vertical grid spacing [km]
if len(sys.argv) > 3:
    outfile = sys.argv[3]
else:
    outfile = os.path.splitext(infile)[0] + '.f0'     # foo.npy -> foo.f0

dz_cm = delta_z_km * 1e5                               # km -> cm

print("info :: reading", infile)
atmin = np.load(infile)
if atmin.ndim != 4:
    print("error :: expected a 4D cube, got shape", atmin.shape)
    sys.exit(1)
shp = atmin.shape

# --- auto-detect the layout and move the quantity axis to the front ----------
if shp[0] in (8, 9):                    # quantity axis first  (Format A)
    Q = shp[0]
    cube = atmin                                       # [Q, NX, NY, NZ]
    print("info :: detected quantity-first layout, shape", shp)
elif shp[-1] in (8, 9):                 # quantity axis last   (Format B)
    Q = shp[-1]
    cube = np.moveaxis(atmin, -1, 0)                   # -> [Q, NX, NY, NZ]
    print("info :: detected quantity-last layout, shape", shp)
else:
    print("error :: cannot find a quantity axis of length 8 or 9 in shape", shp)
    sys.exit(1)

# --- map channels to physical variables (see header) -------------------------
# vz is the VERTICAL/line-of-sight velocity; Bz is the vertical field component.
if Q == 9:      # [T, p, rho, vx, vy, vz, Bx, By, Bz]
    T  = cube[0]; p = cube[1]
    vz = cube[5]
    Bx = cube[6]; By = cube[7]; Bz = cube[8]
    print("info :: 9-channel cube [T, p, rho, vx, vy, vz, Bx, By, Bz]")
else:           # Q == 8 : [T, p, vx, vy, vz, Bx, By, Bz]
    T  = cube[0]; p = cube[1]
    vz = cube[4]
    Bx = cube[5]; By = cube[6]; Bz = cube[7]
    print("info :: 8-channel cube [T, p, vx, vy, vz, Bx, By, Bz] (no density)")

NX, NY, NZ = T.shape
print("info :: horizontal x vertical dimensions:", NX, NY, NZ)

# --- build the generic tau axis and the geometric z axis ---------------------
tau = np.linspace(-6.0, 2.0, NZ)                       # generic placeholder
z   = np.arange(NZ) * dz_cm                            # z[0]=0 at the deep layer
print("info :: delta_z = %.3f km  ->  z spans 0 .. %.4e cm (%.1f km)"
      % (delta_z_km, z[-1], z[-1] / 1e5))

# --- clean up the physical fields (kill interpolation artifacts) -------------
T = np.array(T, dtype=np.float64)
p = np.array(p, dtype=np.float64)

print("info :: T  range before flooring: %.1f .. %.1f K" % (T.min(), T.max()))
n_coldT = int(np.count_nonzero(T < 3200.0))
T[T < 3200.0] = 3200.0                                 # SNAPI opacity floor
if n_coldT:
    print("info :: floored %d points with T < 3200 K" % n_coldT)

print("info :: p  range before flooring: %.3e .. %.3e" % (p.min(), p.max()))
n_lowp = int(np.count_nonzero(p < 1e-4))
p[p < 1e-4] = 1e-4                                     # guard vs negatives/zeros
if n_lowp:
    print("info :: floored %d points with p < 1e-4 (negative/zero artifacts)" % n_lowp)

# --- magnetic field: magnitude, inclination (from vertical Bz), azimuth -------
B     = np.sqrt(Bx**2.0 + By**2.0 + Bz**2.0)
theta = np.arccos(Bz / (B + 1e-3))                     # eps avoids 0/0 in field-free
phi   = np.arctan2(By, Bx)

# --- assemble the 12-quantity SNAPI cube -------------------------------------
atmout = np.zeros([12, NX, NY, NZ])
atmout[0]  = tau[None, None, :]
atmout[1]  = z[None, None, :]
atmout[2]  = T
atmout[3]  = p
atmout[4]  = p * 0.05
atmout[7]  = B
atmout[9]  = vz
atmout[10] = theta
atmout[11] = phi

# --- write, reversing the z-axis so index 0 -> top of atmosphere -------------
# NOTE: pyana.fzwrite segfaults on a large negative-stride VIEW, so we hand it a
# contiguous copy of the reversed cube (the reversal is identical either way).
print("info :: writing", outfile, " shape", atmout.shape)
pyana.fzwrite(outfile, np.ascontiguousarray(atmout[:, :, :, ::-1]), 0, 'muram_npy_to_snapi')
print("info :: done.")
