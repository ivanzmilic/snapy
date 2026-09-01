# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo is

`snapy` is a personal collection of solar-physics data-wrangling scripts (not a packaged library). Its central job is converting 3D atmospheric cubes from various radiative-MHD simulation codes and instruments into the **SNAPI `.f0` atmosphere format** (and back), plus assorted spectral-synthesis, reinterpolation, and visualization helpers. Scripts are run standalone from the command line; there is no build system, test suite, linter, or package manifest.

## Running things

- Everything runs under the user's miniconda Python: `python <script>.py <args...>`. There is no `requirements.txt` — dependencies are whatever is installed in that conda env.
- Scripts take **positional `sys.argv` arguments** (no argparse). To learn a script's interface, read the `sys.argv[...]` assignments near the top; several scripts (e.g. `muram_npy_to_snapi.py`) also print a `usage:` line. Argument order and meaning are not standardized across scripts — always check the specific file.
- There are no automated tests. "Testing" a change means running the script on a real cube and sanity-checking the output (this is what the `PLOT_MAPS`-style switches in `cobold_to_snapi.py` are for).

### Key external/local dependencies
- **`pyana`** — external package (installed in the conda env, *not* in this repo). Provides `pyana.fzread(name)["data"]` and `pyana.fzwrite(name, array, 0, 'comment')` for reading/writing ANA `.f0`/`.f0`-style binary cubes. This is the SNAPI on-disk format.
- **`muram.py`** (local) — reader for MURaM simulation output. Key classes: `MuramSnap`, `MuramSubSnap` (full/sub cubes), `MuramSlice`/`MuramTauSlice`, `MuramCube`. Cube objects expose named attributes: `.Temp`, `.Pres`, `.rho`, `.vx/.vy/.vz`, `.Bx/.By/.Bz` (B already multiplied by `sqrt(4*pi)` → Gauss), `.tau`.
- **`firtez_dz.py`** (local) — reader/writer for FIRTEZ-dz models and Stokes profiles (`frz.read_model`, Fortran-record I/O via `scipy.io.FortranFile`).
- **`lightweaver`** (external) — used in `mini_lw_synth.py` and the `*_synth*` notebooks for actual spectral synthesis.
- Also common: `numpy`, `scipy`, `matplotlib`, `h5py` (CO5BOLD input), `astropy.io.fits`.

## The SNAPI `.f0` atmosphere format (the core data structure)

This is the shared contract that ties most scripts together. A SNAPI atmosphere is a float array of shape **`[12, NX, NY, NZ]`** where the first axis is the physical quantity:

| idx | quantity | notes |
|----|----------|-------|
| 0 | `log(tau)` | independent depth variable; usually a generic placeholder `np.linspace(-6, 2, NZ)` (SNAPI recomputes real optical depth from the structure) |
| 1 | `z` geometric height [cm] | often `arange(NZ) * dz`, dz given in km then ×1e5 |
| 2 | `T` [K] | **floored at 3200 K** (SNAPI opacity floor) |
| 3 | `p` gas pressure [dyn/cm²] | interpolate in `log10` space, not linear |
| 4 | `pe` electron pressure [dyn/cm²] | crude guess `≈ 0.05 * p` when not available; refined by SNAPI |
| 5, 6 | unused (0) | |
| 7 | `|B|` [G] | field magnitude |
| 8 | microturbulence | `firtez_to_snapi.py` stores `sqrt(vmic)`; otherwise 0 |
| 9 | `v_LOS` / `vz` [cm/s] | vertical velocity |
| 10 | `theta` [rad] | field inclination from vertical, `arccos(Bz / (|B| + eps))` |
| 11 | `phi` [rad] | field azimuth, `arctan2(By, Bx)` |

### Conventions that recur and matter
- **z-axis reversal on write.** Simulation cubes usually have z-index 0 = deep/bottom; SNAPI wants index 0 = top of atmosphere. Scripts therefore write `atmout[:, :, :, ::-1]`. `pyana.fzwrite` can **segfault on a large negative-stride view**, so pass a contiguous copy: `np.ascontiguousarray(atmout[..., ::-1])` (see `muram_npy_to_snapi.py`).
- **Divide-by-zero guard** in the field geometry: always `arccos(Bz / (|B| + eps))` with a small eps so field-free points don't blow up.
- **MURaM coordinate remap.** In `muram_to_snapi.py` the MURaM loader's `x` axis is the vertical: it maps `snap.Bx→Bz`, `snap.By→Bx`, `snap.Bz→By`, and `vz = snap.vx`. Don't assume MURaM's named components line up with SNAPI's without checking the specific `type` branch.
- MURaM B fields already include the `sqrt(4*pi)` cgs factor (applied in `muram.py`); scripts reading raw MURaM sometimes re-apply it — check before double-counting.

## Script families

- **`*_to_snapi.py`** — the main converters into SNAPI `.f0`: `cobold_to_snapi.py` (CO5BOLD h5py), `muram_to_snapi.py` (MURaM binary, with `type` = `muram`/`muramsub`/`muramt` branches for different transposes), `muram_npy_to_snapi.py` (pre-extracted `.npy` cubes, auto-detects 8- vs 9-channel and quantity-first vs -last layout), `firtez_to_snapi.py` (FIRTEZ, with Gaussian-smoothing modes), `np_to_snapi.py` (Shah's MURaM `.npy`, B ignored), `mihi_to_snapi.py` (MiHi FITS Stokes cubes → normalized observations + wavelength/mask `.dat`).
- **`muram_to_cube.py`, `muram_to_firtez.py`, `muram_dat_reading.py`, `make_iso_tau_cube.py`** — MURaM readers that emit other formats / iso-tau slices.
- **SNAPI cube manipulation** — `reinterpolate_snapi_atmosphere.py` (re-grid onto a new tau grid, cubic 1D interp, pressure in log space), `rotate_snapi_atm.py` (accurate 3D rotation via `RegularGridInterpolator`, y–z plane at a time), `f0_to_dat.py` (dump a column to text), `snapi_pops_to_firtez_dcs.py` (departure coefficients).
- **Synthesis / training** — `mini_lw_synth.py` (single-line Lightweaver synthesis callable from notebooks), `prep_training.py` (SVD compression + packing cubes for ML), `prepare_lgrid.py`, `pyana_to_fits_series_spectra.py`.
- **Notebooks (`.ipynb`)** — exploratory visualization and prototyping (MURaM 2D/3D, CO5BOLD, tilting cubes, spectra inspection). Treat these as scratch/analysis, not reusable modules.

## Working style in this repo
- Scripts are terse and evolve by copy-paste-and-tweak between the `*_to_snapi` variants. When adding a converter for a new input format, start from the closest existing one and preserve the 12-quantity layout, the flooring/eps guards, and the z-reversal-on-write.
- Hardcoded cuts and magic indices (e.g. a z-range crop like `[..., 105:]`, downsample `skip` factors) are common and expected; they're per-dataset knobs, not bugs.
