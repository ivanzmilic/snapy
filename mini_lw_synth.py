# A piece of code that outsources a synthesis of a single line to a separate function
# which can be called from the notebook.

import numpy as np
import lightweaver as lw 
from lightweaver.rh_atoms import (Al_atom, C_atom, CaII_atom, Fe_simple_atom, H_6_atom, He_9_atom, MgI_atom, N_atom, Na_atom, O_atom, S_atom, Si_atom)


def airtovac(lambda_air):

    s = 1E2/(lambda_air*10.0);
    n = 1.0 + 0.00008336624212083 + 0.02408926869968 / (130.1065924522 - s*s) + 0.0001599740894897 / (38.92568793293 - s*s);
    return lambda_air * n

def synth(atmos, conserve, prd, stokes, wave, mu, actives, lte=False):
    
    '''
    Synthesize a spectral region with given parameters:
    
    Parameters
    ----------
    atmos : lw.Atmosphere - The atmospheric model in which to synthesise the line.
    
    conserve : bool - Whether to start from LTE electron density and conserve charge, or simply use from the electron density present in the atmospheric model.

    prd: bool - whether to use prd or no, most of the time it's no 

    stokes: bool - whether to synth all 4 Stokes parameters or I only - for tracking we start with I only
    
    wave : np.ndarray Array of wavelengths over which to resynthesise the final line profile

    mu : mu angle, gonna use 1.0 most of the times 

    actives: list of active species to synthesize

    Returns
    -------
    ctx : lw.Context -The Context object that was used to compute the equilibrium populations -> Gonna not return this
    
    Iwave : np.ndarray - The intensity at given mu and wave    '''
    
    # Configure the atmospheric angular quadrature - only matters for NLTE. Gonna use 3 for faster calc
    atmos.quadrature(3)#, force3d=True)
    # See if you can force this to 
    # Replace this with atmos.rays ( specify mu ) - let's think how to use Stokes with it
    # ctx.single_stokes_fs
    
    # Configure the set of atomic models to use. Contrary to SNAPI you have to explicitly specify all species
    # Annoying, but since you have all the atoms in the file - that is fine.
    aSet = lw.RadiativeSet([H_6_atom(), C_atom(), O_atom(), Si_atom(), Al_atom(), CaII_atom(), Fe_simple_atom(), He_9_atom(), MgI_atom(), N_atom(), Na_atom(), S_atom()])
    
    # Set actives to the ones you want
    #aSet.set_active(actives)
    aSet.set_active(actives)
    
    # Compute the necessary wavelength dependent information (SpectrumConfiguration).
    spect = aSet.compute_wavelength_grid()
    
    # Calculate electron density in lte, we are never using the electron density from the model
    eqPops = aSet.iterate_lte_ne_eq_pops(atmos, direct=True)

    # Configure the Context which holds the state of the simulation for the  backend, and provides 
    # the python interface to the backend.
    # Feel free to increase Nthreads to increase the number of threads the program will use.
    # I would always stick to Nthreads = 1 as we are looking to mpi this one
    
    ctx = lw.Context(atmos, spect, eqPops, conserveCharge=False, Nthreads=1, formalSolver='piecewise_linear_1d')
    
    # Iterate the Context to convergence (using the iteration function now
    # provided by Lightweaver). Go test this one in order to calculate stuff in LTE!
    
    #lw.iterate_ctx_se(ctx, prd=prd)
    
    # NLTE quiet
    if (lte == False):
        lw.iterate_ctx_se(ctx, prd=prd, quiet=True)

    # LTE:
    else:
        ctx.formal_sol_gamma_matrices()

    
    # Update the background populations based on the converged solution and
    # compute the final intensity for mu=1 on the provided wavelength grid.
    #eqPops.update_lte_atoms_Hmin_pops(atmos)
    
    # Calculate the (stokes) spectru, at the provided wavegrid at the specified mu
    
    Iwave = ctx.compute_rays(wave, [mu], stokes=stokes) 

    # We will want to return some populations or so, at some point (Firtez pipeline)
    #return ctx, Iwave

    # For now we are only returning the intensity:
    if (stokes == False):
        Iwave = Iwave.reshape(1,-1)
    
    return ctx, Iwave


def synth_final(z_scale, atmos, wave, t_wave, lte=False):
    
    ND = atmos.shape[1]
    wavevac = airtovac(wave)
    
    atmos = lw.Atmosphere.make_1d(scale=lw.ScaleType.Geometric, depthScale=z_scale[:], temperature=atmos[0, :],  \
        vlos=atmos[3,:]/1E2, vturb=np.ones(ND)*0.0, Pgas=atmos[1,:]/10.0, convertScales=False)

    ctxtemp, Itemp = synth(atmos, conserve=False, prd=False, stokes=False, wave=wave, mu=1.0, actives='Na', lte=lte)
    
    return ctxtemp, Itemp
    
