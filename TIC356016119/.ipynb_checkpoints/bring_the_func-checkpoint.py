from astropy.modeling import models, fitting, optimizers, statistic, custom_model
import numpy as np
import matplotlib.pyplot as plt
import lightkurve as lk
import emcee
import batman
from lightkurve.lightcurve import LightCurve
from ldtk import (LDPSetCreator,BoxcarFilter)
from ldtk.filters import kepler
import corner as triangle

def Limb_Dark(Teff = 3300,
              Terr = 100,
              log_g = 4.4,
              g_err = 0.1,
              met = 0.0,
              met_err = 0.02,
              filters = [kepler]):
        
    sc = LDPSetCreator(filters=filters,
                   teff=[Teff,Terr],
                   logg=[log_g, g_err],
                   z=[met, met_err])
    
    ps = sc.create_profiles(nsamples=500)
    qc,qe = ps.coeffs_qd(do_mc=True) 
    
    return qc,qe

def BATMAN(Rp,
           t0,
           LD,
           t = None,
           Baseline = 1.0):
    
    
                    
    params = batman.TransitParams()
    params.t0 = t0                       # time of inferior conjunction ()
    params.per = 18.7913430                       # period in days
    params.rp = Rp                       # planet radius (in units of stellar radii)
    params.a = 56.8                         # semi-major axis (in units of stellar radii)
    params.inc = 89.1                     # orbital inclination (in degrees)
    params.ecc = 0.                      # eccentricity
    params.w = 90.                       # longitude of periastron (in degrees)
    params.u = LD                        # limb darkening coefficients [u1, u2]
    params.limb_dark = "quadratic"       # limb darkening model
        
    m = batman.TransitModel(params, t)   # initializes model
    
    flux = m.light_curve(params)*Baseline   # calculates light curve
    
    return flux


def lnprior(lc,
            theta,
            expected_t0):
    
    C1, C2, rp, t_0 = theta
        
    if (0.0<=rp<= 0.1)and(np.min(lc.time)-0.015<=t_0<=np.max(lc.time)):
        
        rp_mean = 0.0375
        rp_sigma = 0.0246
        t0_mean = 2367.818618
        t0_sigma = 0.0139 #days, calculated from propagating the period error

        
        lnprior_rp = - (1.0/2.0)*((rp-rp_mean)/rp_sigma)**2.0
#         lnprior_a = - (1.0/2.0)*((a-a_mean)/a_sigma)**2.0
#         lnprior_i = - (1.0/2.0)*((i-i_mean)/i_sigma)**2.0

        lnprior_t0 = - (1.0/2.0)*((t_0-t0_mean)/t0_sigma)**2.0
        
        return lnprior_t0 + lnprior_rp
    
    return -np.inf


def lnprob(theta,
           lc,
           planet_period,
           LD,
           airmass,
           planet_a,
           expected_t0,
           plot = False):
        
    # Pull out some model parameters
    C1, C2, rp, t_0 = theta
        
    # First we want a model to perform the lnprob calculation with.
    model = (C1 + C2*(airmass-1))*BATMAN(rp, t_0, LD, t = lc.time)
    
    if plot:
        
        hires_times = np.linspace(np.min(lc.time),np.max(lc.time),1000)
        model_to_plot = BATMAN(rp, t_0, LD, t = hires_times)
        
        plt.figure()
        plt.errorbar(lc.time,lc.flux/(C1 + C2*(airmass-1)),yerr=lc.flux_err/(C1 + C2*(airmass-1)),
                    fmt='o',alpha=0.5,color='royalblue',markersize='5')
        plt.plot(hires_times,model_to_plot,label='Best-Fit Model',color='k',zorder=100)
        plt.ylabel('Normalized Flux',fontsize=18)
        plt.xlabel('BJD',fontsize=16)
        
    # We need to make sure the uniform priors are accounted for:
    ln_prior = lnprior(lc, theta, expected_t0)
    
    # This is a Gaussian likelihood, for independent data points
    
    chisq = np.sum((lc.flux - model)**2/(lc.flux_err)**2)
    ln_like = (np.sum(1/np.sqrt(2*np.pi*(lc.flux_err))) - 0.5*chisq)
    
    return ln_prior + ln_like
    
def plot_chain(sampler,
               static_params,
               start=0,
               stop=0):
    
    planet_radius, planet_period, planet_a, planet_i, expected_t0, LD, offset, observatory, airmass, toi = static_params

    C1_all, C2_all, Rp_all, t0_all = sampler.chain.T
    
    ndim = np.shape(sampler.chain.T[:,0,0])[0]
    
    C1, C2, Rp, t0 = sampler.chain[:, start:, :].reshape((-1, ndim)).T
    
    C1_burn, C2_burn, Rp_burn, t0_burn = sampler.chain[:, :start, :].reshape((-1, ndim)).T
    
    plt.figure(figsize=(14,14))
    gs = plt.matplotlib.gridspec.GridSpec(2,2)
    
    #Walker Plots
    
    ax_Rp = plt.subplot(gs[0,0])
    ax_Rp.plot(Rp_all.flatten(),color='black',alpha=0.5)
    ax_Rp.axvspan(start, stop, zorder=-1,alpha=0.3)
    ax_Rp.set_ylabel('Rp')
    
    ax_t0 = plt.subplot(gs[0,1])
    ax_t0.plot(t0_all.flatten(),color='black',alpha=0.5)
    ax_t0.axvspan(start, stop, zorder=-1,alpha=0.3)
    ax_t0.set_ylabel('t0')
    
#     ax_A = plt.subplot(gs[1,0])
#     ax_A.plot(A_all.flatten(),color='black',alpha=0.5)
#     ax_A.axvspan(start, stop, zorder=-1,alpha=0.3)
#     ax_A.set_ylabel('a/R*')
    
#     ax_i = plt.subplot(gs[1,1])
#     ax_i.plot(I_all.flatten(),color='black',alpha=0.5)
#     ax_i.axvspan(start, stop, zorder=-1,alpha=0.3)
#     ax_i.set_ylabel('Inclination')
    
    #Histograms
    
    ax_Rphist = plt.subplot(gs[1,0])
    ax_Rphist.hist(Rp.flatten(),color='black',bins=100)
    ax_Rphist.axvline(planet_radius,zorder=100,color='purple',label='Expected Value')
    ax_Rphist.set_xlabel('Rp')
    ax_Rphist.legend()
    
    ax_t0hist = plt.subplot(gs[1,1])
    ax_t0hist.hist(t0.flatten(),color='black',bins=100)
    ax_t0hist.axvline(expected_t0,zorder=100,color='purple',label='Expected Value')
    ax_t0hist.set_xlabel('t0')
    ax_t0hist.legend()
    
    plt.savefig(f'Figs/TOI_{toi}_plotchain.pdf')
    
def corner(samples,
           labels,
           toi = 'no_toi_given'):
        
    Samples = samples.T
    
    np.random.seed(68)
    figure = triangle.corner(Samples,labels=labels,
                       quantiles=[0.16, 0.5, 0.84],
                       show_titles=True, title_kwargs={"fontsize": 12})
    plt.savefig(f'Figs/TOI_{toi}_Corner.pdf')
    
def light_curve(lc,
                parameters,
                sampler,
                static_params,
                nsteps,
                epoch,
                tic_id):
    
    ########################
    # Read in the parameters
    ########################
    
    planet_radius, planet_period, planet_a, planet_i, expected_t0, LD, offset, observatory, airmass, toi = static_params
    hires_times = np.linspace(np.min(lc.time),np.max(lc.time),1000) #This array is for the model that gets plotted
    C1_best, C2_best, rp_best, t_0_best = parameters
    ndim = np.shape(sampler.chain.T[:,0,0])[0]
    burnin = int(0.25*nsteps)
    C1, C2, Rp, t0 = sampler.chain[:, burnin:, :].reshape((-1, ndim)).T
    
    ########################################################################
    # Calculate the model with the best-fit MCMC params, calculate residuals
    ########################################################################
    
    model_to_plot = BATMAN(rp_best, t_0_best, LD, t = hires_times)
    best_model = (C1_best + C2_best*(airmass-1))*BATMAN(rp_best, t_0_best,LD, t = lc.time)
    residual = (lc.flux-best_model)/lc.flux_err
    
    #########################
    # Plot the best-fit model
    #########################
            
    f, (a0, a1) = plt.subplots(2,1, gridspec_kw = {'height_ratios':[4,1]},
                               figsize=(14,7),sharex=True)
    
    a0.set_title(f'Modeled LCO Light Curve for TOI {toi}',fontsize=20)
    a0.errorbar(lc.time,lc.flux/(C1_best + C2_best*(airmass-1)),
                yerr=lc.flux_err/(C1_best + C2_best*(airmass-1)),
                fmt='o',alpha=0.5,color='royalblue',markersize='5',
                label=observatory+' data')
    a0.plot(hires_times,model_to_plot,label='Best-Fit Model',color='k',zorder=100)
    a0.set_ylabel('Normalized Flux',fontsize=18)
    a0.minorticks_on()
    a0.legend(loc='lower right')
    a0.tick_params(labelbottom=False, labeltop=False, labelleft=True, labelright=False,
                   bottom=True, top=True, left=True, right=True)
    
    #########################################
    # Plot 200 random models from the samples
    #########################################
    
    for j in range(0,200,1):
        i = np.random.randint(low=0,high=(nsteps-burnin)*100)
        sigma_model = BATMAN(Rp = Rp[i], t0 = t0[i],LD = LD, t = hires_times)
        a0.plot(hires_times,sigma_model,color='red',alpha = 0.1,
                linewidth=0.8,zorder=-1000,label='Random Samples')
        
    ####################
    # Plot the residuals
    ####################

    a1.scatter(lc.time,residual,color='royalblue',alpha=0.5)
    a1.axhline(0,color='k')
    a1.set_ylim(0-1.5*np.max(np.abs(residual)),0+1.5*np.max(np.abs(residual)))
    a1.minorticks_on()
    a1.set_ylabel(r'Residuals ($\rm{\sigma}$)',fontsize=15)
    a1.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False,
                   bottom=True, top=True, left=True, right=True)
        
    plt.xlabel(f'BJD - {offset}',fontsize=16)
    plt.savefig(f'Figs/TOI_{toi}_modeled.pdf')
    
def rms_plot(lc,
             parameters,
             static_params):
    
    C1_best, C2_best, rp_best, t_0_best = parameters
    planet_radius, planet_period, planet_a, planet_i, expected_t0, LD, offset, observatory, airmass, toi = static_params
        
    lc_fixed = lc/(C1_best + C2_best*(airmass-1))
    
    best_model = BATMAN(rp_best, t_0_best, LD, t = lc_fixed.time)
    
    residual = (lc_fixed.flux-best_model)
    residual_lc = LightCurve(time = lc_fixed.time, flux = residual)
    
    time_range = (np.amax(lc.time)-np.amin(lc.time))*24.0*60.0
    min_per_exp = time_range/len(lc.time)
    
    bin_size = []
    n = []
    RMS = []
    ERR_RMS = []
    
    for N in range(1,20,1):
        
        i = N*min_per_exp
        
        LC = residual_lc.remove_outliers(sigma=5).bin(binsize=int(N))
        
        rms = np.std(LC.flux)
        err = rms/np.sqrt(2*(len(LC.flux)-1))
        
        if int(10/min_per_exp) == N:
            print('RMS at 10 min binning = {:.7f}'.format(rms))
            print('Error at 10 min binning = {:.7f}'.format(err))
            
        n.append(N)
        RMS.append(rms)
        ERR_RMS.append(err)
        bin_size.append(i)
            
    plt.errorbar(bin_size,RMS,yerr=ERR_RMS)
    plt.plot(bin_size,RMS[0]/np.sqrt(n))
    plt.xlabel('Bin Size (min)')
    plt.ylabel('RMS')
    plt.yscale('log')
    plt.xscale('log')
    plt.title(f'RMS vs Bins for TOI {toi}')