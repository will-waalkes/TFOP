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
              met_err = 0.01,
              filters = [kepler]):
        
    sc = LDPSetCreator(filters=filters,
                   teff=[Teff,Terr],
                   logg=[log_g, g_err],
                   z=[met, met_err])
    
    ps = sc.create_profiles(nsamples=500)
    qc,qe = ps.coeffs_qd(do_mc=True) 
    
    return qc,qe

def BATMAN(P,
           Rp,
           t0,
           inc,
           A,
           LD,
           t = None,
           Baseline = 1.0):
    
    
                    
    params = batman.TransitParams()
    params.t0 = t0                       # time of inferior conjunction ()
    params.per = P                       # period in days
    params.rp = Rp                       # planet radius (in units of stellar radii)
    params.a = A                         # semi-major axis (in units of stellar radii)
    params.inc = inc                     # orbital inclination (in degrees)
    params.ecc = 0.                      # eccentricity
    params.w = 90.                       # longitude of periastron (in degrees)
    params.u = LD                        # limb darkening coefficients [u1, u2]
    params.limb_dark = "quadratic"       # limb darkening model
        
    #m = batman.TransitModel(params, t, exp_time=60.0, supersample_factor = 12)   # initializes model
    m = batman.TransitModel(params, t)   # initializes model
    
    flux = m.light_curve(params)*Baseline   # calculates light curve
    
    return flux

# def lnprior(lc,
#             astropy_model):
        
#     i_min = np.arccos(1/astropy_model.a)*180/np.pi
    
#     n = len(astropy_model.param_names)

#     i = 0
#     for k in astropy_model.param_names:
#         if astropy_model.fixed[k] == False:
#             for m in range(n):
#                 if astropy_model.param_names[m] == str(k):
#                     astropy_model.parameters[m] = params[i]
#                     i += 1

#     if (astropy_model.a <= astropy_model.b):
#         astropy_model.a = astropy_model.b + 0.01
    
#     model = astropy_model(lc.time)

#     if (0.0 <= astropy_model.radius <= 1.0) and (lc.time[0] <= astropy_model.t0 <= lc.time[-1] ) and (1.0 <= astropy_model.a <= 200.0 ) and (i_min <= astropy_model.i <= 90.0):
        
#         return 0.0
    
#     return -np.inf

def lnprior(lc,
            theta,
            planet_a,
            expected_t0):
    
    C1, C2, rp, t_0, a, i = theta
        
    i_min = np.arccos(1/a)*180/np.pi

    if (0.0<=rp<= 1.0)and(expected_t0-0.07<=t_0<=expected_t0+0.07)and(2<=a<=planet_a+10)and(70.0<=i<= 90.0):
        
#         rp_mean = 0.0533
#         rp_sigma = 0.004
#         a_mean = 23.93
#         a_sigma = 2.04
#         i_mean = 89.29
#         i_sigma = 0.51
        
#         lnprior_rp = - (1.0/2.0)*((rp-rp_mean)/rp_sigma)**2.0
#         lnprior_a = - (1.0/2.0)*((a-a_mean)/a_sigma)**2.0
#         lnprior_i = - (1.0/2.0)*((i-i_mean)/i_sigma)**2.0
        
        return 0.0 #lnprior_rp + lnprior_a + lnprior_i
    
    return -np.inf


    

# def lnprob(astropy_model,
#            lc,
#            planet_period,
#            LD,
#            airmass,
#            plot = False):
        
#     # Pull out some model parameters
#     C1, C2, rp, t_0, a, i = theta
        
#     # First we want a model to perform the lnprob calculation with.
#     model = (C1 + C2*(airmass-1))*BATMAN(planet_period, rp, t_0, i, a, LD, t = lc.time)
    
#     if plot:
        
#         hires_times = np.linspace(lc.time[0],lc.time[-1],1000)
#         model_to_plot = BATMAN(planet_period, rp, t_0, i, a, LD, t = hires_times)
        
#         plt.figure()
#         plt.errorbar(lc.time,lc.flux/(C1 + C2*(airmass-1)),yerr=lc.flux_err/(C1 + C2*(airmass-1)),
#                     fmt='o',alpha=0.5,color='royalblue',markersize='5')
#         plt.plot(hires_times,model_to_plot,label='Best-Fit Model',color='k',zorder=100)
#         plt.ylabel('Normalized Flux',fontsize=18)
#         plt.xlabel('BJD',fontsize=16)
        
#     # We need to make sure the uniform priors are accounted for:
#     ln_prior = lnprior(lc, astropy_model)
    
#     # This is a Gaussian likelihood, for independent data points
    
#     chisq = np.sum((lc.flux - model)**2/(lc.flux_err)**2)
#     ln_like = (np.sum(1/np.sqrt(2*np.pi*(lc.flux_err))) - 0.5*chisq)
    
#     return ln_prior + ln_like

def lnprob(theta,
           lc,
           planet_period,
           LD,
           airmass,
           planet_a,
           expected_t0,
           plot = False):
        
    # Pull out some model parameters
    C1, C2, rp, t_0, a, i = theta
        
    # First we want a model to perform the lnprob calculation with.
    model = (C1 + C2*(airmass-1))*BATMAN(planet_period, rp, t_0, i, a, LD, t = lc.time)
    
    if plot:
        
        hires_times = np.linspace(lc.time[0],lc.time[-1],1000)
        model_to_plot = BATMAN(planet_period, rp, t_0, i, a, LD, t = hires_times)
        
        plt.figure()
        plt.errorbar(lc.time,lc.flux/(C1 + C2*(airmass-1)),yerr=lc.flux_err/(C1 + C2*(airmass-1)),
                    fmt='o',alpha=0.5,color='royalblue',markersize='5')
        plt.plot(hires_times,model_to_plot,label='Best-Fit Model',color='k',zorder=100)
        plt.ylabel('Normalized Flux',fontsize=18)
        plt.xlabel('BJD',fontsize=16)
        
    # We need to make sure the uniform priors are accounted for:
    ln_prior = lnprior(lc, theta, planet_a, expected_t0)
    
    # This is a Gaussian likelihood, for independent data points
    
    chisq = np.sum((lc.flux - model)**2/(lc.flux_err)**2)
    ln_like = (np.sum(1/np.sqrt(2*np.pi*(lc.flux_err))) - 0.5*chisq)
    
    return ln_prior + ln_like
    
def plot_chain(sampler,
               static_params,
               start=0,
               stop=0):
    
    planet_radius, planet_period, planet_a, planet_i, expected_t0, LD, offset, observatory, airmass, toi = static_params

    C1_all, C2_all, Rp_all, t0_all, A_all, I_all = sampler.chain.T
    
    ndim = np.shape(sampler.chain.T[:,0,0])[0]
    
    C1, C2, Rp, t0, A, Inc = sampler.chain[:, start:, :].reshape((-1, ndim)).T
    
    C1_burn, C2_burn, Rp_burn, t0_burn, A_burn, Inc_burn = sampler.chain[:, :start, :].reshape((-1, ndim)).T
    
    plt.figure(figsize=(14,14))
    gs = plt.matplotlib.gridspec.GridSpec(4,2)
    
    #Walker Plots
    
    ax_Rp = plt.subplot(gs[0,0])
    ax_Rp.plot(Rp_all.flatten(),color='black',alpha=0.5)
    ax_Rp.axvspan(start, stop, zorder=-1,alpha=0.3)
    ax_Rp.set_ylabel('Rp')
    
    ax_t0 = plt.subplot(gs[0,1])
    ax_t0.plot(t0_all.flatten(),color='black',alpha=0.5)
    ax_t0.axvspan(start, stop, zorder=-1,alpha=0.3)
    ax_t0.set_ylabel('t0')
    
    ax_A = plt.subplot(gs[1,0])
    ax_A.plot(A_all.flatten(),color='black',alpha=0.5)
    ax_A.axvspan(start, stop, zorder=-1,alpha=0.3)
    ax_A.set_ylabel('a/R*')
    
    ax_i = plt.subplot(gs[1,1])
    ax_i.plot(I_all.flatten(),color='black',alpha=0.5)
    ax_i.axvspan(start, stop, zorder=-1,alpha=0.3)
    ax_i.set_ylabel('Inclination')
    
    #Histograms
    
    ax_Rphist = plt.subplot(gs[2,0])
    ax_Rphist.hist(Rp.flatten(),color='black',bins=100)
    ax_Rphist.axvline(planet_radius,zorder=100,color='purple',label='Expected Value')
    ax_Rphist.set_xlabel('Rp')
    ax_Rphist.legend()
    
    ax_t0hist = plt.subplot(gs[2,1])
    ax_t0hist.hist(t0.flatten(),color='black',bins=100)
    ax_t0hist.axvline(expected_t0,zorder=100,color='purple',label='Expected Value')
    ax_t0hist.set_xlabel('t0')
    ax_t0hist.legend()
    
    ax_Ahist = plt.subplot(gs[3,0])
    ax_Ahist.hist(A.flatten(),color='black',bins=100)
    ax_Ahist.axvline(planet_a,zorder=100,color='purple',label='Expected Value')
    ax_Ahist.set_xlabel('a/R*')
    ax_Ahist.legend()
    
    ax_ihist = plt.subplot(gs[3,1])
    ax_ihist.hist(Inc.flatten(),color='black',bins=100)
    ax_ihist.axvline(planet_i,zorder=100,color='purple',label='Expected Value')
    ax_ihist.set_xlabel('Inclination')
    ax_ihist.legend()
    
    plt.savefig(f'Figs/TOI_{toi}_plotchain.pdf')
    
def corner(samples,
           labels,
           toi = 'no_toi_given'):
        
    Samples = samples.T
    
    np.random.seed(41)
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
    hires_times = np.linspace(lc.time[0]-1,lc.time[-1]+1,1000) #This array is for the model that gets plotted
    C1_best, C2_best, rp_best, t_0_best, a_best, i_best = parameters
    ndim = np.shape(sampler.chain.T[:,0,0])[0]
    burnin = int(0.25*nsteps)
    C1, C2, Rp, t0, A, Inc = sampler.chain[:, burnin:, :].reshape((-1, ndim)).T
    
    ########################################################################
    # Calculate the model with the best-fit MCMC params, calculate residuals
    ########################################################################
    
    model_to_plot = BATMAN(planet_period, rp_best, t_0_best, i_best, a_best, LD, t = hires_times)
    best_model = (C1_best + C2_best*(airmass-1))*BATMAN(planet_period, rp_best, t_0_best,
                                                        i_best, a_best, LD, t = lc.time)
    residual = (lc.flux-best_model)/lc.flux_err
    
    #########################
    # Plot the best-fit model
    #########################
            
    f, (a0, a1) = plt.subplots(2,1, gridspec_kw = {'height_ratios':[4,1]},
                               figsize=(14,7),sharex=True)
    
    a0.set_title('Modeled LCO Light Curve for TOI {}'.format(toi),fontsize=20)
    a0.errorbar((lc.time-t_0_best)*24,lc.flux/(C1_best + C2_best*(airmass-1)),
                yerr=lc.flux_err/(C1_best + C2_best*(airmass-1)),
                fmt='o',alpha=0.5,color='royalblue',markersize='5',
                label=observatory+' data')
    a0.plot((hires_times-t_0_best)*24,model_to_plot,label='Best-Fit Model',color='k',zorder=100)
    a0.set_ylabel('Normalized Flux',fontsize=18)
    a0.set_xlim(-3,3)
    a0.minorticks_on()
    a0.legend(loc='lower right')
    a0.tick_params(labelbottom=False, labeltop=False, labelleft=True, labelright=False,
                   bottom=True, top=True, left=True, right=True)
    
    #########################################
    # Plot 200 random models from the samples
    #########################################
    
    for j in range(0,200,1):
        i = np.random.randint(low=0,high=(nsteps-burnin)*100)
        sigma_model = BATMAN(P = planet_period, Rp = Rp[i], t0 = t0[i],
                             inc = Inc[i], A = A[i], LD = LD, t = hires_times)
        a0.plot((hires_times-t_0_best)*24,sigma_model,color='red',alpha = 0.1,
                linewidth=0.8,zorder=-1000,label='Random Samples')
        
    ####################
    # Plot the residuals
    ####################

    a1.scatter((lc.time-t_0_best)*24,residual,color='royalblue',alpha=0.5)
    a1.axhline(0,color='k')
    a1.set_ylim(0-1.5*np.max(np.abs(residual)),0+1.5*np.max(np.abs(residual)))
    a1.minorticks_on()
    a1.set_ylabel(r'Residuals ($\rm{\sigma}$)',fontsize=15)
    a1.set_xlabel('Time from Mid-Transit (hr)')
    a1.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False,
                   bottom=True, top=True, left=True, right=True)
        
    plt.xlabel('BJD - '+str(offset),fontsize=16)
    plt.savefig(f'Figs/TOI_{toi}_modeled.pdf')
    
    ###############################################
    # Download the TESS light curve for this object
    ###############################################
    
#     lcf = lk.search_lightcurvefile('TIC {}'.format(int(tic_id))).download_all()
#     Tlc = lcf[0].PDCSAP_FLUX.remove_outliers(sigma=10).normalize()
#     folded_Tlc = Tlc.fold(period=planet_period,t0=epoch).bin(binsize=10)
    
    ###################################################
    # Plot LCO and TESS data with the same model params
    ###################################################
    
            
#     f, (a0, a1) = plt.subplots(1,2,figsize=(14,5),sharey=True)
    
#     a0.errorbar((Tlc.time-epoch)*24,Tlc.flux,yerr=Tlc.flux_err,
#                 fmt='o',alpha=0.5,color='gray',markersize='5',
#                 label='TESS Data')
#     a0.errorbar((folded_Tlc.time*planet_period*24),folded_Tlc.flux,
#                 yerr=folded_Tlc.flux_err,
#                 fmt='o',alpha=0.8,color='royalblue',markersize='6',
#                 label='Folded and Binned TESS Data')
#     a0.plot((hires_times-t_0_best)*24,model_to_plot,label='Best-Fit Model',
#             color='k',zorder=100)

#     a0.set_title('TESS Light Curve for TOI {}'.format(toi),fontsize=20)
#     a0.set_xlabel('Time from Mid-Transit (hr)',fontsize=18)
#     a0.set_ylabel('Normalized Flux',fontsize=18)
#     a0.set_xlim(-1.5,1.5)
#     a0.set_ylim(np.median(lc.flux/(C1_best + C2_best*(airmass-1)))-8*np.median(lc.flux_err/(C1_best + C2_best*(airmass-1))),
#                 np.median(lc.flux/(C1_best + C2_best*(airmass-1)))+6*np.median(lc.flux_err/(C1_best + C2_best*(airmass-1))))
#     a0.minorticks_on()
#     a0.legend(loc='upper right')
#     a0.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False,
#                    bottom=True, top=True, left=True, right=True)
    
    
#     a1.errorbar((lc.time-t_0_best)*24,lc.flux/(C1_best + C2_best*(airmass-1)),
#                 yerr=lc.flux_err/(C1_best + C2_best*(airmass-1)),
#                 fmt='o',alpha=0.5,color='royalblue',markersize='5',
#                 label=observatory+' data')
#     a1.plot((hires_times-t_0_best)*24,model_to_plot,label='Best-Fit Model',color='k',zorder=100)

#     a1.set_title('Modeled LCO Light Curve for TOI {}'.format(toi),fontsize=20)
#     a1.set_xlabel('Time from Mid-Transit (hr)',fontsize=18)
#     a1.set_xlim(-1.5,1.5)
#     a1.minorticks_on()
#     a1.legend(loc='upper right')
#     a1.tick_params(labelbottom=True, labeltop=False, labelleft=False, labelright=False,
#                    bottom=True, top=True, left=True, right=True)
    
#     for j in range(0,200,1):
#         i = np.random.randint(low=0,high=(nsteps-burnin)*100)
#         sigma_model = BATMAN(P = planet_period, Rp = Rp[i], t0 = t0[i],
#                              inc = Inc[i], A = A[i], LD = LD, t = hires_times)
#         a1.plot((hires_times-t_0_best)*24,sigma_model,color='red',alpha = 0.1,
#                 linewidth=0.8,zorder=-1000,label='Random Samples')
    
#     plt.tight_layout()
#     plt.savefig('Figs/TOI_{}_LCOandTESS.pdf'.format(toi))
    
def rms_plot(lc,
             parameters,
             static_params):
    
    C1_best, C2_best, rp_best, t_0_best, a_best, i_best = parameters
    planet_radius, planet_period, planet_a, planet_i, expected_t0, LD, offset, observatory, airmass, toi = static_params
        
    lc_fixed = lc/(C1_best + C2_best*(airmass-1))
    
    best_model = BATMAN(planet_period, rp_best, t_0_best, i_best, a_best, LD, t = lc_fixed.time)
    
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