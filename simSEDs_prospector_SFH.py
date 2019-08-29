
# coding: utf-8

# In[3]:


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import prospect.io.read_results as pread
import numpy as np
from prospect.transforms import zfrac_to_sfrac, logsfr_ratios_to_sfrs
import sys
import glob, os
import pickle
#import caesar
import track_progens

# In[2]:


"""
Sys args:

1. Intrinsic files (via caesar)
2. Derived files (.h5 prospector output)
"""


# In[53]:


def get_best(res, **kwargs):
    """Get the maximum a posteriori parameters.
    From prospect.utils.plotting
    """
    imax = np.argmax(res['lnprobability'])
    # there must be a more elegant way to deal with differnt shapes
    try:
        i, j = np.unravel_index(imax, res['lnprobability'].shape)
        theta_best = res['chain'][i, j, :].copy()
    except(ValueError):
        theta_best = res['chain'][imax, :].copy()

    theta_names = res.get('theta_labels', res['theta_labels'])
    return theta_names, theta_best



#command line args


def arg_parse(args):
    
    for arg in args:
        if arg.startswith('--snap_dir='):
            snap_dir = arg[11:]
        if arg.startswith('--nonpar_dir='):
            nonpar_dir = arg[12:]
        if arg.startswith('--par_dir='):
            par_dir = arg[10:]
        if arg.startswith('--sfh_type='):
            sfh_type = arg[12:]
        if arg.startswith('--galaxy_ids='):
            galaxy_ids = arg[14:]

    return snap_dir, nonpar_dir, par_dir, sfh_type, galaxy_ids

snap_dir, nonpar_dir, par_dir, sfh_type, galaxy_ids = arg_parse(args=sys.argv[1:])
snap_range = range(54, 306)


for galaxy in galaxy_ids:


    #Get intrinsic SFH from caesar - track_progens written by Romeel Dave
    
    snaptime, intrinsic_sfr = track_progens.sfh_from_progens(snap_dir, galaxy_ids, snap_range, other_props=False)


    #Get derived from prospect output 
    nonpar_file = glob.glob(nonpar_dir+'snap305.galaxy%03d' % (galaxy)+'*.h5')[0]
    par_file = glob.glob(par_dir+'snap305.galaxy%03d' % (galaxy)+'*.h5')[0]
    res_nonpar, _,  _ = pread.results_from(nonpar_file)
    res_par, _, _ = pread.results_from(par_file)


    thetas_par, theta_best_par = get_best(res_par)
    #load parametric SFH
    tau_idx = [i for i, s in enumerate(thetas_par) if 'tau' in s]
    tau = theta_best_par[tau_idx]
    mass_idx = [i for i, s in enumerate(thetas_par) if 'mass' in s]
    stellar_mass = theta_best_par[mass_idx]
    age_idx = [i for i, s in enumerate(thetas_par) if 'age' in s]
    age = theta_best_par[age_idx] #Gyr                                                                                                                                                                                               
    time = np.arange(0, 13.5, 0.5)
    if 'fburst' in thetas_par:
        mfrac_burst_idx = [i for i, s in enumerate(thetas_par) if 'fburst' in s]
        mfrac_burst = theta_best_par[mfrac_burst_idx]
        mfrac_else = 1.0 - mfrac_burst
        fage_burst_idx = [i for i, s in enumerate(thetas_par) if 'fage_burst' in s]
        fage_burst = theta_best_par[fage_burst_idx]
        tburst = tburst_from_fage(age, fage_burst)
        norm_tau = (mfrac_else * stellar_mass) / age
        sfr_par = []
        for i in time:
            sfr_tau = norm_tau * np.exp(-i / tau)
            sfr_par.append(sfr_tau)
        time.append(tburst)
        time.sort()
        idx = test.index(tburst)
        sfr_par = np.sort(sfr_par, idx, mfrac_burst*stellar_mass)
        ifburst = '+burst'
      
    else:
        norm = stellar_mass / age
        sfr_par = []
        for i in time:
            sfr_par.append(norm * np.exp(-i / tau))
        ifburst = ''
    
    #Fancy footwork to get sfr from non parametric...parameters
    thetas, theta_best = get_best(res_nonpar)
    #gotta handle different parameter outputs for each model
    if model == 'continuity':
        mass_idx = [i for i, s in enumerate(thetas) if 'logmass' in s]
        logmass = theta_best[mass_idx]
        sfrfrac_idx = [i for i, s in enumerate(thetas) if 'sfr' in s]
        sfrfrac_bins = []
        for i in sfr_idx:
            sfrfrac_bins.append(theta_best[i])
        time_bins_log = next(item for item in res['model_params'] if item["name"] == "agebins")['init']
        time_bins = []
        for time in time_bins_log:
            time_bins.append(10**time / 1.0e9) #Gyr
        sfr_bins = logsfr_ratios_to_sfrs(logmass, sfrfrac_bins, time_bins)
       
        '''time_bins_unraveled = list(set(np.ravel(time_bins)))
        time_in_each_bin = [y-x for x, y in zip(time_bins_unraveled, time_bins_unraveled[1:])]
        mass_idx = [i for i, s in enumerate(thetas) if 'logmass' in s]
        logmass = theta_best[mass_idx]
        stellar_mass = 10**logmass
        frac_sfr_bins = []
        for i in np.arange(len(sfr_bins)):
            frac_sfr_bins.append((stellar_mass * (10**sfr_bins[i])) / (time_in_each_bin[i] * 1.0e9))'''

    if model == 'dirichlet':
        zfrac_idx = [i for i, s in enumerate(thetas) if 'zfrac' in s]
        zfrac_bins = []
        for i in zfrac_idx:
            zfrac_bins.append(theta_best[i])
        sfr_bins = zfrac_to_sfrac(zfrac_bins)
        time_bins_log = next(item for item in res['model_params'] if item["name"] == "agebins")['init']
        time_bins = []
        for time in time_bins_log:
            time_bins.append(10**time / 1.0e9) #Gyr
    
    derived_sfr = sfr_bins
    model_time = time_bins
    


    #plot
    plt.figure(figsize=(10, 8))
    plt.tick_params(axis='both', which='major', labelsize=13)
    plt.tick_params(axis='both', which='minor', labelsize=10)
    #parametric
    plt.plot(time, sfr_par, color='blue', lw=3, label='Derived SFH_par'+ifburst)
    #non parametric
    for i in range(len(frac_sfr_bins)):
        plt.plot([time_bins[i][0], time_bins[i][1]], [frac_sfr_bins[i], frac_sfr_bins[i]], color='orange', lw=3, label='Derived SFH_'+model)
    #intrinsic 
    plt.plot(snaptime, intrinsic_sfr, color='black', lw=0.5, label='Intrinsic SFH')
    plt.ylabel('SFR [M$_{\odot}$ / yr]', fontsize=25)
    plt.xlabel('Lookback Time [Gyr]', fontsize=25)
    plt.xlim([0, 14])
    plt.legend()
    plt.savefig(prosp_dir+'SFH_galaxy'+str(galaxy)+'.png', dpi=300)






