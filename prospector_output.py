# coding: utf-8
import matplotlib 
matplotlib.use('Agg')
from corner import quantile
import prospect.io.read_results as pread
from prospect.models.transforms import zfrac_to_sfrac, logsfr_ratios_to_sfrs, zfrac_to_sfr, tburst_from_fage
import numpy as np
import sys
import glob, os
from prospector_output_utilities import *
import pandas as pd
from time import time

t = time()

def arg_parser(args):
    for arg in args:
        if arg.startswith('--prosp_dir='):
            prosp_dir = arg[12:]

        elif arg.startswith('--galaxies='):
            galaxies = int(arg[11:])
        
        elif arg.startswith('--pd_dir='):
            pd_dir = arg[9:]
            
    return prosp_dir, pd_dir, galaxies

prosp_dir, pd_dir, galaxies = arg_parser(args=sys.argv[1:])

der_smass = []
int_smass = []
sfr = [] 
sfr_unc16 = []
sfr_unc84 = []
metallicity=[]
gal_list = []
mweighted_age = []
for galaxy_number in range(galaxies):
    res, obs, mod = 0, 0, 0
    galaxy_num = "{:03d}".format(galaxy_number)
    infile = prosp_dir+'/snap*.galaxy'+galaxy_num+'_*_mcmc.h5'
    globfiles = glob.glob(infile)
    if len(globfiles) == 0:
        print('No files matching glob pattern %s' % infile)
        continue
    try:
        for prosp_output in glob.glob(infile):
            print(prosp_output)
            res, obs, mod = pread.results_from(prosp_output)
        pdfile = pd_dir+'/grid_physical_properties.305_galaxy'+galaxy_num+'.npz'
        pd_data = np.load(pdfile)
        int_smass.append(np.sum(pd_data['grid_star_mass']))
    except:
        print('File corrupted')
        continue
    sps = pread.get_sps(res)
    thetas, theta_best = get_best(res)
    mean_mod = mod.mean_model(theta_best, obs, sps)
    mass_frac = mean_mod[-1]
    if any('massmet_1' in s for s in thetas):
        mass_index = [i for i, s in enumerate(thetas) if 'massmet_1' in s]
        der_smass.append(10**theta_best[mass_index[0]] * mass_frac)
        gal_list.append(galaxy_num)
        metal_idx = [i for i, s in enumerate(thetas) if 'massmet_2' in s]
        metallicity.append(theta_best[metal_idx[0]])
    else:
        mass_index = [i for i, s in enumerate(thetas) if 'mass' in s]
        metal_idx = [i for i, s in enumerate(thetas) if 'logzsol' in s]
        if thetas[mass_index[0]] == 'mass' or thetas[mass_index[0]] == 'total_mass':
            der_smass.append(theta_best[mass_index[0]] * mass_frac)
        else: der_smass.append(10**theta_best[mass_index[0]] * mass_frac)
        gal_list.append(galaxy_num)
        metallicity.append(theta_best[metal_idx[0]]) 


   
    if any('tau' in s for s in thetas):
        tau_idx = [i for i, s in enumerate(thetas) if 'tau' in s]
        tau = theta_best[tau_idx]
        mass_idx = [i for i, s in enumerate(thetas) if 'mass' in s]
        stellar_mass = theta_best[mass_idx]
        #fburst also has an fburst_age which screws with the normal recipe, fix below
        age_idx = thetas.index('tage')
        age = theta_best[age_idx] #Gyr
        time = np.arange(0, 13.8, 0.05)
        if 'fburst' in thetas:
            mfrac_burst_idx = [i for i, s in enumerate(thetas) if 'fburst' in s]
            mfrac_burst = theta_best[mfrac_burst_idx]
            mfrac_else = 1.0 - mfrac_burst
            fage_burst_idx = [i for i, s in enumerate(thetas) if 'fage_burst' in s]
            fage_burst = theta_best[fage_burst_idx]
            tburst = tburst_from_fage(age, fage_burst)
            norm_tau = (mfrac_else * stellar_mass) / (age * 1.0e9)
            sfr_onegal = []

            for i in time:
                sfr_tau = norm_tau * (i - age) * np.exp(-(i-age) / tau)
                sfr_onegal.append(sfr_tau)
            time = np.append(time, tburst[0])
            time.sort()
            idx = np.where(time == tburst[0])[0][0]
            norm_burst = (mfrac_burst[0] * stellar_mass) / ((time[idx + 1] - time[idx - 1])*1.0e9) 
            sfr_onegal = np.insert(sfr_onegal, idx, norm_burst)
            sfh_type = 'ParametricTau+Burst'
            sfr.append(sfr_onegal)

            age_mw = tau_massweighted_age(sfr_onegal, time)
            mweighted_age.append(age_mw)
            
        else:
            norm = stellar_mass / (age*1.0e9)
            sfr_onegal = []
            for i in time:
                sfr_now = norm * (i-age) * np.exp(-(i-age) / tau)
                sfr_onegal.append(sfr_now)
            sfr.append(sfr_onegal)
            sfh_type = 'ParametricTau'
                  
            age_mw = tau_massweighted_age(sfr_onegal, time)
            mweighted_age.append(age_mw)

  
            
    if any('z_fraction' in s for s in thetas):
        #mass_idx = [i for i, s in enumerate(thetas) if 'total_mass' in s]
        total_mass = (10**theta_best[mass_index]) * mass_frac
        zfrac_idx = [i for i, s in enumerate(thetas) if 'z_fraction' in s]
        #zfrac_bins = []
        #for i in zfrac_idx:
        #    zfrac_bins.append(theta_best[i])     
        time_bins_log = next(item for item in res['model_params'] if item["name"] == "agebins")['init']
        time = []
        for val in time_bins_log:
            time.append(13.6 - (10**val / 1.0e9)) #Gyr
        #sfr_onegal = zfrac_to_sfr(total_mass, zfrac_bins, time_bins_log)
        #sfr.append(sfr_onegal)

        #get 16th and 84th quantiles for SFR bins
        zfrac_chain = [item[zfrac_idx[0]:zfrac_idx[-1]+1]  for item in res['chain']]
        sfr_chain = []
        for i in range(len(zfrac_chain)):
            sfr_chain.append(zfrac_to_sfr(total_mass, zfrac_chain[i], time_bins_log))

        sfr_unc = []
        for i in range(np.shape(zfrac_chain)[1]+1): 
            sfr_unc.append(quantile([item[i] for item in sfr_chain], [.16, .5, .84]))
        sfr_16 = [item[0] for item in sfr_unc]
        sfr_50 = [item[1] for item in sfr_unc]
        sfr_84 = [item[2] for item in sfr_unc]
        sfr.append(sfr_50)
        sfr_unc16.append(sfr_16)
        sfr_unc84.append(sfr_84)
        sfh_type = 'Dirichlet'


        age = nonpara_massweighted_age(sfr_50, np.unique(np.ravel(time)))
        mweighted_age.append(age)

    if any('logsfr_ratios' in s for s in thetas):
        #mass_idx = [i for i, s in enumerate(thetas) if 'massmet_1' in s]
        logmass = np.log10((10**theta_best[mass_index]) * mass_frac)
        sfrfrac_idx = [i for i, s in enumerate(thetas) if 'sfr' in s]
        #sfrfrac_bins = []
        #for i in sfrfrac_idx:
        #    sfrfrac_bins.append(theta_best[i])
        time_bins_log = next(item for item in res['model_params'] if item["name"] == "agebins")['init']
        time = []
        for val in time_bins_log:
            time.append(10**val / 1.0e9) #Gyr
        #sfr_onegal = logsfr_ratios_to_sfrs(logmass, sfrfrac_bins, time_bins_log)
        #sfr.append(sfr_onegal)

        #get 16th and 84th quantiles for SFR bins
        sfrfrac_chain = [item[sfrfrac_idx[0]:sfrfrac_idx[-1]+1] for item in res['chain']]
        sfr_chain = []
        for i in range(len(sfrfrac_chain)):
            sfr_chain.append(logsfr_ratios_to_sfrs(logmass, sfrfrac_chain[i], time_bins_log))
        sfr_unc = []
        for i in range(np.shape(sfrfrac_chain)[1] + 1):
            sfr_unc.append(quantile([item[i] for item in sfr_chain], [.16, .5, .84]))
        sfr_16 = [item[0] for item in sfr_unc]
        sfr_50 = [item[1] for item in sfr_unc]
        sfr_84 = [item[2] for item in sfr_unc]
        sfr.append(sfr_50)
        sfr_unc16.append(sfr_16)
        sfr_unc84.append(sfr_84)
        sfh_type = 'Continuity' 


        age = nonpara_massweighted_age(sfr_50, np.unique(np.ravel(time)))
        mweighted_age.append(age)
   

if 'Parametric' in sfh_type:
    sfr_list = np.ravel(sfr)
    dmass_list = [item for item in der_smass for repetitions in range(len(time))]
    imass_list = [item for item in int_smass for repetitions in range(len(time))]
    metal_list = [item for item in metallicity for repititions in range(len(time))]
    age_list = [item for item in mweighted_age for repititions in range(len(time))]
    data = {'SFR' : sfr_list, 'Derived Stellar Mass' : dmass_list, 'Intrinsic Stellar Mass' : imass_list, 'Metallicity': metal_list, 'Mass-weighted Age': age_list}
else: 
    sfr_list = [item for item in np.ravel(sfr) for repititions in range(2)]
    dmass_list = [item for item in der_smass for repetitions in range(2*len(time))]
    imass_list = [item for item in int_smass for repetitions in range(2*len(time))]
    metal_list = [item for item in metallicity for repititions in range(2*len(time))]
    sfr_16list = [item for item in np.ravel(sfr_unc16) for repetitions in range(2)]
    sfr_84list = [item for item in np.ravel(sfr_unc84) for repetitions in range(2)]
    age_list = [item for item in mweighted_age for repititions in range(2*len(time))]
    data = {'SFR' : sfr_list, 'SFR_16' : sfr_16list, 'SFR_84' : sfr_84list, 'Derived Stellar Mass' : dmass_list, 'Intrinsic Stellar Mass' : imass_list, 'Metallicity' : metal_list, 'Mass-weighted Age': age_list}

labels = [gal_list, np.ravel(time)]

index = pd.MultiIndex.from_product(labels, names=['Galaxy', 'Time [Gyr]'])

s = pd.DataFrame(data, index=index)

s.to_pickle(prosp_dir+'/prospectorsnap305_properties_'+sfh_type+'_'+str(t)+'.pkl')




