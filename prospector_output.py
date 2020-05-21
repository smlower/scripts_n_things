# coding: utf-8
from corner import quantile
import prospect.io.read_results as pread
from prospect.models.transforms import zfrac_to_sfrac, logsfr_ratios_to_sfrs, zfrac_to_sfr, tburst_from_fage, zfrac_to_masses, logsfr_ratios_to_masses
import numpy as np
import sys
import glob, os
from prospector_output_utilities import *
import pandas as pd
from time import time





t = time()
print('loaded stuff')
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

print('read args')


def get_best(res, **kwargs):
    imax = np.argmax(res['lnprobability'])
    theta_best = res['chain'][imax, :].copy()

    return theta_best



#appending path of prosp run script to import model and sps
sys.path.append(prosp_dir)
from run_prosp import build_model, build_sps



smass_map = []
smass_50 = []
smass_16 = []
smass_84 = []
mass_frac_list = []
int_smass = []
sfr_map = []
sfr_50 = [] 
sfr_unc16 = []
sfr_unc84 = []
Z_med = []
gal_list = []
mweighted_age = []
time_list = []
for galaxy_number in [galaxies]:
    print('in loop')
    res, obs, mod = 0, 0, 0
    galaxy_num = "{:03d}".format(galaxy_number)
    infile = prosp_dir+'/galaxy'+galaxy_num+'.h5'
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
    sps = build_sps()
    mod = build_model()
    thetas = mod.theta_labels()
    theta_best = get_best(res)
    mean_mod = mod.mean_model(theta_best, obs, sps)
    mass_frac = mean_mod[-1]
    mass_frac_list.append(mass_frac)
    if any('massmet_1' in s for s in thetas):
        mass_index = [i for i, s in enumerate(thetas) if 'massmet_1' in s]
        mass_map = (10**theta_best[mass_index[0]]) * mass_frac
        mass_chain = [item[mass_index] for item in res['chain']]
        mass_quantiles = quantile(mass_chain, [.16, .5, .84])
        mass = mass_quantiles[1] #median mass
        smass_50.append(mass_quantiles[1])
        smass_16.append(mass_quantiles[0])
        smass_84.append(mass_quantiles[2])
        smass_map.append(mass_map)
        gal_list.append(galaxy_num)
        metal_idx = [i for i, s in enumerate(thetas) if 'massmet_2' in s]
        Z_chain = [item[metal_idx] for item in res['chain']]
        Z_quantiles = quantile(Z_chain, [.16, .5, .84])
        Z_med.append(Z_quantiles[1])
        
    else:
        mass_index = [i for i, s in enumerate(thetas) if 'mass' in s]
        metal_idx = [i for i, s in enumerate(thetas) if 'logzsol' in s]
        if thetas[mass_index[0]] == 'mass' or thetas[mass_index[0]] == 'total_mass':
            mass_map = theta_best[mass_index[0]] * mass_frac
            mass_chain = [item[mass_index] for item in res['chain']]
            mass_quantiles = quantile(mass_chain, [.16, .5, .84])
            mass = mass_quantile[1]
            smass_50.append(mass_quantiles[1])
            smass_16.append(mass_quantiles[0])
            smass_84.append(mass_quantiles[2])
            smass_map.append(mass_map)
            
        else:
            mass_map = (10**theta_best[mass_index[0]]) * mass_frac
            mass_chain = [item[mass_index] for item in res['chain']]
            mass_quantiles = quantile(mass_chain, [.16, .5, .84])
            mass = mass_quantile[1]
            smass_50.append(mass_quantiles[1])
            smass_16.append(mass_quantiles[0])
            smass_84.append(mass_quantiles[2])
            smass_map.append(mass_map)
        
        gal_list.append(galaxy_num)
        Z_chain = [item[metal_idx] for item in res['chain']]
        Z_quantiles = quantile(Z_chain, [.16, .5, .84])
        Z_med.append(Z_quantiles[1])

#now break down routines for different SFH models
   
    if any('tau' in s for s in thetas):
        tau_idx = [i for i, s in enumerate(thetas) if 'tau' in s]
        tau_chain = [item[tau_index] for item in res['chain']]
        tau_quan= quantile(tau_chain, [.16, .5, .84])
        tau50 = tau_quan[1]
        tau16 = tau_quan[0]
        tau84 = tau_quan[2]
        stellar_mass = 10**(mass)  #normalization is total mass formed
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
            norm_burst = (mfrac_burst[0] * stellar_mass) / (1e8)#((time[idx + 1] - time[idx - 1])*1.0e9) 
            sfr_onegal = np.insert(sfr_onegal, idx, norm_burst)
            sfh_type = 'ParametricTau+Burst'
            sfr_50.append(np.ravel(sfr_onegal))
            time_list.append(time)
            age_mw = tau_massweighted_age(np.ravel(sfr_onegal), time)
            mweighted_age.append(age_mw)
            
        else:
            norm = stellar_mass / (age*1.0e9)
            sfr_onegal = []
            for i in time:
                sfr_now = norm * (i-age) * np.exp(-(i-age) / tau)
                sfr_onegal.append(sfr_now)
            sfr_50.append(sfr_onegal)
            sfh_type = 'ParametricTau'
            print('SFH:', np.ravel(sfr_onegal))
            age_mw = tau_massweighted_age(np.ravel(sfr_onegal), time)
            mweighted_age.append(age_mw)
            print('mwage:', age_mw)


    elif any('tage' in s for s in thetas):
        print('got to SFH')
        stellar_mass = mass / mass_frac
        age_idx = thetas.index('tage')
        age = theta_best[age_idx] #Gyr                                                                                 
        norm = stellar_mass / (age*1.0e9)
        sfr_50.append(norm)
        sfh_type = 'Constant'
        time = np.arange(0, 13.8, 0.05)
        age_mw = const_massweighted_age(np.full(shape=len(time), fill_value=norm), time)
        mweighted_age.append(age_mw)
        
            
    elif any('z_fraction' in s for s in thetas):
        total_mass = (10**mass)
        zfrac_idx = [i for i, s in enumerate(thetas) if 'z_fraction' in s]
        zfrac_bins = []
        for i in zfrac_idx:
            zfrac_bins.append(theta_best[i])     
        time_bins_log = next(item for item in res['model_params'] if item["name"] == "agebins")['init']
        time = []
        for val in time_bins_log:
            time.append(13.8 - (10**val / 1.0e9)) #Gyr
        #sfr_onegal = zfrac_to_sfr(total_mass, zfrac_bins, time_bins_log)
        #sfr.append(sfr_onegal)

        #get 16th and 84th quantiles for SFR bins
        zfrac_chain = [item[zfrac_idx[0]:zfrac_idx[-1]+1]  for item in res['chain']]
        sfr_chain = []
        mass_chain = []
        for i in range(len(zfrac_chain)):
            sfr_chain.append(zfrac_to_sfr(total_mass, zfrac_chain[i], time_bins_log))
            mass_chain.append(zfrac_to_masses(total_mass, zfrac_chain[i], time_bins_log))
        sfr_unc = []
        mass_unc = []
        for i in range(np.shape(zfrac_chain)[1]+1): 
            sfr_unc.append(quantile([item[i] for item in sfr_chain], [.16, .5, .84]))
            mass_unc.append(quantile([item[i] for item in mass_chain], [.16, .5, .84]))
        sfr_16 = [item[0] for item in sfr_unc]
        sfr_med = [item[1] for item in sfr_unc]
        sfr_84 = [item[2] for item in sfr_unc]
        mass_50 = [item[1] for item in mass_unc]
        sfr_50.append(sfr_med)
        sfr_unc16.append(sfr_16)
        sfr_unc84.append(sfr_84)
        sfh_type = 'Dirichlet'


        age = nonpara_massweighted_age(mass_50, np.unique(np.ravel(time_bins_log)), total_mass)
        mweighted_age.append(age)

    elif any('logsfr_ratios' in s for s in thetas):
        #mass_idx = [i for i, s in enumerate(thetas) if 'massmet_1' in s]
        logmass = np.log10(mass)
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
        mass_chain = []
        for i in range(len(sfrfrac_chain)):
            sfr_chain.append(logsfr_ratios_to_sfrs(logmass, sfrfrac_chain[i], time_bins_log))
            mass_chain.append(logsfr_ratios_to_masses(logmass, sfrfrac_chain[i], time_bins_log))
        sfr_unc = []
        mass_unc = []
        for i in range(np.shape(sfrfrac_chain)[1] + 1):
            sfr_unc.append(quantile([item[i] for item in sfr_chain], [.16, .5, .84]))
            mass_unc.append(quantile([item[i] for item in mass_chain], [.16, .5, .84]))
        sfr_16 = [item[0] for item in sfr_unc]
        sfr_50 = [item[1] for item in sfr_unc]
        sfr_84 = [item[2] for item in sfr_unc]
        mass_50 = [item[1] for item in mass_unc]
        sfr.append(sfr_50)
        sfr_unc16.append(sfr_16)
        sfr_unc84.append(sfr_84)
        sfh_type = 'Continuity' 


        age = nonpara_massweighted_age(mass_50, np.unique(np.ravel(time_bins_log)), 10**logmass)
        mweighted_age.append(age)


if 'Constant' in sfh_type:
    sfr_list = np.ravel(sfr)
    galaxy_list = [item for item in gal_list]
    dmass_list = [item for item in der_smass]
    imass_list = [item for item in int_smass ]
    massfrac_list = [item for item in mass_frac_list]
    metal_list = [item for item in metallicity ]
    age_list = [item for item in mweighted_age ]
    data = {'Galaxy': galaxy_list, 'SFR' : sfr_list, 'Derived Stellar Mass' : dmass_list, 'Intrinsic Stellar Mass' : imass_list, 'Mass Fraction': massfrac_list, 'Metallicity': metal_list, 'Mass-weighted Age': age_list}



if 'Parametric' in sfh_type:
    sfr_list = np.ravel(sfr_50)
    timelist = np.ravel(time_list)
    galaxy_list = [item for item in gal_list for repetitions in range(len(sfr_onegal))]
    dmass_list = [item for item in smass_50 for repetitions in range(len(sfr_onegal))]
    imass_list = [item for item in int_smass for repetitions in range(len(sfr_onegal))]
    mass16_list = [item for item in smass_16 for repetitions in range(len(time))]
    mass84_list = [item for item in smass_84 for repetitions in range(len(time))]
    massfrac_list = [item for item in mass_frac_list for repetitions in range(len(sfr_onegal))]
    metal_list = [item for item in Z_med for repititions in range(len(sfr_onegal))]
    age_list = [item for item in mweighted_age for repititions in range(len(sfr_onegal))]
    data = {'Galaxy': galaxy_list, 'Time': timelist, 'SFR' : sfr_list,'Mass_16' : mass16_list, 'Mass_84': mass84_list, 'Derived Stellar Mass' : dmass_list, 'Intrinsic Stellar Mass' : imass_list, 'Mass Fraction': massfrac_list, 'Metallicity': metal_list, 'Mass-weighted Age': age_list}
else: 
    map_mass_list = [item for item in smass_map for repetitions in range(2*len(time))]
    sfr_list = [item for item in np.ravel(sfr_50) for repititions in range(2)]
    dmass_list = [item for item in smass_50 for repetitions in range(2*len(time))]
    mass16_list = [item for item in smass_16 for repetitions in range(2*len(time))]
    mass84_list = [item for item in smass_84 for repetitions in range(2*len(time))]
    imass_list = [item for item in int_smass for repetitions in range(2*len(time))]
    massfrac_list = [item for item in mass_frac_list for repetitions in range(2*len(time))]
    metal_list = [item for item in Z_med for repititions in range(2*len(time))]
    sfr_16list = [item for item in np.ravel(sfr_unc16) for repetitions in range(2)]
    sfr_84list = [item for item in np.ravel(sfr_unc84) for repetitions in range(2)]
    age_list = [item for item in mweighted_age for repititions in range(2*len(time))]
    galaxy_list = [item for item in gal_list for repetitions in range(2*len(time))]
    timelist = np.ravel(time)
    data = {'Galaxy': galaxy_list, 'Time': timelist, 'SFR' : sfr_list, 'SFR_16' : sfr_16list, 'SFR_84' : sfr_84list, 'MAP Mass' : map_mass_list, 'Derived Stellar Mass' : dmass_list, 'Mass_16' : mass16_list, 'Mass_84': mass84_list, 'Intrinsic Stellar Mass' : imass_list, 'Mass Fraction' : massfrac_list, 'Metallicity' : metal_list, 'Mass-weighted Age': age_list}

#labels = [gal_list, np.ravel(time)]

#index = pd.MultiIndex.from_product(labels, names=['Galaxy', 'Time [Gyr]'])

#s = pd.DataFrame(data, index=index)

#s.to_pickle(prosp_dir+'/prospectorsnap305_properties_'+sfh_type+'_'+str(t)+'.pkl')


np.savez(prosp_dir+'/galaxy_'+str(galaxies)+'.npz', data=data)



