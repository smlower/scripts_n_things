# coding: utf-8
import prospect.io.read_results as pread
from prospect.models.transforms import zfrac_to_sfrac, logsfr_ratios_to_sfrs, zfrac_to_sfr, tburst_from_fage
import numpy as np
import sys
from glob import glob 
import os
import pandas as pd
from prospector_output_utilities import *



def get_best_mod(res, mod, **kwargs):
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

    theta_names = res.get('theta_labels', mod.theta_labels())
    return theta_names, theta_best



def build_model(attn_mod, free_dust, **kwargs):

    from prospect.models.templates import TemplateLibrary
    from prospect.models import priors, sedmodel


    model_params = TemplateLibrary["parametric_sfh"]
    model_params.update(TemplateLibrary["dust_emission"])

    #fixed delayed-tau SFH                                                                                 
    model_params['tau']['isfree'] = True
    model_params["lumdist"] = {"N": 1, "isfree": False,
                         "init": 1.0e-5, "units": "Mpc"}
    model_params['tage']['prior'] = priors.TopHat(mini=0.0, maxi=13.8)
    model_params['mass']['prior'] = priors.TopHat(mini=1e7, maxi=1e13)
    model_params['logzsol']['isfree'] = True
    model_params['logzsol']['prior'] = priors.ClippedNormal(mini=-1.8, maxi=0.2, mean=0.0, sigma=0.3)

    #dust emission                                                                                         
    
    model_params['duste_gamma']['isfree'] = free_dust
    model_params['duste_qpah']['isfree'] = free_dust
    model_params['duste_umin']['isfree'] = free_dust
    model_params = model_library(model_params, attn_mod, True)
    

    model = sedmodel.SedModel(model_params)


    return model


def model_library(model_params, attenuation_model, dust1=False):

    from prospect.models import priors
    if attenuation_model == 0:
        name = 'MW_fixed'
        model_params["dust_type"] = {"N": 1, "isfree": False, "init": 1, "units": name}
        model_params["mrw"] = {"N": 1, "isfree": False,
                           "init": 3.1, "units": "The ratio of total to selective absorption which characterizes the MW\
 extinction curve", 'prior' : priors.TopHat(mini=2.0, maxi=5.0)}
        model_params['uvb'] = {'N': 1, 'isfree': False,
                           'init': 1.0, 'units':"characterizing the strength of the 2175A extinction feature with respe\
ct to the standard Cardelli et al. determination for the MW",
                           'prior': priors.TopHat(mini = 0.1, maxi=3.0)}
        model_params["dust1"]  = {"N": 1, "isfree": dust1,
                     "init": 0.0, "prior" : priors.TopHat(mini=0.0, maxi=1.5),
                     "units": "optical depth towards young stars"}

    elif attenuation_model == 1:
        name = 'SMC'
        model_params["dust_type"] = {"N": 1, "isfree": False, "init": 0, "units": name} #power-law           
        model_params["dust2"]["prior"] = priors.TopHat(mini=0.0, maxi=3.0)
        model_params["dust1"]  = {"N": 1, "isfree": dust1,
                     "init": 0.0, "prior" : priors.TopHat(mini=0.0, maxi=1.5),
                     "units": "optical depth towards young stars"}
        model_params["dust_index"] = {"N": 1, "isfree": False,
                     "init": -1.19, "units": "power-law multiplication of Calzetti",
                     "prior": priors.TopHat(mini=-2.0, maxi=0.5)}
        model_params["dust_index1"] = {"N": 1, "isfree": False,
                     "init": -1.0, "units": "power-law multiplication of Calzetti",
                     "prior": priors.TopHat(mini=-2.0, maxi=0.5)}


    elif attenuation_model == 2:
        name = 'Calzetti_fixed'
        model_params["dust_type"] = {"N": 1, "isfree": False, "init": 2, "units": name} #Calzetti            
        model_params["dust2"]["prior"] = priors.TopHat(mini=0.0, maxi=3.0)
        model_params["dust1"] = {"N": 1, "isfree": False,
                         "init": 0.0, "units": "optical depth towards young stars"}

    return model_params



def build_sps(zcontinuous=1, compute_vega_mags=False, **extras):
    from prospect.sources import CSPSpecBasis
    sps = CSPSpecBasis(zcontinuous=zcontinuous,
                       compute_vega_mags=compute_vega_mags)
    return sps



free_dust = int(sys.argv[1])




bf_mass_smc = []
bf_mass_cal = []
sfr_smc = []
sfr_cal = []
true_mass = []
true_sfr = []
galaxy_id = [] 

smc_dir = '/ufrc/narayanan/s.lower/atten_test/simba_runs/SMC/'
cal_dir = '/ufrc/narayanan/s.lower/atten_test/simba_runs/calzetti/'
pd_dir = '/ufrc/narayanan/s.lower/pd_runs/simba_m25n512/snap305_boxtest/snap305'
caesar_dir = '/ufrc/narayanan/s.lower/caesar_sfr50.npz'

if free_dust == 0:
    smc_dir = smc_dir+'fixed_dust/'
    cal_dir = cal_dir+'fixed_dust/'
    #dust1 = 0


print('Building sps object')
sps = build_sps()

for galaxy in range(0, 500):
    galaxy_num = "{:03d}".format(galaxy)
    caesar = np.load(caesar_dir, allow_pickle=True)
    smc, cal, pdy = 0, 0, 0
    for infile in glob(smc_dir+'snap*_galaxy'+galaxy_num+'_*.h5'):
        smc, obs, _ = pread.results_from(infile)
    for infile in glob(cal_dir+'snap*_galaxy'+galaxy_num+'_*.h5'):
        cal, _, _ = pread.results_from(infile)
    for infile in glob(pd_dir+'/grid_physical_properties.305_galaxy'+galaxy_num+'.npz'):
        pdy = np.load(infile, allow_pickle=True)

    if isinstance(smc, int) == True or isinstance(cal, int) == True or\
 isinstance(pdy, int) == True:
        #print('smc: ',isinstance(smc, int))                           
        #print('cal: ',isinstance(cal, int))                           
        #print('pd: ',isinstance(pd, int))                             
        print('galaxy %s corrupted'%(galaxy))
        continue
    else:
        print('onto galaxy %s'%(galaxy))
        true_mass.append(np.sum(pdy['grid_star_mass'])/1.989e33)
        true_sfr.append(np.log10(caesar['sfr_50'][galaxy]))
        mod_smc = build_model(1, free_dust)
        mod_cal = build_model(2, free_dust)
        #sps = build_sps()
        smc_thetas, best_smc = get_best_mod(smc, mod_smc)
        cal_thetas, best_cal = get_best_mod(cal, mod_cal)
        mean_mod_smc = mod_smc.mean_model(best_smc, obs, sps)
        mean_mod_cal = mod_cal.mean_model(best_cal, obs, sps)
        mass_frac_smc = mean_mod_smc[-1]
        mass_frac_cal = mean_mod_cal[-1]
        mass_index = [i for i, s in enumerate(smc_thetas) if 'mass' in s]
        #metal_idx = [i for i, s in enumerate(best_smc) if 'logzsol' in s]
        mass_smc = best_smc[mass_index[0]] * mass_frac_smc
        mass_cal = best_cal[mass_index[0]] * mass_frac_cal
        bf_mass_smc.append(mass_smc)
        bf_mass_cal.append(mass_cal)
    
        galaxy_id.append(galaxy_num)
        #metallicity.append(theta_best[metal_idx[0]]) 


   
   
        tau_idx = smc_thetas.index('tau')
        tau_smc = best_smc[tau_idx]
        tau_cal = best_cal[tau_idx]
        age_idx = smc_thetas.index('tage')
        age_smc = best_smc[age_idx] #Gyr
        age_cal = best_cal[age_idx]
        time = 13.8 #np.arange(0, 13.8, 0.05)
        norm_smc = (mass_smc / mass_frac_smc) / (age_smc * 1.0e9)
        norm_cal = (mass_cal / mass_frac_cal) / (age_cal * 1.0e9)
        sfr_smc_now = norm_smc * (time-age_smc) * np.exp(-(time-age_smc) / tau_smc)
        sfr_cal_now = norm_cal * (time - age_cal) * np.exp(-(time-age_cal) / tau_cal)
        sfr_smc.append(sfr_smc_now)
        sfr_cal.append(sfr_cal_now)
        #age_mw = tau_massweighted_age(sfr_onegal, time)
        #mweighted_age.append(age_mw)

    


d = {'mass_calzetti': bf_mass_cal, 'mass_smc' : bf_mass_smc, 'log(sfr)_cal' : sfr_cal,
     'log(sfr)_smc' : sfr_smc, 'true_mass' : true_mass, 'true_log(sfr)' : true_sfr}

a = pd.DataFrame(data=d, index=galaxy_id)

a.to_pickle('/ufrc/narayanan/s.lower/atten_test/simba_runs/simba_bestfit_freedust'+str(free_dust)+'.pkl')





