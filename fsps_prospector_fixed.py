#!/usr/bin/env python
# coding: utf-8

import numpy as np
from sedpy.observate import load_filters
from sedpy.observate import Filter
from prospect.fitting import fit_model
from prospect.io import write_results as writer
import os, sys

run_params = {'verbose':False,
              'debug':False,
              'output_pickles': False,
              # dynesty Fitter parameters
              'nested_bound': 'multi', # bounding method
              'nested_sample': 'auto', # sampling method
              'nested_nlive_init': 300,
              'nested_nlive_batch': 300,
              'nested_bootstrap': 0,
              'nested_dlogz_init': 0.05,
              'nested_weight_kwargs': {"pfrac": 1.0},
              }



# --------------
# Model Definition
# --------------

def build_model(**kwargs):

    from prospect.models.templates import TemplateLibrary
    from prospect.models import priors, sedmodel


    model_params = TemplateLibrary["parametric_sfh"]
    model_params.update(TemplateLibrary["dust_emission"])
    
    #fixed delayed-tau SFH
    model_params['tau']['isfree'] = False
    model_params['tau']['init'] = 1.0
    model_params['sf_start'] = {"N": 1, "isfree": False,
                           "init": 1.0, 
                            "units": "start of SF, Gyr"}
    model_params["lumdist"] = {"N": 1, "isfree": False,
                         "init": 1.0e-5, "units": "Mpc"}
    model_params['tage']['prior'] = priors.TopHat(mini=1.0, maxi=10.0)
    model_params['tage']['init'] = 5.0
    model_params['mass']['prior'] = priors.TopHat(mini=1e7, maxi=1e13)
    model_params['logzsol']['init'] = 0.2
    model_params['logzsol']['isfree'] = True
    model_params['logzsol']['prior'] = priors.ClippedNormal(mini=-1.5, maxi=0.5, mean=0.0, sigma=0.3)                            
    
    #dust emission
    model_params['duste_gamma']['init'] = 0.01
    model_params['duste_qpah']['init'] = 3.5
    model_params['duste_umin']['init'] = 1.0                                                                                                 
    model_params = model_library(model_params, int(sys.argv[1]), False)

    model = sedmodel.SedModel(model_params)
    

    return model


# --------------
# OBS
# --------------
galex = ['galex_FUV']
hst_wfc3_uv  = ['wfc3_uvis_f336w', 'wfc3_uvis_f475w','wfc3_uvis_f555w',  'wfc3_uvis_f814w']
hst_wfc3_ir = ['wfc3_ir_f125w', 'wfc3_ir_f160w']
irac = ['spitzer_irac_ch1']
herschel_pacs = ['herschel_pacs_70', 'herschel_pacs_100', 'herschel_pacs_160']
herschel_spire = ['herschel_spire_250', 'herschel_spire_350', 'herschel_spire_500']

filternames = (galex + hst_wfc3_uv + hst_wfc3_ir + irac + herschel_pacs + herschel_spire)


def build_obs(**kwargs):
    sps = build_sps()
    mod = build_model()
    fake_obs = {'filters': load_filters(filternames), 'wavelength': None}

    
    mod.params['dust_type'] = 0 # non-KC attenuation curve
    #mod.params['dust1'] = 1.0
    mod.params['dust_index'] = -0.7
    initial_theta = mod.initial_theta.copy()
    initial_theta[mod.theta_labels().index('dust2')] = 1.5
    initial_theta[mod.theta_labels().index('logzsol')] = 0.0
    initial_theta[mod.theta_labels().index('tage')] = 2.0
    spec, mags, stellar_mass = mod.mean_model(initial_theta,sps=sps,obs=fake_obs)

    obs = {}
    obs['maggies'] = mags  
    obs['filters'] = load_filters(filternames)
    obs['maggies_unc'] = mags * 0.1
    obs['phot_mask'] = np.isfinite(mags)
    obs['wavelength'] = None
    obs['spectrum'] = None
    return obs


def build_sps(zcontinuous=1, compute_vega_mags=False, **extras):
    from prospect.sources import CSPSpecBasis
    sps = CSPSpecBasis(zcontinuous=zcontinuous,
                       compute_vega_mags=compute_vega_mags)
    return sps

def build_noise(**extras):
    return None, None

def build_all(**kwargs):
    return (build_obs(**kwargs), build_model(**kwargs),
            build_sps(**kwargs), build_noise(**kwargs))






def model_library(model_params, attenuation_model, dust1=False, **kwargs):
    print("In library: ",attenuation_model)
    from prospect.models import priors
    if attenuation_model == 0:
        name = 'MW'
        model_params["dust_type"] = {"N": 1, "isfree": False, "init": 1, "units": name}
        model_params["mrw"] = {"N": 1, "isfree": False,
                           "init": 3.1, "units": "The ratio of total to selective absorption which characterizes the MW extinction curve", 'prior' : priors.TopHat(mini=2.0, maxi=5.0)}
        model_params['uvb'] = {'N': 1, 'isfree': True,
                           'init': 1.0, 'units':"characterizing the strength of the 2175A extinction feature with respect to the standard Cardelli et al. determination for the MW",
                           'prior': priors.TopHat(mini = 0.1, maxi=3.0)}
        model_params["dust1"]  = {"N": 1, "isfree": dust1,
                     "init": 0.0, "prior" : priors.TopHat(mini=0.0, maxi=1.5),
                     "units": "optical depth towards young stars"}

    elif attenuation_model == 1:
        name = 'PL'
        model_params["dust_type"] = {"N": 1, "isfree": False, "init": 0, "units": name} #power-law          
        model_params["dust2"]["prior"] = priors.TopHat(mini=0.0, maxi=2.0)
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
        name = 'Calzetti'
        print('got here')
        model_params["dust_type"] = {"N": 1, "isfree": False, "init": 2, "units": name} #Calzetti
        model_params["dust2"]["prior"] = priors.TopHat(mini=0.0, maxi=2.0)
        model_params["dust1"] = {"N": 1, "isfree": False,
                         "init": 0.0, "units": "optical depth towards young stars"}

    elif attenuation_model == 3:
        name = 'KC'
        model_params["dust_type"] = {"N": 1, "isfree": False, "init": 4, "units": name} #Conroy & Kriek dust atten.                                                           
        model_params["dust2"]["prior"] = priors.TopHat(mini=0.0, maxi=2.0)
        model_params["dust1"] = {"N": 1, "isfree": dust1,
                     "init": 0.0, "units": "optical depth towards young stars",
                     "prior": priors.TopHat(mini=0.0, maxi=2.0)}
        model_params["dust_index"] = {"N": 1, "isfree": True,
                     "init": -0.7, "units": "power-law multiplication of Calzetti",
                     "prior": priors.TopHat(mini=-2.5, maxi=0.7)}

    return model_params





if __name__ == '__main__':

    print('Building model, obs, and sps')
    obs, model, sps, noise = build_all(**run_params)
    name = model.init_config['dust_type']['units']
    print('Attenuation model:',name)
    run_params["sps_libraries"] = sps.ssp.libraries
    hfile = "/ufrc/narayanan/s.lower/atten_test/inputPowerLaw_fixedmodel%s.h5"%(name)
    print(hfile)
    print('Running fits')
    output = fit_model(obs, model, sps, noise, **run_params)
    print('Done. Writing now')
    writer.write_hdf5(hfile, run_params, model, obs,
                  output["sampling"][0], output["optimization"][0],
                  tsample=output["sampling"][1],
                  toptimize=output["optimization"][1])

    hfile.close()



