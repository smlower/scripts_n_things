import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from decimal import Decimal
import prospect.io.read_results as pread
import sys, os
import sedpy
from glob import glob
import pandas as pd
import numpy as np

gal = int(sys.argv[1])



def load_obs(sed_file):

    from hyperion.model import ModelOutput
    from astropy.cosmology import Planck15
    from astropy import units as u
    from astropy import constants
    m = ModelOutput(sed_file)

    wav,flux = m.get_sed(inclination='all',aperture=-1)
    wav  = np.asarray(wav)*u.micron #wav is in micron  
    wav = wav.to(u.AA)
    flux = np.asarray(flux)*u.erg/u.s
    dl = 10.0*u.pc
                                                                                                                              
    dl = dl.to(u.cm)
    flux /= (4.*3.14*dl**2.)
    nu = constants.c.cgs/(wav.to(u.cm))
    nu = nu.to(u.Hz)

    flux /= nu
    flux = flux.to(u.Jy)
    maggies = flux[0] / 3631.

    return maggies.value, wav



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

    theta_names = res.get('theta_labels', res['model'].theta_labels())
    return theta_names, theta_best



#tau_model = '/ufrc/narayanan/s.lower/simSEDs/simbam25n512_newfof/snap305_parametricSFH/tau/snap305.galaxy'+'{0:03}'.format(gal)+'*.h5'
#burst_model = '/ufrc/narayanan/s.lower/simSEDs/simbam25n512_newfof/snap305_parametricSFH/burst/snap305.galaxy'+'{0:03}'.format(gal)+'*_mcmc.h5'
cont_model = '/ufrc/narayanan/s.lower/simSEDs/simbam25n512_newfof/snap305_nonparametricSFH/filtered_cont/snap305.galaxy'+'{0:03}'.format(gal)+'*.h5'
dir_model = '/ufrc/narayanan/s.lower/simSEDs/simbam25n512_newfof/snap305_nonparametricSFH/filtered_dirichlet/snap305.galaxy'+'{0:03}'.format(gal)+'*.h5'
observed = '/ufrc/narayanan/s.lower/pd_runs/simba_m25n512/snap305_boxtest/snap305/snap305.galaxy'+'{0:03}'.format(gal)+'.rtout.sed'


tau_pkl = '/ufrc/narayanan/s.lower/simSEDs/simbam25n512_newfof/snap305_parametricSFH/tau/prospectorsnap305_properties_ParametricTau.pkl'

tau_data = pd.read_pickle(tau_pkl)

#dust1 = '/ufrc/narayanan/s.lower/simbam25n512_newfof/snap305_nonparametricSFH/dust1_toggle/snap305.galaxy'+'{0:03}'.format(gal)+'*_mcmc.h5'
dust1 = '/ufrc/narayanan/s.lower/simSEDs/simbam25n512_newfof/snap305_nonparametricSFH/dust1_toggle/snap305.galaxy'+'{0:03}'.format(gal)+'*_mcmc.h5'

obs_spec, obs_wave = load_obs(observed)
for file in glob(dust1):
    res_t, obs_t, mod_t = pread.results_from(file)
dust_sps = pread.get_sps(res_t)

#for file in glob(burst_model):
#    res_b, obs_b, mod_b = pread.results_from(file)
#burst_sps = pread.get_sps(res_b)

for file in glob(cont_model):
    res_c, obs_c, mod_c = pread.results_from(file)
cont_sps = pread.get_sps(res_c)

for file in glob(dir_model):
    res_d, obs_d, mod_d = pread.results_from(file)
dir_sps = pread.get_sps(res_d)


dust_best = get_best(res_t)[1]
#burst_best = get_best(res_b)[1]
cont_best = get_best(res_c)[1]
dir_best = get_best(res_d)[1]

phot_wavelengths = [x.wave_mean for x in obs_t['filters']]
spec_wavelengths = dust_sps.wavelengths

dust_spec, dust_phot, dust_massfrac = mod_t.mean_model(dust_best, obs_t, sps=dust_sps)
#burst_spec, burst_phot, burst_massfrac = mod_b.mean_model(burst_best, obs_b, sps=burst_sps)
cont_spec, cont_phot, cont_massfrac = mod_c.mean_model(cont_best, obs_c, sps=cont_sps)
dir_spec, dir_phot, dir_massfrac = mod_d.mean_model(dir_best, obs_d, sps=dir_sps)


plt.figure(figsize=(10, 8))
plt.tick_params(axis='both', which='major', labelsize=15)
plt.tick_params(axis='both', which='minor')

plt.scatter(phot_wavelengths, obs_t['maggies'], color='black', s=100, marker='+', zorder=10, label='Observed Phot')
plt.plot(obs_wave, obs_spec, color='gray', lw=1, label='True Spec')

plt.plot(spec_wavelengths, dust_spec, color='midnightblue' , lw=2, label=r'Dirichlet+dust1=0 Best Fit Spec')

#plt.plot(spec_wavelengths, burst_spec, color='turquoise' , ls='--', lw=2, label=r'$\tau$+burst Best Fit Spec')

plt.plot(spec_wavelengths, cont_spec, color='crimson' , lw=2, label='Cont Best Fit Spec')

plt.plot(spec_wavelengths, dir_spec, color='darkorange' , ls='--', lw=2, label='Dirichlet Best Fit Spec')

plt.annotate('Galaxy '+str(gal)+'\nM* ='+"{:.3e}".format(Decimal(str(tau_data.loc[['{0:03}'.format(gal)]]['Intrinsic Stellar Mass'][0] / 1.989e33)))+' M$_{\odot}$',xy=(1e1, 1e8), fontsize=15)


plt.loglog()
plt.ylim([1e3, 1e11])
#plt.xlim([1e2, 1e8])
plt.legend(loc='lower right', fontsize=15)
plt.ylabel('Flux [maggies]', fontsize=23)
plt.xlabel('Wavelength [$\AA$]', fontsize=23)

plt.savefig('/ufrc/narayanan/s.lower/SEDs/dust1/SED_'+str(gal)+'_OG.png', fpi=300, bbox_inches='tight')
