import fsps
import numpy as np
import prospect.io.read_results as pread
import sys, os, glob
from prospector_output_utilities import *
from prospect.models import transforms
from astropy.cosmology import Planck15
from astropy import units as u


def find_nearest(array,value):
    idx = (np.abs(np.array(array)-value)).argmin()
    return idx

galaxy = int(sys.argv[1])

gal_num = "{:03d}".format(galaxy)

res, obs, mod = 0, 0, 0

infile = '/ufrc/narayanan/s.lower/simSEDs/simbam25n512_newfof/snap305_ultra_mist/dirichlet/sfh3/snap305.galaxy'+gal_num+'_*.h5'

#prosp_prosp = '/ufrc/narayanan/s.lower/simSEDs/simbam25n512_newfof/snap305_ultra_mist/

for filee in glob.glob(infile):
    print(filee)
    res, obs, mod = pread.results_from(filee)

thetas, best = get_best(res)

sps = pread.get_sps(res)

agebins = res['model'].params['agebins']

Z_idx = [i for i, s in enumerate(thetas) if 'massmet_2' in s]
zfrac_idx = [i for i, s in enumerate(thetas) if 'z_fraction' in s]
zfracs = best[zfrac_idx]
Z = best[Z_idx]
mass_idx = [i for i, s in enumerate(thetas) if 'massmet_1' in s]
mass = best[mass_idx]

print('Z', Z)
print('zfracs', zfracs)
print('mass', mass)

spec, phot, mass_frac = mod.mean_model(best, obs, sps)

#st_mass = (10**mass) * mass_frac

agelims = np.array([(10**x) / 1e9 for x in np.unique(np.ravel(agebins))])
sfrs = transforms.zfrac_to_sfr(mass, zfracs, agebins)

#sfrs.insert(0,sfrs[0])
sfrs = np.insert(sfrs, 0, sfrs[0])
#agelims1 = agelims[::-1]
#sfrs1 = sfrs[::-1]


sp_nodust = fsps.StellarPopulation(zcontinuous=1, add_dust_emission=False, logzsol=Z, sfh=3, dust_type=2, dust2=0)


sp_nodust.set_tabular_sfh(agelims, sfrs)

wav, lum_hz_permass = sp_nodust.get_spectrum(tage=agelims.max())

lum_hz = lum_hz_permass * (10**mass)

W_hz = lum_hz * 3.827e26  #now W/Hz                                                    
z = 0.01
dl = Planck15.luminosity_distance(z)
dl = dl.to(u.m)

conversion = (1./(4.*np.pi*(dl.value)**2)) * (1.0e26)
to_flux = W_hz * conversion #now Jy
to_maggies = to_flux / 3631.

stellar_spec = to_maggies
obs_spec = spec

vband_idx = find_nearest(wav,3000)

vband_extinction = obs_spec[vband_idx]/stellar_spec[vband_idx]


extinction = np.zeros(len(spec))
for i in range(spec.shape[0]):
    extinction[i] = obs_spec[i] / stellar_spec[i]


tau = -1 * np.log(extinction)
tau_v = -1 * np.log(vband_extinction)


np.savez('/ufrc/narayanan/s.lower/simSEDs/simbam25n512_newfof/snap305_ultra_mist/dirichlet/sfh3/dir_atten_galaxy'+gal_num+'.npz', tau=tau, tau_v=tau_v, wave=wav)
