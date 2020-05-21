from corner import quantile
import prospect.io.read_results as pread
from prospect.models.transforms import zfrac_to_sfrac, logsfr_ratios_to_sfrs, zfrac_to_sfr, tburst_from_fage, zfrac_to_masses, logsfr_ratios_to_masses
import numpy as np
import sys
import glob, os


def nonpara_massweighted_age(sfr_in_each_bin, time_bins):
    top = 0.0
    bottom = 0.0
    time = (10**time_bins) / 1e9

    for bin_ in range(len(sfr_in_each_bin)):
        top += np.abs(time[bin_] * sfr_in_each_bin[bin_])
        bottom += np.abs(sfr_in_each_bin)
    return top / bottom

prosp_dir = '/ufrc/narayanan/s.lower/simSEDs//simbam25n512_newfof/experiments/dust_screen/dirichlet/'
pd_dir = '//ufrc/narayanan/s.lower/pd_runs/simba_m25n512/snap305_dustscreen/snap305/'
galaxy = int(sys.argv[1])

#appending path of prosp run script to import model and sps
sys.path.append(prosp_dir)
from run_prosp import build_model, build_sps



sfr_50 = [] 
sfr_16 = []
sfr_84 = []

print('now reading files')
galaxy_num = "{:03d}".format(galaxy)
infile = prosp_dir+'/galaxy'+str(galaxy)+'.h5'
globfiles = glob.glob(infile)

for prosp_output in glob.glob(infile):
    print(prosp_output)
    res, obs, _ = pread.results_from(prosp_output)
pdfile = pd_dir+'/grid_physical_properties.305_galaxy'+str(galaxy)+'.npz'
pd_data = np.load(pdfile)
int_mass = np.sum(pd_data['grid_star_mass'])

print('sps and model')
sps = build_sps()
mod = build_model()
thetas = mod.theta_labels()

thetas_50 = []
thetas_16 = []
thetas_84 = []
print('quantiles for all thetas')
for theta in thetas:
    idx = thetas.index(theta)
    chain = [item[idx] for item in res['chain']]
    quan = quantile(chain, [.16, .5, .84])
    thetas_50.append(quan[1])
    thetas_16.append(quan[0])
    thetas_84.append(quan[2])


mod_50 = mod.mean_model(thetas_50, obs, sps)
massfrac_50 = mod_50[-1]
mod_16 = mod.mean_model(thetas_16, obs, sps)
massfrac_16 = mod_16[-1]
mod_84 = mod.mean_model(thetas_84, obs, sps)
massfrac_84 = mod_84[-1]

print('mass and Z')
mass = thetas_50[thetas.index('massmet_1')]
mass_50 = thetas_50[thetas.index('massmet_1')]
mass_16 = thetas_16[thetas.index('massmet_1')]
mass_84 = thetas_84[thetas.index('massmet_1')]
Z_50 = thetas_50[thetas.index('massmet_2')]
Z_16 = thetas_16[thetas.index('massmet_2')]
Z_84 = thetas_84[thetas.index('massmet_2')]

print('sfr')
total_mass = 10**mass
time_bins_log = next(item for item in res['model_params'] if item["name"] == "agebins")['init']
print('time: ',len(time_bins_log))
print(time_bins_log)
zfrac_idx = [i for i, s in enumerate(thetas) if 'z_fraction' in s]
#print(zfrac_idx)
zfrac_chain = [item[zfrac_idx[0]:zfrac_idx[-1]+1]  for item in res['chain']]
#print(np.shape(zfrac_chain))
sfr_chain = []
for i in range(len(zfrac_chain)):
    sfr_chain.append(zfrac_to_sfr(total_mass, zfrac_chain[i], time_bins_log))
print(np.shape(sfr_chain))
sfr_quan = []
for i in range(np.shape(zfrac_chain)[1]+1): 
    sfr_quan.append(quantile([item[i] for item in sfr_chain], [.16, .5, .84]))

sfr_50.append([item[1] for item in sfr_quan])
sfr_16.append([item[0] for item in sfr_quan])
sfr_84.append([item[2] for item in sfr_quan])
time = []
#for val in time_bins_log:
#    time.append((10**val / 1.0e9)) #Gyr


print('sfr:',len(sfr_50[0]))
print(sfr_50[0][::-1])

data = {'Galaxy' : galaxy,  'SFR_50' : sfr_50[0][::-1], 'SFR_16' : sfr_16[0][::-1],
        'SFR_84' : sfr_84[0][::-1], 'True Mass' : int_mass, 'Mass_50' : mass_50, 
        'Mass_16' : mass_16, 'Mass_84' : mass_84, 'Z_50' : Z_50, 'Z_16' : Z_16,
        'Z_84' : Z_84, 'Massfrac' : massfrac_50, 'Massfrac16' : massfrac_16, 'Massfrac84' : massfrac_84,}


np.savez(prosp_dir+'/galaxy_'+str(galaxy)+'.npz', data=data)
