import numpy as np
import prospect.io.read_results as pread
import fsps
from prospect.models.transforms import zfrac_to_sfr
from corner import quantile
import sys, os
import yt

#prosp_dir = '/ufrc/narayanan/s.lower/simSEDs//simbam25n512_newfof/experiments/dust_screen/dirichlet/'
#pd_dir = '//ufrc/narayanan/s.lower/pd_runs/simba_m25n512/snap305_dustscreen/snap305/'

#prosp_dir = '/ufrc/narayanan/s.lower/simSEDs//simbam25n512_newfof/prod_runs/tau/final2/'
prosp_dir = '/ufrc/narayanan/s.lower/simSEDs/simbam25n512_newfof/prod_runs/tau/no_agb/'
galaxy = int(sys.argv[1])
#appending path of prosp run script to import model and sps                                                                                                 

print('loading model and sps')
sys.path.append(prosp_dir)
from run_prosp import build_model, build_sps

mod = build_model()
sps = build_sps()

output_file = np.load(prosp_dir+'/galaxy_'+str(galaxy)+'.npz', allow_pickle=True)

dirich_timelist = np.unique(np.array([ 1.38000000e+01,  1.37000000e+01,  1.37000000e+01,  1.34688689e+01,
        1.34688689e+01,  1.32487758e+01,  1.32487758e+01,  1.28823933e+01,
        1.28823933e+01,  1.22724872e+01,  1.22724872e+01,  1.12571945e+01,
        1.12571945e+01,  9.56706658e+00,  9.56706658e+00,  6.75356055e+00,
        6.75356055e+00,  2.07, 0.]))

#sfr50 = output_file['data'][()]['SFR_50']
#print(sfr50)
#massfrac = output_file['data'][()]['Massfrac']
#mass50 = output_file['data'][()]['Mass_50']
true_mstar = output_file['data'][()]['True Mass'] / 1.989e33        



print('loading prosp results')
res, _, _ = pread.results_from(prosp_dir+'/galaxy'+str(galaxy)+'.h5')
thetas = mod.theta_labels()
print(thetas)
imax = np.argmax(res['lnprobability'])
theta_max = res["chain"][imax, :]

smass_idx = [i for i, s in enumerate(thetas) if 'massmet_1' in s][0]
smass = theta_max[smass_idx]
duste_gamma_idx = [i for i, s in enumerate(thetas) if 'duste_gamma' in s][0]
gamma = theta_max[duste_gamma_idx]#[item[duste_gamma_idx] for item in res['chain']]
duste_umin_idx = [i for i, s in enumerate(thetas) if 'duste_umin' in s][0]
umin = theta_max[duste_umin_idx]#[item[duste_umin_idx] for item in res['chain']]
#duste_qpah_idx = [i for i, s in enumerate(thetas) if 'duste_qpah' in s][0]
#pah = theta_max[duste_qpah_idx]#[item[duste_qpah_idx] for item in res['chain']]
dust2_idx = [i for i, s in enumerate(thetas) if 'dust2' in s][0]
d2 = theta_max[dust2_idx]#[item[dust2_idx] for item in res['chain']]
dust_index_idx = [i for i, s in enumerate(thetas) if 'dust_index' in s][0]
dindex = theta_max[dust_index_idx]#[item[dust_index_idx] for item in res['chain']]
logzsol_idx = [i for i, s in enumerate(thetas) if 'massmet_2' in s][0]
Z = theta_max[logzsol_idx]#[item[logzsol_idx] for item in res['chain']]


tau_idx = [i for i, s in enumerate(thetas) if 'tau' in s][0]
tage_idx = [i for i, s in enumerate(thetas) if 'tage' in s][0]
tau = theta_max[tau_idx]
tage = theta_max[tage_idx]


'''
Z = quantile(logzsol, [.5])
dindex = quantile(dust_index, [.5])
d2 = quantile(dust2, [.5])
umin = quantile(duste_umin, [.5])
pah = quantile(duste_qpah, [.5])
gamma = quantile(duste_gamma, [.5])
'''




#print('SFH')
'''
zfrac_idx = [i for i, s in enumerate(thetas) if 'z_fraction' in s]
zfrac_chain = [item[zfrac_idx[0]:zfrac_idx[-1]+1]  for item in res['chain']]
total_mass = (10**smass)
time_bins_log = next(item for item in res['model_params'] if item["name"] == "agebins")['init']
sfr_chain = []
for i in range(len(zfrac_chain)):
    sfr_chain.append(zfrac_to_sfr(total_mass, zfrac_chain[i], time_bins_log))
    
sfr50 = []
for i in np.arange(9, -1, -1):
    sfr_quantiles = quantile([item[i] for item in sfr_chain], [.5])
    sfr50.append(sfr_quantiles)
'''


#print(sfr50)
    
print('generating SPS')
sp = fsps.StellarPopulation(zcontinuous=1, dust_type=4, 
                            add_agb_dust_model=False, sfh=4, tau=tau, tage=tage)
'''
timelist = []
for i in range(len(dirich_timelist)-1):
    timelist.append((dirich_timelist[i] + dirich_timelist[i+1])/2)

#print(sfr50, timelist)

dirich_timelist = np.unique(np.array([ 1.38000000e+01,  1.37000000e+01,  1.37000000e+01,  1.34688689e+01,
        1.34688689e+01,  1.32487758e+01,  1.32487758e+01,  1.28823933e+01,
        1.28823933e+01,  1.22724872e+01,  1.22724872e+01,  1.12571945e+01,
        1.12571945e+01,  9.56706658e+00,  9.56706658e+00,  6.75356055e+00,
        6.75356055e+00,  2.07, 0.]))

dir_binned = np.array([])
time_binned = np.array([])
for i in range(9):
    dir_binned = np.append(dir_binned, np.full(10,sfr50[i]))
    time_binned = np.append(time_binned, np.linspace(dirich_timelist[i], dirich_timelist[i+1]-0.01, 10))
'''
#print(time_binned)
#print(dir_binned)

#print(np.all(time_binned[1:] > time_binned[:-1]))

#sp.set_tabular_sfh(time_binned, dir_binned)


sp.params['duste_gamma'] = gamma
sp.params['duste_umin'] = umin
sp.params['duste_qpah'] = 5.86
sp.params['dust2'] = d2
sp.params['dust1'] = 0.0
sp.params['dust_index'] = dindex
sp.params['logzsol'] = Z


#print(sp.params)
#print(sp.params.iteritems())


inferred_dust_mass = sp.dust_mass * (10**(smass))


print('loading snapshot')
filt_snaps = '/orange/narayanan/s.lower/simba/desika_filtered_snaps/snap305/'
ds = yt.load(filt_snaps+'/galaxy_'+str(galaxy)+".hdf5")
data = ds.all_data()
print('getting dust mass')
dust_mass = data.ds.arr(data[("PartType0", "Dust_Masses")].value, 'code_mass')
true_dust_mass = np.sum(dust_mass.in_units('Msun').value)

print('inferred dust mass', inferred_dust_mass)
print('true dust mass', true_dust_mass)


np.savez(prosp_dir+'/galaxy'+str(galaxy)+'_dustmasses.npz', galaxy=galaxy, prosp_mass=inferred_dust_mass, true_mass=true_dust_mass, stellar_mass=true_mstar)
