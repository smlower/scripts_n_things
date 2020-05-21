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

model = sys.argv[2]

#prosp_dir = '/ufrc/narayanan/s.lower/simSEDs//simbam25n512_newfof/experiments/dust_screen/'+model+'/'
#pd_dir = '/ufrc/narayanan/s.lower/pd_runs/simba_m25n512/snap305_dustscreen/snap305/'
galaxy = int(sys.argv[1])


prosp_dir = '/ufrc/narayanan/s.lower/simSEDs//simbam25n512_newfof/prod_runs/'+model+'/no_agb/'
pd_dir = '/ufrc/narayanan/s.lower/pd_runs/simba_m25n512/snap305/mist_pd/snap305/'


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

try:
    for prosp_output in glob.glob(infile):
        print(prosp_output)
        res, obs, _ = pread.results_from(prosp_output)
    pdfile = pd_dir+'/grid_physical_properties.305_galaxy'+str(galaxy_num)+'.npz'
    pd_data = np.load(pdfile)
    int_mass = np.sum(pd_data['grid_star_mass'])
except:
    print('file not found')

print('sps and model')
sps = build_sps()
mod = build_model()
thetas = mod.theta_labels()
#print(mod)
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

mass_50 = thetas_50[thetas.index('massmet_1')]
mass_16 = thetas_16[thetas.index('massmet_1')]
mass_84 = thetas_84[thetas.index('massmet_1')]
Z_50 = thetas_50[thetas.index('massmet_2')]
Z_16 = thetas_16[thetas.index('massmet_2')]
Z_84 = thetas_84[thetas.index('massmet_2')]


if model == 'tau':
    age = [item[thetas.index('tage')] for item in res['chain']]
    tau = [item[thetas.index('tau')] for item in res['chain']]
    mass = [10**item[thetas.index('massmet_1')] for item in res['chain']]
    time = np.arange(0, 13.8, 0.1)
    norm = np.asarray(mass) / (np.asarray(age)*1.0e9)
    sfr = []
    for j in range(len(mass)):
        sfr_current = []
        for i in time:
            sfr_current.append(norm[j] * (i-age[j]) * np.exp(-(i-age[j]) / tau[j]))
        sfr.append(sfr_current)
    sfr_quan = []
    for i in range(len(time)):
        sfr_quan.append(quantile([item[i] for item in sfr], [.16, .5, .84]))
    sfr_50.append([item[1] for item in sfr_quan])
    sfr_16.append([item[0] for item in sfr_quan])
    sfr_84.append([item[2] for item in sfr_quan])


elif model == 'burst':
    print('made it to burst sfh section')
    
    imax = np.argmax(res['lnprobability'])
    theta_max = res["chain"][imax, :]
    age = theta_max[thetas.index('tage')]
    tau = theta_max[thetas.index('tau')]
    mass = theta_max[thetas.index('massmet_1')]
    print(mass)
    fburst = theta_max[thetas.index('fburst')]
    fage = theta_max[thetas.index('fage_burst')]
    sfr = []
    time = np.arange(0, 13.8, 0.1)
    tburst = tburst_from_fage(age, fage)
    norm_tau = ((1 - fburst) * 10**mass) / (age * 1.0e9)
    print(norm_tau)
    for j in time:
        if j < age:
            sfr.append(0)
        else:
            sfr.append(norm_tau * (j - age) * np.exp(-(j-age) / tau))
    time = np.append(time, tburst)
    time.sort()
    idx = np.where(time == tburst)[0][0]
    norm_burst = (fburst * mass) / (1e8)#((time[idx + 1] - time[idx - 1])*1.0e9) 
    sfr.insert(idx, norm_burst)
    #print(sfr)
    #print('sfr :', np.shape(sfr))
    #print('time: ',np.shape(time))
    #print(time)
    #sfr_quan = []
    #for i in range(len(time)):
    #    sfr_quan.append(quantile([item[i] for item in sfr], [.16, .5, .84]))
    #print(sfr_quan)
    #sfr_50.append([item[1] for item in sfr_quan])
    #sfr_16.append([item[0] for item in sfr_quan])
    #sfr_84.append([item[2] for item in sfr_quan])
    sfr_84, sfr_50, sfr_16 = [0], sfr, [0]

    print(sfr_50)
else:
    age = [item[thetas.index('tage')] for item in res['chain']]
    mass = [10**item[thetas.index('massmet_1')] for item in res['chain']]
    sfr = []
    for i in range(len(mass)):
        sfr.append(mass[i] / (age[i]*1.0e9))
    
    
    sfr_quan = quantile(sfr, [.16, .5, .84])
    sfr_50 = sfr_quan[1]
    sfr_16 = sfr_quan[0]
    sfr_84 = sfr_quan[2]

if (model == 'tau'):

    data = {'Galaxy' : galaxy, 'Time': np.unique(np.ravel(time)), 'SFR_50' : sfr_50[0], 'SFR_16' : sfr_16[0],
        'SFR_84' : sfr_84[0], 'True Mass' : int_mass, 'Mass_50' : mass_50, 
        'Mass_16' : mass_16, 'Mass_84' : mass_84, 'Z_50' : Z_50, 'Z_16' : Z_16,
            'Z_84' : Z_84, 'Massfrac' : massfrac_50, 'Massfrac16' : massfrac_16, 'Massfrac84' : massfrac_84}
elif model == 'burst':
    data = {'Galaxy' : galaxy, 'Time': time, 'SFR_50' : sfr_50, 'SFR_16' : sfr_16[0],
        'SFR_84' : sfr_84[0], 'True Mass' : int_mass, 'Mass_50' : mass_50,
        'Mass_16' : mass_16, 'Mass_84' : mass_84, 'Z_50' : Z_50, 'Z_16' : Z_16,
            'Z_84' : Z_84, 'Massfrac' : massfrac_50, 'Massfrac16' : massfrac_16, 'Massfrac84' : massfrac_84}
else:

    data = {'Galaxy' : galaxy, 'SFR_50' : sfr_50, 'SFR_16' : sfr_16,
        'SFR_84' : sfr_84, 'True Mass' : int_mass, 'Mass_50' : mass_50,
        'Mass_16' : mass_16, 'Mass_84' : mass_84, 'Z_50' : Z_50, 'Z_16' : Z_16,
            'Z_84' : Z_84, 'Massfrac' : massfrac_50, 'Massfrac16' : massfrac_16, 'Massfrac84' : massfrac_84}

np.savez(prosp_dir+'/galaxy_'+str(galaxy)+'.npz', data=data)
