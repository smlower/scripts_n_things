import prospect.io.read_results as pread
import numpy as np
import pandas
from glob import glob


smc_dir = '/ufrc/narayanan/s.lower/atten_test/simba_runs/SMC/'
cal_dir = '/ufrc/narayanan/s.lower/atten_test/simba_runs/calzetti/'
pd_dir = '/ufrc/narayanan/s.lower/pd_runs/simba_m25n512/snap305_boxtest/snap305'
caesar_dir = '/ufrc/narayanan/s.lower/caesar_sfr50.npz'

smc_mass = []
cal_mass = []
sfr_smc = []
sfr_cal = []
true_mass = []
true_sfr = []
galaxy_id = []




for galaxy in range(0, 1001):
    galaxy_num = "{:03d}".format(galaxy)
    caesar = np.load(caesar_dir, allow_pickle=True)
    smc, cal, pd = 0, 0, 0
    for infile in glob(smc_dir+'snap*_galaxy'+galaxy_num+'_*.h5'):
        smc, _, _ = pread.results_from(infile)
    for infile in glob(cal_dir+'snap*_galaxy'+galaxy_num+'_*.h5'):
        cal, _, _ = pread.results_from(infile)
    for infile in glob(pd_dir+'/grid_physical_properties.305_galaxy'+galaxy_num+'.npz'):
        pd = np.load(infile, allow_pickle=True)
    
    if isinstance(smc, int) == True or isinstance(cal, int) == True or isinstance(pd, int) == True:
        #print('smc: ',isinstance(smc, int))
        #print('cal: ',isinstance(cal, int))
        #print('pd: ',isinstance(pd, int))
        print('galaxy %s corrupted'%(galaxy))
        continue
    else:
        print('onto galaxy %s'%(galaxy))
        true_mass.append(np.sum(pd['grid_star_mass'])/1.989e33)
        true_sfr.append(np.log10(caesar['sfr_50'][galaxy]))
        
        msimba1 = [item[0] for item in smc['chain']]
        msimba2 = [item[0] for item in cal['chain']]
        
        smc_mass.append(np.median(msimba1))
        cal_mass.append(np.median(msimba2))

        tausimba1 = [item[4] for item in smc['chain']]
        tausimba2 = [item[4] for item in cal['chain']]

        agesimba1 = [item[3] for item in smc['chain']]
        agesimba2 = [item[3] for item in cal['chain']]
        
        sfrsmc = []
        for i in range(len(msimba1)):
            sfrsmc.append(np.log10((msimba1[i] / (agesimba1[i]*1.e9)) * (13.8-agesimba1[i]) *np.exp(-(13.8-agesimba1[i]) / tausimba1[i])))
 
        sfrcal = []

        for i in range(len(msimba2)):
            sfrcal.append(np.log10((msimba2[i] / (agesimba2[i]*1.e9)) * (13.8-agesimba2[i]) *np.exp(-(13.8-agesimba2[i]) / tausimba2[i])))
  

        sfr_smc.append(np.median(sfrsmc))
        sfr_cal.append(np.median(sfrcal))
        
        galaxy_id.append(galaxy)



d = {'mass_calzetti': cal_mass, 'mass_smc' : smc_mass, 'log(sfr)_cal' : sfr_cal, 
     'log(sfr)_smc' : sfr_smc, 'true_mass' : true_mass, 'true_log(sfr)' : true_sfr}

a = pandas.DataFrame(data=d, index=galaxy_id)

a.to_pickle('/ufrc/narayanan/s.lower/atten_test/simba_runs/simba_medians_freedust.pkl')
