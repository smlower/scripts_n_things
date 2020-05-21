import numpy as np
import prospect.io.read_results as pread
import pandas as pd
import tqdm
import matplotlib
matplotlib.use('Agg')
from corner import quantile


def find_nearest(array,value):
    idx = (np.abs(np.array(array)-value)).argmin()
    return idx



#global attn curve params
dd63=6300.00
lamv=3000.0
dlam=350.0
lamuvb=2175.0

def Kriek_Conroy(lam, dust2, dust_index): 
    w63 = find_nearest(lam,dd63)
    cal00 = np.zeros(np.shape(lam)[0])
    for i in range(w63, np.shape(lam)[0]):
        cal00[i] = 1.17*( -1.857+1.04*(1.0e4/lam[i])) + 1.78 
    for i in range(0, w63):
        cal00[i]= 1.17*(-2.156+1.509*(1.0e4/lam[i]) -0.198*(1.0e4/lam[i])**2 + 0.011*(1.0e4/lam[i])**3) + 1.78
    cal00 = (cal00/0.44/4.05)
    eb = 0.85 - (1.9 * dust_index)  #KC13 Eqn 3                                                                                               
    #Drude profile for 2175A bump                                                                                                             
    drude = eb*(lam*dlam)**2 / ( (lam**2-lamuvb**2)**2 + (lam*dlam)**2 )

    attn_curve = dust2*(cal00+drude/4.05)*(lam/lamv)**dust_index

    return attn_curve

wav = np.linspace(100, 15000, 10000)

atten_ratio = []
atten_gal = []
for i in tqdm.tqdm(range(1900)):
    res_dir = 0
    try:
        dir_file = '/ufrc/narayanan/s.lower/simSEDs/simbam25n512_newfof/prod_runs/dirichlet/fixed_PAH/galaxy'+"{:03d}".format(i)+'.h5'
        res_dir, _, _ = pread.results_from(dir_file)
        pd_dat = np.load('/ufrc/narayanan/s.lower/pd_runs/simba_m25n512/snap305/attenuation_curves/042020/attenuation_curve.305_galaxy'+str(i)+'.npz')
    
    except: 
        continue
    
    atten_gal.append(i)
    dust2_chain = [item[res_dir['theta_labels'].index('dust2')]  for item in res_dir['chain']]
    #dustindex_chain = [item[res_dir['theta_labels'].index('dust_index')]  for item in res_dir['chain']]
    dust2_50 = quantile(dust2_chain, [.5])[0]
    #dustindex_50 = quantile(dustindex_chain, [.5])[0]
    #curve_dir = Kriek_Conroy(wav, dust2_50, dustindex_50)
    #filt_wav = pd_dat['wav_rest'] 
    #filt_tau = pd_dat['tau']
    filt_tauv = pd_dat['tau_v']
    #real_curve = (2.5*np.log10(np.exp(filt_tau[0])))
    #real_5500 = real_curve[find_nearest(filt_wav, 5500)]
    #prosp_5500 = curve_dir[find_nearest(wav, 5500)]
    atten_ratio.append(dust2_50 - filt_tauv)


np.savez('//ufrc/narayanan/s.lower/pd_runs/simba_m25n512/snap305/attenuation_curves/042020/dust2_diff_04222020.npz', galaxy=atten_gal, dat=atten_ratio)
