#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
from hyperion.model import ModelOutput
from astropy import constants
import astropy.units as u
import prospect.io.read_results as pread
import glob

# In[2]:


def find_nearest(array,value):
    idx = (np.abs(np.array(array)-value)).argmin()
    return idx

# In[3]:


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

# In[6]:


#attenuation from prospector

gal = 287

dir_file = '/Volumes/ufrc/narayanan/s.lower/simSEDs/simbam25n512_newfof/snap305_nonparametricSFH/filtered_dirichlet/snap305.galaxy'+"{:03d}".format(gal)+'_*.h5'
dirmm_file = '/Volumes/ufrc/narayanan/s.lower/simSEDs/simbam25n512_newfof/snap305_massmetal/dirichlet/snap305.galaxy'+"{:03d}".format(gal)+'_*.h5'
#contmm_file = '/Volumes/ufrc/narayanan/s.lower/simSEDs/simbam25n512_newfof/snap305_massmetal/continuity/snap305.galaxy'+"{:03d}".format(gal)+'_*.h5'
#cont_file = '/Volumes/ufrc/narayanan/s.lower/simSEDs/simbam25n512_newfof/snap305_nonparametricSFH/filtered_cont/snap305.galaxy'+"{:03d}".format(gal)+'_*.h5'
windowz_file = '/Volumes/ufrc/narayanan/s.lower/simSEDs/simbam25n512_newfof/snap305_metal_test/snap305.galaxy'+"{:03d}".format(gal)+'_*.h5'
dust1_file = '/Volumes/ufrc/narayanan/s.lower/simSEDs/simbam25n512_newfof/snap305_nonparametricSFH/dust1_toggle/snap305.galaxy'+"{:03d}".format(gal)+'_*.h5'
tau_file = '/Volumes/ufrc/narayanan/s.lower/simSEDs/simbam25n512_newfof/snap305_parametricSFH/tau/snap305.galaxy'+"{:03d}".format(gal)+'_*.h5'


for file in glob.glob(dir_file):
    res_dir, _, _ = pread.results_from(file)

for file in glob.glob(windowz_file):
    res_windowz, _, _ = pread.results_from(file)

for file in glob.glob(dirmm_file):
    res_dirmm, _, _ = pread.results_from(file)
    
#for file in glob.glob(contmm_file):
#    res_contmm, _, _ = pread.results_from(file)
    
for file in glob.glob(tau_file):
    res_tau, _, _ = pread.results_from(file)
    
for file in glob.glob(dust1_file):
    res_dust1, _, _ = pread.results_from(file)

# In[35]:


pd_atten = np.load('/Volumes/ufrc/narayanan/s.lower/pd_runs/simba_m25n512/snap305_boxtest/attenuation_curves/attenuation_curve.305_galaxy12.npz')

# In[36]:


pd_wav = pd_atten['wav_rest']
pd_tau = pd_atten['tau']

pd_V = find_nearest(pd_wav, 0.55)

# In[37]:


plt.figure(figsize=(10, 8))
plt.tick_params(axis='both', which='major', labelsize=15)
plt.tick_params(axis='both', which='minor')
plt.plot(pd_wav*1.0e4, (1.086*pd_tau[0])/(1.086*pd_tau[0][find_nearest(pd_wav, 0.55)]), color='black', lw=1., label='True Attenuation')



#plt.xscale('log')
plt.xlim([900, 7000])
plt.xlabel('$\lambda$ [$\AA$]', fontsize=23)
plt.ylim([0, 5])
plt.ylabel('A$_{\lambda}$ / A$_{V}$', fontsize=23)
plt.legend(loc='best', fontsize=15)
plt.title('Galaxy '+str(gal)+': KC Attenuation', fontsize=23)


#plt.savefig('/Users/sidneylower/Documents/prosp_final_plots/attenuation_curve'+str(gal)+'.png', dpi=300)

# In[9]:


sps_dir = pread.get_sps(res_dir)
thetas_dir, theta_best_dir = get_best(res_dir)
dust2_idx_dir = [i for i, s in enumerate(thetas_dir) if 'dust2' in s]
dust2_dir = theta_best_dir[dust2_idx_dir[0]]
dustindex_idx_dir = [i for i, s in enumerate(thetas_dir) if 'dust_index' in s]
dust_index_dir = theta_best_dir[dustindex_idx_dir[0]]


# In[10]:


sps_dust1 = pread.get_sps(res_dust1)
thetas_dust1, theta_best_dust1 = get_best(res_dust1)
dust2_idx_dust1 = [i for i, s in enumerate(thetas_dust1) if 'dust2' in s]
dust2_dust1 = theta_best_dust1[dust2_idx_dust1[0]]
dustindex_idx_dust1 = [i for i, s in enumerate(thetas_dust1) if 'dust_index' in s]
dust_index_dust1 = theta_best_dust1[dustindex_idx_dust1[0]]

# In[11]:


sps_windowz = pread.get_sps(res_windowz)
thetas_windowz, theta_best_windowz = get_best(res_windowz)
dust2_idx_windowz = [i for i, s in enumerate(thetas_windowz) if 'dust2' in s]
dust2_windowz = theta_best_dir[dust2_idx_windowz[0]]
dustindex_idx_windowz = [i for i, s in enumerate(thetas_windowz) if 'dust_index' in s]
dust_index_windowz = theta_best_windowz[dustindex_idx_windowz[0]]


# In[49]:


sps_cont = pread.get_sps(res_contmm)
thetas_cont, theta_best_cont = get_best(res_cont)
dust2_idx_cont = [i for i, s in enumerate(thetas_cont) if 'dust2' in s]
dust2_cont = theta_best_cont[dust2_idx_cont[0]]
dustindex_idx_cont = [i for i, s in enumerate(thetas_cont) if 'dust_index' in s]
dust_index_cont = theta_best_cont[dustindex_idx_cont[0]]

# In[12]:


sps_dirmm = pread.get_sps(res_dirmm)
thetas_dirmm, theta_best_dirmm = get_best(res_dirmm)
dust2_idx_dirmm = [i for i, s in enumerate(thetas_dirmm) if 'dust2' in s]
dust2_dirmm = theta_best_dirmm[dust2_idx_dirmm[0]]
dustindex_idx_dirmm = [i for i, s in enumerate(thetas_dirmm) if 'dust_index' in s]
dust_index_dirmm = theta_best_dirmm[dustindex_idx_dirmm[0]]

# In[51]:


sps_contmm = pread.get_sps(res_contmm)
thetas_contmm, theta_best_contmm = get_best(res_contmm)
dust2_idx_contmm = [i for i, s in enumerate(thetas_contmm) if 'dust2' in s]
dust2_contmm = theta_best_contmm[dust2_idx_contmm[0]]
dustindex_idx_contmm = [i for i, s in enumerate(thetas_contmm) if 'dust_index' in s]
dust_index_contmm = theta_best_contmm[dustindex_idx_contmm[0]]

# In[13]:


sps_tau = pread.get_sps(res_tau)
thetas_tau, theta_best_tau = get_best(res_tau)
dust2_idx_tau = [i for i, s in enumerate(thetas_tau) if 'dust2' in s]
dust2_tau = theta_best_tau[dust2_idx_tau[0]]
dustindex_idx_tau = [i for i, s in enumerate(thetas_tau) if 'dust_index' in s]
dust_index_tau = theta_best_tau[dustindex_idx_tau[0]]

# In[14]:


#global attn curve params
dd63=6300.00
lamv=5500.0
dlam=350.0
lamuvb=2175.0

# In[15]:


def Kriek_Conroy(lam, dust2, dust_index): 
    w63 = find_nearest(lam,dd63)
    cal00 = np.zeros(np.shape(lam)[0])
    for i in range(w63, np.shape(lam)[0]):
        cal00[i] = 1.17*( -1.857+1.04*(1.0e4/lam[i])) + 1.78 
    for i in range(0, w63):
        cal00[i]= 1.17*(-2.156+1.509*(1.0e4/lam[i]) -0.198*(1.0e4/lam[-1])**2 + 0.011*(1.0e4/lam[-1])**3) + 1.78
    #R=4.05 NB: I'm not sure I have this normalization correct...                                                                                            
    cal00 = (cal00/0.44/4.05)

    eb = 0.85 - (1.9 * dust_index)  #KC13 Eqn 3                                                                                                           

    #Drude profile for 2175A bump                                                                                                                            
    drude = eb*(lam*dlam)**2 / ( (lam**2-lamuvb**2)**2 + (lam*dlam)**2 )

    attn_curve = dust2*(cal00+drude/4.05)*(lam/lamv)**dust_index

    return attn_curve



# In[16]:


def Calzetti(lam, dust2):
    w63 = find_nearest(lam,dd63)
    cal00 = np.zeros(np.shape(lam)[0])
    for i in range(w63, np.shape(lam)[0]):
        cal00[i] = 1.17*( -1.857+1.04*(1.0e4/lam[i])) + 1.78
    for i in range(0, w63):
        cal00[i]= 1.17*(-2.156+1.509*(1.0e4/lam[i]) -0.198*(1.0e4/lam[-1])**2 + 0.011*(1.0e4/lam[-1])**3) + 1.78
    #R=4.05 NB: I'm not sure I have this normalization correct...                                                                                            
    cal00 = cal00/0.44/4.05
    
    attn_curve = cal00 * dust2
    
    return attn_curve

# In[17]:


wav = np.linspace(100, 15000, 10000)

# In[18]:


curve_dir = Kriek_Conroy(wav, dust2_dir, dust_index_dir)
curve_dirmm = Kriek_Conroy(wav, dust2_dirmm, dust_index_dirmm)
curve_windowz = Kriek_Conroy(wav, dust2_windowz, dust_index_windowz)
#curve_contmm = Kriek_Conroy(wav, dust2_contmm, dust_index_contmm)
#curve_dust1 = Kriek_Conroy(wav, dust2_dust1, dust_index_dust1)

curve_tau = Kriek_Conroy(wav, dust2_tau, dust_index_tau)



# In[22]:


plt.figure(figsize=(10, 8))
plt.tick_params(axis='both', which='major', labelsize=15)
plt.tick_params(axis='both', which='minor')
plt.plot(wav, curve_dir, color='darkorange', lw=3, label='Dirichlet')
plt.plot(wav, curve_dirmm, color='darkorange', alpha=0.7, ls='--', lw=3, label='M*-Z + Dirichlet')
plt.plot(wav, curve_windowz, color='darkorange', alpha=0.7, lw=3, ls=':', label='Pseudo-fixed Z + Dirichlet')
#plt.plot(wav, curve_dust1, color='darkorange', alpha=0.7, lw=3, ls='-.',label='Dust1=0 + Dirichlet')
plt.plot(wav, curve_tau, color='midnightblue', alpha=0.7, lw=3, label=r'$\tau$')
plt.plot(pd_wav*1.0e4, (1.086*pd_tau[0])/(1.086*pd_tau[0][find_nearest(pd_wav, 0.55)]), color='black', lw=1., label='True Attenuation')



#plt.xscale('log')
plt.xlim([900, 7000])
plt.xlabel('$\lambda$ [$\AA$]', fontsize=23)
plt.ylim([0, 5])
plt.ylabel('A$_{\lambda}$ / A$_{V}$', fontsize=23)
plt.legend(loc='best', fontsize=15)
plt.title('Galaxy '+str(gal)+': KC Attenuation', fontsize=23)


#plt.savefig('/Users/sidneylower/Documents/prosp_final_plots/attenuation_curve'+str(gal)+'.png', dpi=300)

# In[91]:


dust2_contmm, dust2_dirmm

# In[39]:


from hyperion.model import ModelOutput

# In[61]:


comp_sed = ModelOutput('/Volumes/ufrc/narayanan/s.lower/pd_runs/simba_m25n512/snap305_boxtest/mono_seds/snap305/snap305.galaxy005.rtout.sed')
wav_rest_sed,dum_lum_obs_sed = comp_sed.get_sed(inclination='all',aperture=-1)

# In[56]:


plt.loglog(wav_rest_sed, dum_lum_obs_sed[0])

# In[ ]:



