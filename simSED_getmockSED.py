
# coding: utf-8

# In[6]:


import numpy as np
from hyperion.model import ModelOutput
from astropy.cosmology import Planck15
from astropy import units as u
from astropy import constants
import sys, os

np.set_printoptions(precision=2)


# In[4]:


#========================================================

#z = 0.0249 --> right now, snap305 is redshift 0
z = 1.007
run = sys.argv[1] + sys.argv[2]


#========================================================


# In[21]:

#print('On galaxy: ', sys.argv[2])
#try:
m = ModelOutput(run)
wav,flux = m.get_sed(inclination='all',aperture=-1)
#except: print('error loading galaxy: ',sys.argv[2])
wav  = np.asarray(wav)*u.micron #wav is in micron
wav *= (1.+z)

flux = np.asarray(flux)*u.erg/u.s
#dl = 10.0*u.pc
dl = Planck15.luminosity_distance(z)
dl = dl.to(u.cm)

    
flux /= (4.*3.14*dl**2.)
    
nu = constants.c.cgs/(wav.to(u.cm))
nu = nu.to(u.Hz)

flux /= nu
flux = flux.to(u.mJy)




#make a Prospector readable csv from the fluxes @ different filter positions
#error will be 5-10% of flux

filter_wavs = np.array([.15, .271, 0.33519, .475, .555, .606, .7499, .814, 1.058, 1.2516, 1.3969, 1.52359, 23.2096, 70., 100., 160., 250., 350., 500.]) #um 
flx = [None] * len(filter_wavs)
flxe = [None] * len(filter_wavs)
for i in range(0, len(filter_wavs)):
    flx[i] = flux[0][(np.abs(wav.value - filter_wavs[i])).argmin()].value
    flxe[i] = 0.01* flx[i]
flx = np.asarray(flx)
flxe = np.asarray(flxe)


# In[24]:


phot = np.insert(flxe, np.arange(len(flx)), flx)


# In[25]:

galex = ['galex_FUV']
hst_wfc3_uv  = ['wfc3_uvis_f275w', 'wfc3_uvis_f336w', 'wfc3_uvis_f475w','wfc3_uvis_f555w', 'wfc3_uvis_f606w', 'wfc3_uvis_f814w']
sdss = ['sdss_i0']
hst_wfc3_ir = ['wfc3_ir_f105w', 'wfc3_ir_f125w', 'wfc3_ir_f140w', 'wfc3_ir_f160w']
spitzer_mips = ['spitzer_mips_24']
#wise = ['wise_w4']                                                                                                                                                                  
herschel_pacs = ['herschel_pacs_70', 'herschel_pacs_100', 'herschel_pacs_160']
herschel_spire = ['herschel_spire_250', 'herschel_spire_350', 'herschel_spire_500']

filternames = (galex + hst_wfc3_uv + sdss + hst_wfc3_ir + spitzer_mips + herschel_pacs + herschel_spire)



#header info
#filter_names = ['galex_FUV', 'galex_FUV_err', 'wfc3_uvis_f275w', 'wfc3_uvis_f275w_err', 'wfc3_uvis_f336w', 'wfc3_uvis_f336w_err', 'wfc3_uvis_f475w', 
#                'wfc3_uvis_f475w_err','wfc3_uvis_f555w', 'wfc3_uvis_f555w_err', 'wfc3_uvis_f606w', 'wfc3_uvis_f606w_err', 'sdss_i0', 
#                'sdss_i0_err', 'wfc3_uvis_f814w',  'wfc3_uvis_f814w_err', 'wfc3_ir_f105w', 'wfc3_ir_f05w_err', 'wfc3_ir_f125w', 'wfc3_ir_f125w_err', 
#                'wfc3_ir_f140w', 'wfc3_ir_f140w_err', 'wfc3_ir_f160w', 'wfc3_ir_f160w_err',  'spitzer_mips_24', 'spitzer_mips_24_err',
#                'herschel_pacs_70', 'herschel_pacs_70_err', 'herschel_pacs_100', 'herschel_pacs_100_err', 
#                'herschel_pacs_160', 'herschel_pacs_160_err', 'herschel_spire_250', 'herschel_spire_250_err', 'herschel_spire_350',
#               'herschel_spire_350_err', 'herschel_spire_500', 'herschel_spire_500_err']

#id_ = str(os.path.splitext(sys.argv[2])[0])
id_ = sys.argv[2].split(os.extsep, 2)[1]

redshift = z


# In[26]:


line1 = ['id']
line1.extend(filternames)
line2 = [id_]
line2.extend(flx)
dat = np.stack((line1, line2))


# In[27]:


np.savetxt(sys.argv[1]+id_+'.csv', dat, delimiter=",", fmt="%s") 

