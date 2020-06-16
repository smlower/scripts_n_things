#;from __future__ import print_function
#attenuation_curve_gen.py

#generate npz attenuation curves so that we don't have to go through
#the slow reading in of rtout files but once


import numpy as np
from hyperion.model import ModelOutput
from astropy.cosmology import Planck13
from astropy import units as u
from astropy import constants
import h5py
from hyperion_input_seds import get_input_seds
from glob import glob

import sys, os



#========================================================
#MODIFIABLE HEADER (make this a function later with argv)
#z = 0.001


galaxy_num = sys.argv[1]
snap_str = '305'

#composite SED
sed_directory = '/orange/narayanan/s.lower/simba/pd_runs/snap305/'
#stellar only SED
sources_directory = sed_directory


output_directory = '/orange/narayanan/s.lower/simba/pd_runs/attenuation_curves/snap305/'
#========================================================



def find_nearest(array,value):
    idx = (np.abs(np.array(array)-value)).argmin()
    return idx

def sed_downsample(source_lam,lum_source,wav_rest):
    source_lam_downsampled = np.zeros(len(wav_rest))
    lum_source_downsampled = np.zeros(len(wav_rest))
    for i in range(len(wav_rest)):
        idx = find_nearest(source_lam.value,wav_rest[i].value)
        source_lam_downsampled[i] = source_lam[idx].value
        lum_source_downsampled[i] = lum_source[idx].value

    return source_lam_downsampled,lum_source_downsampled


#check if path exists - if not, create it
#if not os.path.exists(output_directory): os.makedirs(output_directory)


sed_file = glob(sed_directory+'/*galaxy'+str(galaxy_num)+'.rtout.sed')[0]
#print(sed_file)

stellar_file = glob(sources_directory+'/*.galaxy'+str(galaxy_num)+'.rtout.sed')[0]




comp_sed = ModelOutput(sed_file)
wav_rest_sed,dum_lum_obs_sed = comp_sed.get_sed(inclination='all',aperture=-1)
wav_rest_sed =wav_rest_sed* u.micron #wav is in micron 
nu_rest_sed = constants.c.cgs/wav_rest_sed.cgs
lum_obs_sed = dum_lum_obs_sed
lum_obs_sed = lum_obs_sed * u.erg/u.s
nu_rest_sed = constants.c.cgs/(wav_rest_sed.to(u.cm))
fnu_obs_sed = lum_obs_sed.to(u.Lsun)
fnu_obs_sed /= nu_rest_sed.to(u.Hz)
fnu_obs_sed = fnu_obs_sed.to(u.Lsun/u.Hz)

#stellar_sed = np.load(stellar_file)
#nu_rest_stellar = stellar_sed['nu'] #Hz
#fnu_rest_stellar = stellar_sed['fnu'] #Lsun/Hz
#fnu_rest_stellar = fnu_rest_stellar * u.Lsun/u.Hz
#nu_rest_stellar = nu_rest_stellar * u.Hz
#lum_rest_stellar = fnu_rest_stellar * nu_rest_stellar
#lam_rest_stellar = constants.c.cgs/(nu_rest_stellar)

#source_lam_downsampled,lum_source_downsampled = sed_downsample(lam_rest_stellar.to(u#.micron),lum_rest_stellar.to(u.Lsun),wav_rest_sed.to(u.micron))
#source_lam_downsampled = source_lam_downsampled * u.micron
#lum_source_downsampled  = lum_source_downsampled * u.Lsun



source_nu,source_fnu = get_input_seds(stellar_file)
source_nu = source_nu * u.Hz
source_lam = constants.c.cgs/source_nu
source_lam = source_lam.to(u.micron)
source_fnu = source_fnu * u.Lsun/u.Hz
lum_source = source_fnu * source_nu


source_lam_downsampled,lum_source_downsampled = sed_downsample(source_lam.to(u.micron),lum_source.to(u.Lsun),wav_rest_sed.to(u.micron))
source_lam_downsampled = source_lam_downsampled * u.micron                                                                 
lum_source_downsampled  = lum_source_downsampled * u.Lsun



#stellar_sed = ModelOutput(stellar_file)
#wav_rest_stellar,dum_lum_stellar_sed = stellar_sed.get_sed(inclination='all',aperture=-1)
#wav_rest_stellar =wav_rest_stellar* u.micron #wav is in micron                                                     
#nu_rest_stellar = constants.c.cgs/wav_rest_stellar.cgs
#lum_obs_stellar = dum_lum_stellar_sed
#lum_obs_stellar = lum_obs_stellar * u.erg/u.s
#nu_rest_stellar = constants.c.cgs/(wav_rest_stellar.to(u.cm))
#fnu_obs_stellar= lum_obs_stellar.to(u.Lsun)
#fnu_obs_stellar/= nu_rest_stellar.to(u.Hz)
#fnu_obs_stellar= fnu_obs_stellar.to(u.Lsun/u.Hz)


#idx = np.where((lam_rest_stellar.to(u.micron).value > 0.1) & (lam_rest_stellar.to(u.micron).value < 1.0))

#stellar_wav = lam_rest_stellar[idx]
#stellar_flux = lum_obs_stellar[idx]



#find location of vband in lam array
vband_idx = find_nearest(wav_rest_sed.to(u.angstrom).value,3000)
vband_extinction = np.zeros(fnu_obs_sed.shape[0])
for i in range(fnu_obs_sed.shape[0]): vband_extinction[i] = lum_obs_sed[i,vband_idx].cgs/lum_source_downsampled[vband_idx].cgs




#calculate extinction values
extinction = np.zeros([lum_obs_sed.shape[0],len(wav_rest_sed)])
for i in range(lum_obs_sed.shape[0]):
    extinction[i,:] = lum_obs_sed[i,:].cgs/lum_source_downsampled[:].cgs


#tau is calculated by: e^-tau = I/I_0
tau = -1.*np.log(extinction)
tau_v = -1.*np.log(vband_extinction)

inverse_wavelength = 1./wav_rest_sed.to(u.micron)


outfile = output_directory+'attenuation_curve.'+snap_str+'_galaxy'+galaxy_num+'.npz'

np.savez(outfile,wav_rest=wav_rest_sed.value,inverse_wavelength = inverse_wavelength.value,tau=tau,tau_v=tau_v,fnu_obs = fnu_obs_sed.value)






