import yt
import numpy as np
import pandas as pd
import fsps
from hyperion.model import ModelOutput
from sedpy.observate import load_filters
import sys, os, glob
import astropy.constants as constants
import astropy.units as u



galex = ['galex_FUV', 'galex_NUV']
hst_wfc3_uv  = ['wfc3_uvis_f275w', 'wfc3_uvis_f336w', 'wfc3_uvis_f475w','wfc3_uvis_f555w', 'wfc3_uvis_f606w', 'wfc3_uvis_f814w']
sdss = ['sdss_i0']
hst_wfc3_ir = ['wfc3_ir_f105w', 'wfc3_ir_f125w', 'wfc3_ir_f140w', 'wfc3_ir_f160w']
spitzer_irac = ['spitzer_irac_ch1']
spitzer_mips = ['spitzer_mips_24']
herschel_pacs = ['herschel_pacs_70', 'herschel_pacs_100', 'herschel_pacs_160']
herschel_spire = ['herschel_spire_250', 'herschel_spire_350', 'herschel_spire_500']
jwst_miri = ['jwst_f560w', 'jwst_f770w', 'jwst_f1000w', 'jwst_f1130w', 'jwst_f1280w', 'jwst_f1500w', 'jwst_f1800w']
jwst_nircam = ['jwst_f070w', 'jwst_f090w', 'jwst_f115w', 'jwst_f150w', 'jwst_f200w', 'jwst_f277w']
filternames = (galex + hst_wfc3_uv +  hst_wfc3_ir + jwst_miri + jwst_nircam + spitzer_irac + spitzer_mips + herschel_pacs + herschel_spire)
filters_unsorted = load_filters(filternames)
waves_unsorted = [x.wave_mean for x in filters_unsorted]
filters = [x for _,x in sorted(zip(waves_unsorted,filters_unsorted))]

print('Loading SP')
fsps_ssp = fsps.StellarPopulation(sfh=0,
                zcontinuous=1,
                imf_type=2,
                zred=0.0, add_dust_emission=False)
solar_Z = 0.0196


snaps_dir = '/orange/s.lower/eagle/filtered_snapshots/snap028/'
pd_dir = '/orange/s.lower/eagle/pd_runs/snap28/'
snaps = glob.glob(snaps_dir+'/galaxy*.hdf5')

for galaxy in range(len(snaps)):
    
    try:
        m = ModelOutput(pd_dir+'/snap28.galaxy'+str(galaxy)+'.rtout.sed')
        ds = yt.load(snaps_dir+'/galaxy_'+str(galaxy)+'.hdf5')
    except:
        continue
    
    
    #get mock photometry
    wav, flux = m.get_sed(inclination='all', aperture=-1)
    wav  = np.asarray(wav)*u.micron #wav is in micron                                                                                                        
    wav = wav.to(u.AA)
    flux = np.asarray(flux)*u.erg/u.s
    dl = (10. * u.pc).to(u.cm)
    flux /= (4.*3.14*dl**2.)
    nu = constants.c.cgs/(wav.to(u.cm))
    nu = nu.to(u.Hz)
    flux /= nu
    flux = flux.to(u.Jy)
    for i in range(len(filters)):
        flux_range = []
        wav_range = []
        for j in filters[i].wavelength:
            flux_range.append(flux[find_nearest(wav.value,j)].value)
            wav_range.append(wav[find_nearest(wav.value,j)].value)
        a = np.trapz(wav_range * filters[i].transmission* flux_range, wav_range, axis=-1)
        b = np.trapz(wav_range * filters[i].transmission, wav_range)
        flx.append(a/b)
        flxe.append(0.03* flx[i])

    flx = np.asarray(flx)
    flxe = np.asarray(flxe)


    #get properties: M*, SFR, Z, Mdust
    ad = ds.all_data()
    star_masses = ad[('PartType4', 'Masses')]
    scalefactor = dd[("PartType4", "StellarFormationTime")]
    formation_z = (1.0 / scalefactor) - 1.0
    stellar_formation_age = ds.cosmology.t_from_z(formation_z).in_units("Gyr")
    simtime = ds.cosmology.t_from_z(ds.current_redshift).in_units("Gyr")
    stellar_ages = (simtime - stellar_formation_age).in_units("Gyr")


    
