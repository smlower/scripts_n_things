import fsps
import glob
import yt
import caesar
import numpy as np
import os, sys
import tqdm


def ssfr_relation(mstar):
    return (10**-13) * (10**mstar)



snapshot = sys.argv[1]

snap_dir = '/orange/narayanan/desika.narayanan/gizmo_runs/simba/m25n512/output/'
cfile = glob.glob(snap_dir+'/Groups/caesar_0'+"{:03d}".format(int(snapshot))+'_*.hdf5')
ids = np.load('/orange/narayanan/s.lower/simba/snap'+str(snapshot)+'_galaxy.npz', allow_pickle=True)
gal = ids['galid']
#gal_file = np.load('/orange/narayanan/s.lower/simba/_galaxy.npz')['galids']

print('loading yt snap')
ds = yt.load(snap_dir+'snapshot_'+snapshot+'.hdf5')


print('quick loading caesar')
obj = caesar.quick_load(cfile[0])
obj.yt_dataset = ds
dd = obj.yt_dataset.all_data()
print('loading particle data')
scalefactor = dd[("PartType4", "StellarFormationTime")]
star_masses = dd[("PartType4", "Masses")]
star_metal = dd[("PartType4", 'metallicity')]
print('now loading fsps')
fsps_ssp = fsps.StellarPopulation(sfh=0,
                zcontinuous=1,
                imf_type=2,
                zred=0.0, add_dust_emission=False)
solar_Z = 0.0196


# Compute the age of all the star particles from the provided scale factor at creation                               
formation_z = (1.0 / scalefactor) - 1.0

# Compute the age of all the star particles from the provided scale factor at creation                              
#yt_cosmo = yt.utilities.cosmology.Cosmology()
yt_cosmo = yt.utilities.cosmology.Cosmology(hubble_constant=0.68)
stellar_formation_age = yt_cosmo.t_from_z(formation_z).in_units("Gyr")
# Age of the universe right now                                                                                     
simtime = yt_cosmo.t_from_z(ds.current_redshift).in_units("Gyr")
stellar_ages = (simtime - stellar_formation_age).in_units("Gyr")


sfr_50 = []
for INDEX in tqdm.tqdm(gal):
    
    star_ages = stellar_ages[obj.galaxies[INDEX].slist]
    w50 = np.where(star_ages.in_units('Myr').value < 100)[0]
    #print('gal '+str(INDEX)+': ',len(w50))
    
    if len(w50) == 0:
        #sfr_50myr = 1.0e-10
        #sfr_50.append(sfr_50myr)
        #sfr_50.append(ssfr_relation(np.log10(obj.galaxies[INDEX].masses['stellar'].in_units('Msun').value)))
        sfr_50.append(0.0)
    else:
        initial_mass = 0.0
        for star in w50:
            
            current_mass = star_masses[star].in_units('Msun')
            fsps_ssp.params["logzsol"] = np.log10(star_metal[star] / solar_Z)
            mass_remaining = fsps_ssp.stellar_mass
            initial_mass += current_mass / np.interp(np.log10(stellar_ages[star]*1.e9),fsps_ssp.ssp_ages,mass_remaining)  
        
        sfr_50myr = initial_mass/100.e6
        sfr_50.append(sfr_50myr)




np.savez('/ufrc/narayanan/s.lower/caesar_sfr100_'+str(snapshot)+'_llim0.npz', sfr_50=sfr_50)
