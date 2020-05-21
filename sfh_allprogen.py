#!/usr/bin/env python
# coding: utf-8

import yt
import caesar
import numpy as np
import fsps
import sys, os
from prospector_output_utilities import *

#old_caesar = pd.read_pickle('/Users/sidneylower/Documents/snap305_dirichlet/caesar_galaxy_properties_evo2000.pkl')
caesar_snap = '/orange/narayanan/desika.narayanan/gizmo_runs/simba/m25n512/output/Groups/caesar_0305_z0.000.hdf5'
snapshot = '/orange/narayanan/desika.narayanan/gizmo_runs/simba/m25n512/output/snapshot_305.hdf5'

gal = int(sys.argv[1])

obj = caesar.quick_load(caesar_snap)
solar_Z = 0.0196

ds = yt.load(snapshot)
obj.yt_dataset = ds
dd = obj.yt_dataset.all_data()

scalefactor = dd[("PartType4", "StellarFormationTime")]
star_masses = dd[("PartType4", "Masses")]
star_Z = dd[("PartType4", "metallicity")]
formation_z = (1.0 / scalefactor) - 1.0
yt_cosmo = yt.utilities.cosmology.Cosmology(hubble_constant=0.68)
stellar_formation_age = yt_cosmo.t_from_z(formation_z).in_units("Gyr")
simtime = yt_cosmo.t_from_z(ds.current_redshift).in_units("Gyr")
stellar_ages = (simtime - stellar_formation_age).in_units("Gyr")

#fsps_years = np.linspace(0.01, 14, 200)
#fsps_years = np.arange(0.01, simtime.in_units("Gyr").value, 0.1)
#sfr_bins = fsps_years
#sfr_times = [(j+i)/2 for i, j in zip(sfr_bins[:-1], sfr_bins[1:])] #lookback time

fsps_ssp = fsps.StellarPopulation(sfh=0,
                zcontinuous=1,
                imf_type=2,
                zred=0.0, add_dust_emission=False)

fsps_years = (10**fsps_ssp.ssp_ages[:-3]) / 1.0e9


sfr_list = []
Z_list = []
time = []
for galaxy in [gal]:
    print('On galaxy: ',galaxy)
    sfr_onegal = []
    #gal_list.append(galaxy)
    stellar_mass = star_masses[obj.galaxies[galaxy].slist].in_units('Msun')
    stellar_metallicity = star_Z[obj.galaxies[galaxy].slist]
    star_ages = stellar_ages[obj.galaxies[galaxy].slist]
    smass = stellar_mass
    sage = star_ages

    for i in range(len(fsps_years)-1):
        time.append(fsps_years[-1] - fsps_years[i])
        age_bin = np.where((sage >= (fsps_years[i])) & (sage < (fsps_years[i+1])))[0]
        if len(age_bin) == 0:
            #print(sfr_bins[i])
            sfr_current_bin = 1.0e-10
            sfr_onegal.append(sfr_current_bin)
            if fsps_years[i] < 1.:
                Z_list.append(-6)
            else:
                Z_list.append(Z_list[i-1])
        else:
            current_mass = smass[age_bin]
            initial_mass = np.array([])
            Z = np.array([])
            print('on age bin: ',age_bin)
            for star in age_bin:
                fsps_ssp.params["logzsol"] = np.log10(stellar_metallicity[star] / solar_Z)
                Z = np.append(Z, np.log10(stellar_metallicity[star] / solar_Z))
                mass_remaining = fsps_ssp.stellar_mass
                initial_mass = np.append(initial_mass, smass[star] / np.interp(np.log10(sage[star]*1.e9),fsps_ssp.ssp_ages,mass_remaining))  

            Z_list.append(np.sum(Z * initial_mass) / np.sum(initial_mass))
            sfr_current_bin = sum(initial_mass)/np.abs(fsps_years[i+1]*1e9 - fsps_years[i]*1e9)
            sfr_onegal.append(sfr_current_bin)

    sfr_list = sfr_onegal


np.savez('/orange/narayanan/s.lower/simba/sfh_Z_galaxy'+str(gal)+'_highres.npz', Galaxy=gal, SFR=sfr_list, time=time, Z=Z_list)







