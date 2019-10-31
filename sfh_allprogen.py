#!/usr/bin/env python
# coding: utf-8

# In[1]:


import yt
import caesar
import matplotlib.pyplot as plt
import numpy as np
import fsps
import pandas as pd
from astropy.cosmology import FlatLambdaCDM
from readgadget.modules import header as HEAD
from readgadget.modules.common import *
from readgadget.modules.names import *
from readgadget.modules.hdf5 import *
import sys, os
import h5py


# In[24]:


def fsps_isochrone_data(mass, metallicity, age):
    
    metal_dict = np.array([0.000041,0.00013,0.000231,0.000411,0.000731,0.0013,
                  0.002312,0.004111,0.00731,0.013,0.023118,0.04111])
    
    iso_Z = metal_dict[np.abs(metallicity - metal_dict).argmin()]
    
    isochrone_file = '/Users/sidneylower/fsps/MIST_interp_isochrones/MIST_z'+str(iso_Z)+'.dat'
    
    isochrone = pd.read_csv(isochrone_file, delim_whitespace=True, skiprows=10, comment='#', usecols=['log10_isochrone_age_yr', 'initial_mass', 'star_mass'])
    
    #initial_mass = isochrone['initial_mass'][np.where(np.abs(np.log10(age*1.0e9) - isochrone['log10_isochrone_age_yr']).min() & np.abs(mass - isochrone['star_mass']).min())]
    
    
    age_idx = np.where(isochrone['log10_isochrone_age_yr'] == isochrone['log10_isochrone_age_yr'][np.abs(np.log10(age*1.0e9) - isochrone['log10_isochrone_age_yr']).idxmin()])[0]
    print(age_idx)
    mass_idx = np.where(isochrone['star_mass'] == isochrone['star_mass'][np.abs(mass - isochrone['star_mass']).idxmin()])[0]
    print(mass_idx)
    full_idx = list(set(age_idx) & set(mass_idx))
    print(full_idx)
    initial_mass = isochrone['initial_mass'][full_idx]
    
    return initial_mass


# In[2]:


#snapshot = '/Volumes/ufrc/narayanan/s.lower/pd_runs/simba_m25n512/snap305_boxtest/filtered_snaps/snap305_galaxy003_filtered.hdf5'
caesar_snap = '/Volumes/orange/narayanan/desika.narayanan/gizmo_runs/simba/m25n512/output/Groups/caesar_0305_z0.000.hdf5'


# In[3]:


obj = caesar.quick_load(caesar_snap)
sim = obj
        


# In[12]:


def readsnap(snap,data,ptype,**kwargs):
    h   = HEAD.Header(snap,0,kwargs)
    d,p = pollOptions(h,kwargs,data,ptype)
    h.reading = d

    #print 'reading %s.%s' % (snap,h.extension)

    f = h.f
    initUnits(h)
    
    #print 'reading %d files...' % (h.nfiles)

    for i in range(0,h.nfiles):
        if i > 0:
            h = HEAD.Header(snap,i,kwargs)
            f = h.f
            h.reading = d
            initUnits(h)

        if h.npartThisFile[p] == 0:
            if h.nfiles > 1:
                continue
            print('no %s particles present!' % pNames[p])
            sys.exit()

        if h.fileType == 'hdf5':
            arr = hdf5.hdf5_read(f,h,p)
        elif h.fileType == 'tipsy':
            arr = tipsy.tipsy_read(f,h,p)
        elif h.fileType == 'gadget1':
            arr = gadget1.gadget_read(f,h,p,d)
        elif h.fileType == 'gadget2':
            arr = gadget2.gadget_type2_read(f,h,p)

        f.close()

        ## return nth value
        if h.nth > 1:
            arr = arr[0::h.nth]

        ## put arrays together
        if i > 0:
            if len(arr) > 0:
                return_arr = np.concatenate((return_arr,arr))
        else:
            return_arr = arr
            gadgetPrinter(h,d,p)
            if h.nth > 1 and not h.suppress:
                print('selecting every %d particles' % h.nth)

        ## if requesting a single file, break out of loop
        if h.singleFile:
            break

    if h.double and h.reading != 'pid' and h.reading != 'ParticleIDs':
        return_arr = return_arr.astype(np.float64)
                                    
    return return_arr


# In[13]:


def readhead(snap,data,**kwargs):
    h = HEAD.Header(snap,0,kwargs)
    h.f.close()
    pollHeaderOptions(h,data)
    if data == 'header' or data == 'Header':
        return h.header_vals
    return h.header_vals[headerTypes[data]]


# In[14]:


def tage(cosmo,thubble,a):
    return thubble-cosmo.age(1./a-1).value

def get_tage(cosmo,thubble,tform):
    sa = np.asarray([tage(cosmo,thubble,i) for i in tform])
    return sa

def t_elapsed():
    return np.round(time.time()-TINIT,2)


# In[4]:


params = {'fsps_sfh' : 0,
          'fsps_zcontinuous' : 1,
          'fsps_imf_type' : 2,
          'fsps_zred' : 0.0,
          'star_metal_index' : 0
        }

cosmo = FlatLambdaCDM(H0=100*sim.simulation.hubble_constant, Om0=sim.simulation.omega_matter, Ob0=sim.simulation.omega_baryon,Tcmb0=2.73)
solar_Z = 0.0196


# In[5]:


fsps_ssp = fsps.StellarPopulation(sfh=params['fsps_sfh'],
                zcontinuous=params['fsps_zcontinuous'],
                imf_type=params['fsps_imf_type'],
                zred=params['fsps_zred'], add_dust_emission=False)


# In[8]:


old_caesar = pd.read_pickle('/Users/sidneylower/Documents/snap305_dirichlet/caesar_galaxy_properties_evo2000.pkl')
ultra = pd.read_pickle('/Volumes/ufrc/narayanan/s.lower/simSEDs/simbam25n512_newfof/snap305_ultra/reparam_noll/prospectorsnap305_properties_Dirichlet_1570543762.1712768.pkl')


# In[154]:


galaxy = 56


# In[155]:


snapshot = '/Volumes/ufrc/narayanan/s.lower/pd_runs/simba_m25n512/snap305_boxtest/filtered_snaps/snap305_galaxy{:03d}'.format(galaxy)+'_filtered'


# In[116]:


h = sim.simulation.hubble_constant
sm = readsnap(snapshot,'mass','star',units=1,suppress=1)#/h
tform = readsnap(snapshot,'age','star',units=1,suppress=1)  # expansion factor of formation
smetal = readsnap(snapshot,'metallicity','star',units=1,suppress=1)
thubble = cosmo.age(sim.simulation.redshift).value
sa = get_tage(cosmo,thubble,tform)


# In[117]:


#fsps_years = np.arange(5.0,10.3, 0.05)
#sfr_bins = (10**fsps_years)/1.e9

sfr_bins = np.linspace(1e-4, 13.8, 130)


# In[156]:


sfr = []
for i in range(len(sfr_bins)-1):
    #stellar_mass = part_mass
    #stellar_metallicity = star_metallicity[obj.galaxies[galaxy].slist]
    #star_ages = stellar_ages[obj.galaxies[galaxy].slist]
    #plist = np.array(obj.galaxies[galaxy].slist)  
    smass = sm
    sage = sa
    age_bin = np.where((sage > (sfr_bins[i])) & (sage < (sfr_bins[i+1])))[0]
    #print(age_bin)
    if len(age_bin) == 0:
        sfr_current_bin = 1.0e-10
        sfr.append(sfr_current_bin)
    else:
        current_mass = smass[age_bin]
        #current_Z = stellar_metallicity[age_bin]
        initial_mass = 0.0
        for star in age_bin:
            fsps_ssp.params["logzsol"] = np.log10(smetal[star] / solar_Z)
            mass_remaining = fsps_ssp.stellar_mass
            initial_mass += smass[star] / np.interp(np.log10(sage[star]*1.e9),fsps_ssp.ssp_ages,mass_remaining)  
            #initial_mass += fsps_isochrone_data(current_mass[star], current_Z[star].value, star_ages[star].in_units('Gyr').value)
        sfr_current_bin = initial_mass/np.abs(sfr_bins[i+1]*1e9 - sfr_bins[i]*1e9)
        sfr.append(sfr_current_bin)


# In[157]:


c_timelist = list(dict.fromkeys(old_caesar.index.get_level_values('Time [Gyr]')))
c_sfrlist = old_caesar.loc[[galaxy]]['SFR']


# In[158]:


ultra_timelist = list(dict.fromkeys(ultra.index.get_level_values('Time [Gyr]')))
ultra_sfrlist = list(dict.fromkeys(ultra.loc[['{0:03}'.format(galaxy)]]['SFR']))
ultra_sfr16 = list(dict.fromkeys(ultra.loc[['{0:03}'.format(galaxy)]]['SFR_16']))
ultra_sfr84 = list(dict.fromkeys(ultra.loc[['{0:03}'.format(galaxy)]]['SFR_84']))


# In[121]:


simtime = obj.simulation.time.in_units('Gyr').value


# In[161]:


plt.figure(figsize=(10, 8))
plt.rc('axes', linewidth=2)
plt.tick_params(axis='both', which='major', labelsize=15)
plt.tick_params(axis='both', which='minor')
plt.plot(simtime - sfr_bins[:len(sfr)], sfr, lw=1, label='New', zorder=1)
plt.plot(c_timelist, c_sfrlist, lw=1, label='Old', zorder=0)


for i in range(12):
    plt.plot([ultra_timelist[i], ultra_timelist[i+1]], [ultra_sfrlist[i], ultra_sfrlist[i]], color='lightgreen', lw=4, alpha=0.7, zorder=5)
    plt.fill_between([ultra_timelist[i], ultra_timelist[i+1]], y1=ultra_sfr16[i], y2=ultra_sfr84[i], color='lightgreen', alpha=0.5, zorder=5)

plt.plot([ultra_timelist[0],ultra_timelist[1]], [ultra_sfrlist[0], ultra_sfrlist[0]], color='lightgreen', lw=3, label='Ultra', zorder=5)


#plt.scatter(sfr_bins[:999], sfr)
plt.yscale('log')
plt.ylim([1e-2, 1e3])
plt.xlim([0, 14])
plt.ylabel('SFR', fontsize=23)
plt.xlabel('Time [Gyr]', fontsize=23)
plt.title('SFR: Galaxy '+str(galaxy), fontsize=23)
plt.legend(loc='best', fontsize=15)
#plt.savefig('/Users/sidneylower/Documents/prosp_final_plots/sfr_comp'+str(galaxy)+'.png', dpi=300, bbox_inches='tight')


# In[273]:


fsps_ssp.params["logzsol"] = np.log10(0.00041 / solar_Z)
mass_remaining = fsps_ssp.stellar_mass


# In[270]:


mass_remaining


# In[274]:


mass_remaining


# In[37]:


np.interp(np.log10(sage[star]*1.e9),fsps_ssp.ssp_ages,mass_remaining)


# In[114]:


np.log10(smetal / solar_Z)


# In[131]:


mstar0 = np.asarray([i.masses['stellar'] for i in sim.galaxies if i.central <= 1])


# In[134]:


len(mstar0)


# In[135]:


len(sim.galaxies.masses['stellar'])


# In[152]:


list(i.GroupID for i in sim.galaxies if i.central == 0)


# In[153]:


sim.galaxies[3].central


# In[149]:





# In[162]:


np.log10(sage[0]*1.e9),fsps_ssp.ssp_ages


# In[ ]:




