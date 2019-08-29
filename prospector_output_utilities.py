
# coding: utf-8

'''
Ultimate Prospector Output Script. 
This baby has everything: parameter estimates, mass derivations, and SFH generation.
'''


import astropy.constants as const
from astropy.cosmology import WMAP9
from astropy import units as u
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
#import track_progens


def B_nu(nu,T):
    c, k_B, h = const.c.value, const.k_B.value, const.h.value
    return 2.0 * h * nu**3. * c**-2. / (np.exp(h*nu/(k_B * T)) - 1.0)


# In[2]:


def get_dustmass(Snu,Td,z,nu=345.0):
    """
    Quick dust mass calculator, using Greve+12 formula (Eq. 2)
    
    Parameters:
    Snu: float or array of floats or Quantity's
        Flux density at observed frequency nu, assumed mJy if no unit
    z: float
        Redshift of object
    nu: float or Quantity; default 345GHz
        Observed frequency, assumed GHz if no unit
    
    Returns:
    Md: float
        Apparent dust mass in solar masses
    """
    
    cosmo=WMAP9
    
    # Do some unit handling
    if hasattr(Snu,'unit'):
        Snu = Snu.to('mJy').decompose().value
    else: Snu = u.Quantity(Snu,'mJy').decompose().value
    if hasattr(nu,'unit'):
        nu = (nu.to('GHz')).value

    nu_rest = nu * (1. + z)
    beta = 2.0         # Dust emissivity index
    # Dust opacity; see Hildebrand83,Kruegel&Siebenmorgen94, in m**2 / kg
    # Below value used with beta = 2 in Greve+12
    kappa_nu = 0.045 * (nu_rest/250.)**beta
    # Updated emissivity, reproduces MW dust, see Bianchi+13, beta = 1.5
    #kappa_nu = 0.34 * (nu_rest/1200.)**beta
    # Scoville+14 version, beta = 1.8, from Planck maps, d_GDR = 100
    #kappa_nu = 4.84e-2 * (nu_rest/345.)**beta
    # From Draine&Li07 models / Li&Draine01; beta = 1.7
    #kappa_nu = 0.0431 * (nu_rest/352.)**beta
    #kappa_nu = 3.13 * (nu_rest/3000.)**beta
    
    DL = u.Quantity(cosmo.luminosity_distance(z),'Mpc').to('m').value
    
    denom = B_nu(1e9*nu_rest,Td) - B_nu(1e9*nu_rest,cosmo.Tcmb(z).value)
    Md = u.Quantity(DL**2. * Snu / ((1.+z) * kappa_nu * denom), 'kg')
    Md = Md.to(u.M_sun).value
    
    return Md


# In[3]:


def get_td(flux_array, wavelengths_array, z):
    
    """
    Assuming angstrom for wavelengths, in observer's frame. Shift to restframe with given redshifts. 
    Does not matter what units flux is in. 
    
    """
    
    
    #print(wavelengths_array)
    wav_rest = wavelengths_array / (1.0 + z)
    #print(wav_rest)
    #print(np.where(wav_rest <= 5.0e5)[0])
    dust_peak = np.where(wav_rest >= 5.0e5)[0]
    peak_wavelength = wav_rest[dust_peak][flux_array[dust_peak].argmax()]
    peak_wavelength *= 1.0e-8 #to cm for cgs version of Wien's constant
    tdust = const.b_wien.cgs.value/peak_wavelength
    
    return tdust


# In[4]:


def find_nearest(array,value):
    idx = (np.abs(array - value)).argmin()
    return idx


# In[5]:


def get_lir_lum(flux, wavelengths, z):
    
    """
    Luminosity from spectrum flux, assuming Jy for flux, angstrom for wavelengths.
    
    TO DO: write so that it produces LIR and FUV luminosity --> future thing for when pd dumps that?
    
    """
    
    #wavelengths_to_m = wavelengths / 1e-10 #meter
    
    cosmo = WMAP9
    
    dl = cosmo.luminosity_distance(z).to(u.m)
    
    
    lum = ((4*np.pi*(dl)**2) * (flux*u.Jy) * (const.c.to(u.m/u.s))) / (wavelengths*u.angstrom).to(u.m)
    
    lum_W = lum.to(u.W)
    
    model_spec_Jy = model_spec[i]  * 3631.
    lum = get_luminosity(model_spec_Jy, model_wave[i], zs[i])
    
    model_wave_m = model_wave[i] / 1.0e10
    
    iw1 = find_nearest(model_wave_m,8.0e-6)
    iw2 = find_nearest(model_wave_m,1000.0e-6)
    
    lir = np.trapz(lum[iw1:iw2]/model_wave_m[iw1:iw2],model_wave_m[iw1:iw2])
    print(lir)
    l_ir[i] = lir / 3.828e26

    return lum_W.value
    


# In[8]:


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


# In[ ]:


#plot SED

def plot_SED(phot_wav, spec_wav, photometry, mod_photometry, spectrum, phot_err, directory, gal_num):
    plt.figure(figsize=(12, 10))

    plt.tick_params(axis='both', which='major', labelsize=13)
    plt.tick_params(axis='both', which='minor', labelsize=10)


    plt.loglog(spec_wav, spectrum, label='Best Fit Model',
           lw=1.5, color='black')
    plt.errorbar(phot_wav, mod_photometry, label='Model photometry',
             marker='s', markersize=10, ls='', lw=3, 
             markerfacecolor='none', markeredgecolor='orange', 
             markeredgewidth=3)
    plt.errorbar(phot_wav, photometry, yerr=phot_err, 
             label='Observed photometry', ecolor='blue', 
             marker='o', markersize=10, ls='', lw=3, 
             markerfacecolor='none', markeredgecolor='blue', 
             markeredgewidth=3)
    plt.xlabel('Wavelength [$\AA$]', fontsize=25)
    plt.ylabel('Flux [maggies]', fontsize=25)
    plt.legend(loc='best', fontsize=15, frameon=False)
    plt.ylim([np.min(photometry) - 0.5*np.min(photometry), np.max(photometry) + 0.5*np.max(photometry)])
    plt.title('Snap305 Galaxy '+str(gal_num), fontsize=20)
    plt.savefig(directory+'snap305.galaxy_'+str(gal_num)+'_SED.png', dpi=300)
    plt.close()
    return


# In[ ]:


def physical_properties_plot(i_smass, d_smass, directory, para_data=None, add_smass=True, add_lir=False, add_dmass=False, add_sfr=False):
    #stellar mass, etc. plot
    p_data = pd.read_csv(para_data, delimiter=',')
    if add_smass is True:
        figure = plt.figure(figsize=(10, 8))
        plt.tick_params(axis='both', which='major', labelsize=14)
        plt.tick_params(axis='both', which='minor', labelsize=10)
        plt.hexbin(i_smass, d_smass, xscale='log', yscale='log', bins='log', extent=[6, 13, 6, 13])
        plt.plot([1e6, 1e13],[1e6, 1e13], color='gray', lw=0.9, ls='-.')
        plt.scatter(p_data['int_smass'], p_data['der_smass'], marker='+', c='gray', s=30)
        #plt.yscale('log')
        #plt.xscale('log')
        plt.ylim([1e6, 1e13])
        plt.xlim([1e6, 1e13])
        plt.ylabel('Derived M$_*$ [M$_{\odot}$]', fontsize=18)
        plt.xlabel('Intrinsic M$_*$ [M$_{\odot}$]', fontsize=18)
        plt.savefig(directory+'stellarMass.png', dpi=300)
        plt.close()

        ratio = np.array(i_smass) / np.array(d_smass)
        figure2 = plt.figure(figsize=(10, 8))
        plt.tick_params(axis='both', which='major', labelsize=14)
        plt.tick_params(axis='both', which='minor', labelsize=10)
        plt.scatter(i_smass, ratio, marker='o', c='orange', s=100)
        
        plt.plot([1e6, 1e13],[1., 1.], color='black', lw=0.6, ls='-.')
        
        plt.xscale('log')
        plt.ylabel(r'$\frac{M_{*,derived}}{M_{*,intrinsic}}$', fontsize=18)
        plt.xlabel('Intrinsic M$_*$ [M$_{\odot}$]', fontsize=18)
        plt.savefig(directory+'stellarMass_ratio.png', dpi=300)
        plt.close()



    if add_lir is True:
        print('currently under construction')

    if add_dmass is True:
        print('currently under construction')

    if add_sfr is True:
        print('currently under construction')

    return


def scatter_distribution_plot(cont_smasses, dirilect_smasses, para_smasses, paraburst_smasses):
    dist_cont = []
    dist_dirilect = []
    dist_para = []
    dist_paraburst = []

    for i in range(len(cont_data['int_smass'])):
        dist_cont.append(np.log10(cont_data['int_smass'][i]) - np.log10(cont_data['der_smass'][i]))
    
    for i in range(len(dirilect_data['int_smass'])):
        dist_dirilect.append(np.log10(dirilect_data['int_smass'][i]) - np.log10(dirilect_data['der_smass'][i]))

    for i in range(len(para_data['int_smass'])):
        dist_para.append(np.log10(para_data['int_smass'][i]) - np.log10(para_data['der_smass'][i]))

    for i in range(len(paraburst_data['int_smass'])):
        dist_paraburst.append(np.log10(paraburst_data['int_smass'][i]) - np.log10(paraburst_data['der_smass'][i]))

    plt.figure(figsize=(10, 8))
    n, bins, _ = plt.hist(dist_cont, bins=50, edgecolor='darkorange', color='None', label='Continuity SFH')
    n, bins, _ = plt.hist(dist_dirilect, bins=50, edgecolor='blue', color='None', label='Dirichlet SFH')
    n, bins, _ = plt.hist(dist_para, bins=50, edgecolor='green', color='None', label='Delayed Tau SFH')
    n, bins, _ = plt.hist(dist_paraburst, bins=50, edgecolor='purple', color='None', label='Delayed Tau + Burst SFH')
    plt.xlabel('Log(M$_{intrinsic}$ - M$_{derived}$', fontsize=23)
    plt.ylabel('# Galaxies', fontsize=23)
    plt.savefig(directory+'stellarmass_offset.png', dpi=300)
    plt.close()


    return








def SFH_plot(par_time, par_sfr, non_timebins, non_sfrbins, snaptime, intrinsic_sfr, save_directory, galaxy, SFHmodel, ifburst=False):
    plt.figure(figsize=(10, 8))
    plt.tick_params(axis='both', which='major', labelsize=13)
    plt.tick_params(axis='both', which='minor', labelsize=10)

    #parametric                                                                                                                            
    plt.plot(time, sfr_par, color='blue', lw=3, label='Derived SFH_par'+ifburst)
    #non parametric                                                                                                                        
    for i in range(len(non_sfrbins)):
        plt.plot([non_timebins[i][0], non_timebins[i][1]], [non_sfrbins[i], non_sfrbins[i]], color='orange', lw=3, label='Derived SFH_'+SFHmodel)
    #intrinsic                                                                                                                         
    plt.plot(snaptime, intrinsic_sfr, color='black', lw=0.5, label='Intrinsic SFH')
    plt.ylabel('SFR [M$_{\odot}$ / yr]', fontsize=25)
    plt.xlabel('Lookback Time [Gyr]', fontsize=25)
    plt.xlim([0, 14])
    plt.legend()
    plt.savefig(save_directory+'SFH_galaxy'+str(galaxy)+'.png', dpi=300)
    
    return 


def nonpara_massweighted_age(sfh, time):
    top = 0.0
    bottom = 0.0
    for bin_ in range(len(sfh)):
        top += sfh[bin_] * time[bin_]
        bottom += sfh[bin_]
    return top / bottom



def true_massweighted_age(sfh, time):
    from scipy.integrate import simps
    sfh = np.asarray(sfh)
    sfh_time = []
    for i in range(len(time)):
        sfh_time.append(sfh[i]*time[i])
    
    top = simps(sfh_time, dx=0.01)
    bottom = simps(sfh, dx=0.01)
    
    return top/bottom



def tau_massweighted_age(sfh, time):
    from scipy.integrate import simps
    sfh = np.asarray(sfh)
    sfh_time = []
    for i in range(len(time)):
        sfh_time.append(sfh[i]*time[i])

    top = simps(sfh_time, dx=0.01)
    bottom = simps(sfh, dx=0.01)

    return top/bottom
