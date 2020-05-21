from corner import quantile
import prospect.io.read_results as pread
import numpy as np
import sys
import glob, os
import tqdm
import ipdb

#prosp_dir = '/ufrc/narayanan/s.lower/simSEDs//simbam25n512_newfof/prod_runs/dirichlet/final/'
#pd_dir = '/ufrc/narayanan/s.lower/pd_runs/simba_m25n512/snap305/mist_pd/snap305/'

galaxy = int(sys.argv[1])
model = sys.argv[2]
prosp_dir = '/ufrc/narayanan/s.lower/simSEDs//simbam25n512_newfof/experiments/dust_screen/'+model+'/'
pd_dir = '//ufrc/narayanan/s.lower/pd_runs/simba_m25n512/snap305_dustscreen/snap305/'


sys.path.append(prosp_dir)
from run_prosp import build_model, build_sps
print('sps and model')
sps = build_sps()
mod = build_model()
thetas = mod.theta_labels()


def find_nearest(array,value):
    idx = (np.abs(np.array(array)-value)).argmin()
    return idx



def sample_posteriors(galaxy, num_samples=1000):

    print('now reading files')
    galaxy_num = "{:03d}".format(galaxy)
    infile = prosp_dir+'/galaxy'+galaxy_num+'.h5'
    globfiles = glob.glob(infile)

    try:
        for prosp_output in glob.glob(infile):
            print(prosp_output)
            res, obs, _ = pread.results_from(prosp_output)
    except:
        print('file not found')

    print('quantiles for all thetas')
    rand_draws = np.random.randint(0, np.shape(res['chain'])[0], 1000)

    #print(np.shape(res['chain']))

    theta_chains = []
    for theta in thetas:
        idx = thetas.index(theta)
        chain = [item[idx] for item in res['chain']]
        chain_draws = []
        print('taking random draws from chain for',theta)
        for draw in rand_draws:
            chain_draws.append(chain[draw])
        theta_chains.append(chain_draws)

    return obs, sps, mod, theta_chains


def sample_photometry(galaxy, num_samples=1000):
    print('generating photometry from theta chains')
    obs, sps, mod, theta_chains = sample_posteriors(galaxy, num_samples)
    phot_chain = []
    spec_chain = []

    theta_chains_rearranged = []
    for i in range(len(theta_chains[0])):
        theta_chains_rearranged.append([item[i] for item in theta_chains])

    for param_draw in tqdm.tqdm(theta_chains_rearranged):
        spec, phot,x = mod.mean_model(param_draw, obs, sps)
        phot_chain.append(phot)
        spec_chain.append(spec)


    phot_quan = []
    spec_quan = []
    for i in range(len(obs['maggies'])):
        phot_quan.append([quantile([item[i] for item in phot_chain], [.16, .5, .84])])
        spec_quan.append([quantile([item[i] for item in spec_chain], [.16, .5, .84])])

    return phot_quan, spec_quan, obs


def best_fit_phot(galaxy):
    

    print('now reading files')
    galaxy_num = "{:03d}".format(galaxy)
    infile = prosp_dir+'/galaxy'+str(galaxy)+'.h5'
    globfiles = glob.glob(infile)

    try:
        for prosp_output in glob.glob(infile):
            print(prosp_output)
            res, obs, _ = pread.results_from(prosp_output)
    except:
        print('file not found')
    #ipdb.set_trace()
    imax = np.argmax(res['lnprobability'])
    theta_max = res["chain"][imax, :]
    #obs_maggies = res['obs']['maggies']

    #print(theta_max)
    #print(sps)
    #print(obs_maggies)
    spec, phot, x = mod.mean_model(theta_max, obs,sps)

    return phot




galex = ['galex_FUV', 'galex_NUV']
hst_wfc3_uv  = ['wfc3_uvis_f275w', 'wfc3_uvis_f336w', 'wfc3_uvis_f475w','wfc3_uvis_f555w', 'wfc3_uvis_f606w', 'wfc3_uvis_f814w']
sdss = ['sdss_u0', 'sdss_g0', 'sdss_r0', 'sdss_i0', 'sdss_z0']
hst_wfc3_ir = ['wfc3_ir_f105w', 'wfc3_ir_f110w', 'wfc3_ir_f125w', 'wfc3_ir_f140w', 'wfc3_ir_f160w']
spitzer_mips = ['spitzer_mips_24']
herschel_pacs = ['herschel_pacs_70', 'herschel_pacs_100', 'herschel_pacs_160']
herschel_spire = ['herschel_spire_250', 'herschel_spire_350', 'herschel_spire_500']
filternames = (galex + hst_wfc3_uv + sdss + hst_wfc3_ir +  spitzer_mips + herschel_pacs + herschel_spire)



def true_phot(galaxy):    
    from hyperion.model import ModelOutput
    from astropy import units as u
    from astropy import constants

    from sedpy.observate import load_filters


    #print('galaxy: ',sys.argv[1])
    m = ModelOutput(pd_dir+"/snap305.galaxy"+str(galaxy)+".rtout.sed")
    wav,flux = m.get_sed(inclination=0,aperture=-1)
    wav  = np.asarray(wav)*u.micron #wav is in micron                                                                                                        
    wav = wav.to(u.AA)
    flux = np.asarray(flux)*u.erg/u.s
    dl = (10. * u.pc).to(u.cm)
    flux /= (4.*3.14*dl**2.)
    nu = constants.c.cgs/(wav.to(u.cm))
    nu = nu.to(u.Hz)
    flux /= nu
    flux = flux.to(u.Jy)
    maggies = flux / 3631.

    filters_unsorted = load_filters(filternames)
    waves_unsorted = [x.wave_mean for x in filters_unsorted]
    filters = [x for _,x in sorted(zip(waves_unsorted,filters_unsorted))]
    flx = []
    flxe = []

    for i in range(len(filters)):
        flux_range = []
        wav_range = []
        for j in filters[i].wavelength:
            flux_range.append(maggies[find_nearest(wav.value,j)].value)
            wav_range.append(wav[find_nearest(wav.value,j)].value)
        a = np.trapz(wav_range * filters[i].transmission* flux_range, wav_range, axis=-1)
        b = np.trapz(wav_range * filters[i].transmission, wav_range)
        flx.append(a/b)
        flxe.append(0.03* flx[i])

    flx = np.asarray(flx)
    flxe = np.asarray(flxe)
    flux_mag = flx
    unc_mag = flxe



    return flux_mag



if __name__ == '__main__':


    
    prosp_phot = best_fit_phot(galaxy)
    true_phot = true_phot(galaxy)
    
    residuals = (true_phot - prosp_phot) / true_phot
    print(residuals)

    np.savez(prosp_dir+'/phot_resid_galaxy'+str(galaxy)+'.npz', residuals=residuals, prosp_phot=prosp_phot, true_phot=true_phot, galaxy=galaxy)
