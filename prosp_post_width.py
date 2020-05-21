import prospect.io.read_results as pread
import numpy as np
from glob import glob
import tqdm
import sys
from corner import quantile


prosp_dir = '/ufrc/narayanan/s.lower/simSEDs/simbam25n512_newfof/prod_runs/dirichlet/old/'

sys.path.append(prosp_dir)
from run_prosp import build_model, build_sps


sps = build_sps()
mod = build_model()

files_ = glob(prosp_dir+'/galaxy*.h5')


mass_sigma = []
mass_chains = []
for galaxy in tqdm.tqdm(files_):
    res, obs, _ = pread.results_from(galaxy)
    thetas = mod.theta_labels()

    mass = [item[thetas.index('massmet_1')] for item in res['chain']]
    mass_chains.append(mass)
    mass_std = np.std(mass)
    mass_sigma.append(mass_std)
    
    

np.savez(prosp_dir+'mass_posterior_widths.npz', sigma=mass_sigma, posterior=mass_chains)
