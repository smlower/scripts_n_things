import numpy as np
import prospect.io.read_results as pread
import sys, os
import glob
from prospector_output_utilities import *


prosp_dir = sys.argv[1]

galaxies=sys.argv[2]

for galaxy_number in galaxies:
    print(galaxy_number)
    catch = 0
    res, obs, mod = 0, 0, 0
    galaxy_num = "{:03d}".format(int(galaxy_number))
    infile = prosp_dir+'/snap305.galaxy'+galaxy_num+'_*_mcmc.h5'
    for prosp_output in glob.glob(infile):
        res, obs, mod = pread.results_from(prosp_output)
        print(prosp_output)
    
    sps = pread.get_sps(res)
    print('got sps')
    test = pread.get_model(res)
    print(test)
    thetas, theta_best = get_best(res)
    print('got thetas')
    print(galaxy_num)
    mean_mod = mod.mean_model(theta_best, obs, sps)

    mass_frac = mean_mod[-1]
