import matplotlib.pyplot as plt
import numpy as np






#stellar mass stuff from simSEDs_output script

dist = []
for i in range(len(int_smass)):
    dist.append(np.log10(int_smass[i] - der_smass[i]))


plt.figure(figsize=(16, 12))
n, bins, _ = plt.hist(dist, bins=50, edgecolor='darkorange', color='None')
plt.xlabel('Log(M$_{intrinsic}$ - M$_{derived}$', fontsize=23)
plt.ylabel('# Galaxies', fontsize=23)
plt.savefig(directory+'stellarmass_offset.png', dpi=300)
plt.close()
