import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.family'] = 'serif'
rcParams['font.sans-serif'] = ['Times New Roman']
rcParams['font.size'] = 14
rcParams['figure.figsize'] = [8, 6]

fig, ax = plt.subplots()

with open('clusters_to_plot.txt') as f:
    clusters = f.readlines()

for cluster in clusters:
    cluster = cluster.strip('\n')

    mcmc_fit_file = f'/Users/laurelwhite/xray-analysis/results/{cluster}/results/ICM/ICM_best_fits.npz'
    mcmc_fit = np.load(mcmc_fit_file)

    mcmc_x = mcmc_fit['r']
    mcmc_y = mcmc_fit['te']

    cl_dir = f'/Users/laurelwhite/xray-analysis/results/{cluster}/results/cluster/'
    lsq_T_prof_file = cl_dir + 'T_prof_fit.npz'
    lsq_T_data = np.load(lsq_T_prof_file)

    lsq_x = lsq_T_data['fitx']
    lsq_y = lsq_T_data['fity']

    param_file = f'/Users/laurelwhite/xray-analysis/params/param_{cluster}.txt'

    with open(param_file) as f:
        lines = f.readlines()

    for line in lines:
        if line.split()[0] == 'R500':
            R500 = float(line.split()[-1])

    scaled_lsq_x = lsq_x / R500


#    best_fit_param_file = f'/Users/laurelwhite/xray-analysis/results/{cluster}/results/MCMC_kT/best_fit_params.npy'
#    best_fit_params = np.load(best_fit_param_file)
#    print(best_fit_params)

    ## Scaling relation from Vikhlinin '06
    ## M = M5 (T / 5 keV)** alpha
    M5 = 3.32 * 10**14
    alpha = 1.47

    all_param_file = '/Users/laurelwhite/xray-analysis/params/all_clusters.txt'

    with open(all_param_file) as f:
        lines = f.readlines()

    for line in lines:
        if line.split()[0] == f'{cluster}':
            M500 = float(line.split()[-1]) * 10**14

    kT500 = (M500 / M5) ** (1 / alpha) * 5

#    scaled_lsq_y = lsq_y / kT500
    scaled_lsq_y = lsq_y / np.max(lsq_y)

    ax.plot(scaled_lsq_x, scaled_lsq_y, color = 'steelblue', alpha = 0.2)
#    ax.plot(mcmc_x / R500, mcmc_y / np.max(mcmc_y), color = 'orange', alpha = 0.5)

ax.set_xlabel('Radius / R500')
#ax.set_ylabel('Temperature / kT500')
ax.set_ylabel('Temperature / max(temp)')
ax.set_xscale('log')
ax.set_xlim(0.1, 2)
ax.set_ylim(0, 1.25)
fig.savefig('profiles.png')

