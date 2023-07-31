import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams

rcParams['font.family'] = 'serif'
rcParams['font.sans-serif'] = ['Times New Roman']
rcParams['font.size'] = 14
rcParams['figure.figsize'] = [8, 6]

## Get lit data
x = []
y = []

with open('centroid_shifts.txt', 'r') as f:
    lines = f.readlines()

for line in lines:
    x.append(float(line.split(',')[0]))
    y.append(float(line.split(',')[1]))

## Get my data
clusters = []
centroids = []

with open('clusters_to_plot.txt', 'r') as f:
    lines = f.readlines()

for line in lines:
    clusters.append(line.strip('\n'))

for cluster in clusters:
    with open('results/{}/results/centroid/w.txt'.format(cluster), 'r') as f:
        centroids.append(float(f.readlines()[0].strip('\n')))

norm = np.full_like(np.array(centroids), 1.0/len(clusters))

plt.plot(x, y, label='REXCESS', color='steelblue', linewidth=3)
plt.hist(centroids, bins=np.logspace(np.log10(0.001), np.log10(0.1), 10), histtype='step', weights=norm, label='CEREAL', color='lightsteelblue', linewidth=3)
#plt.xlim(0, 0.1)
plt.xlabel('Centroid shift')
plt.ylabel('Fraction')
plt.legend()
plt.xscale('log')
plt.savefig('centroid_hist.png')

