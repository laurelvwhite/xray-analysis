import matplotlib.pyplot as plt
import numpy as np

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

with open('aas_hist.txt', 'r') as f:
    lines = f.readlines()

for line in lines:
    clusters.append(line.strip('\n'))

for cluster in clusters:
    with open('results/{}/results/centroid/w.txt'.format(cluster), 'r') as f:
        centroids.append(float(f.readlines()[0].strip('\n')))

norm = np.full_like(np.array(centroids), 1.0/len(clusters))

plt.rcParams.update({'font.size': 12})

plt.plot(x, y, label='REXCESS', color='purple', linewidth=3)
plt.hist(centroids, bins=np.logspace(np.log10(0.001), np.log10(0.1), 10), histtype='step', weights=norm, label='CEREAL', color='orange', linewidth=3)
#plt.xlim(0, 0.1)
plt.xlabel('Centroid shift')
plt.ylabel('Fraction')
plt.legend()
plt.xscale('log')
plt.savefig('centroid_hist.png')

