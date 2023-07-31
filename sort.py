import numpy as np
import subprocess as sp

sort = 'concentration'
#sort = 'centroid'

clusters = []
sorted_clusters = []
vals = []

with open('clusters_to_plot.txt', 'r') as f:
    lines = f.readlines()

for line in lines:
    clusters.append(line.strip('\n'))

for cluster in clusters:
    if sort == 'concentration':
        with open('results/{}/results/concentration/conc.txt'.format(cluster), 'r') as f:
            vals.append(float(f.readlines()[0].strip('\n')))
    elif sort == 'centroid':
        with open('results/{}/results/centroid/w.txt'.format(cluster), 'r') as f:
            vals.append(float(f.readlines()[0].strip('\n')))

minimum = 10000000

inds = np.argsort(vals)

for ind in inds:
    sorted_clusters.append(clusters[ind])

command = 'montage '

clus = ''

for cluster in sorted_clusters:
    sp.call(['bash', 'crop.sh', cluster])
    command += 'results/{}/results/figures/labeled_no_border.png '.format(cluster)
    clus += 'results/{}/results/figures/cropped_smoothed.img -log -scalelims 0.000000005 0.00001 -scale log exp 100000 -cmap cool -cmap 2.5 0.5 -colorbar no -zoom TO FIT '.format(cluster)

command += '-geometry +0+0 -tile 10x10 sorted_{}.jpg'.format(sort)

#print(command)
print(clus)

