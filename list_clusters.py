import os

for cluster in os.listdir('results/'):
    if os.path.exists(f'results/{cluster}/results/centroid/w.txt') and os.path.exists(f'results/{cluster}/results/concentration/conc.txt'):
        print(cluster)

