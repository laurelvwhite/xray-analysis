import matplotlib.pyplot as plt
import numpy as np

## Get lit data
x = []
y = []

with open('concentrations.txt', 'r') as f:
    lines = f.readlines()

for line in lines:
    x.append(float(line.split(',')[0]))
    y.append(round(float(line.split(',')[1])))

## Get my data
clusters = []
concentrations = []

with open('aas_hist.txt', 'r') as f:
    lines = f.readlines()

for line in lines:
    clusters.append(line.strip('\n'))

mod = 0
strong = 0

for cluster in clusters:
    with open('results/{}/results/concentration/conc.txt'.format(cluster), 'r') as f:
        conc = float(f.readlines()[0].strip('\n'))
        concentrations.append(conc)
        if conc >= 0.075 and conc < 0.155:
            mod += 1
        elif conc >= 0.155:
            strong += 1

norm = np.full_like(np.array(concentrations), 1.0/len(clusters))

print('Moderate cool cores: {}/{}'.format(mod,len(clusters)))
print('Strong cool cores: {}/{}'.format(strong,len(clusters)))

plt.rcParams.update({'font.size': 12})

plt.plot(x, np.array(y)/47, label='CCCP low-z', color='purple', linewidth=3)
plt.hist(concentrations, weights=norm, histtype='step', label='CEREAL', color='orange', linewidth=3)
plt.xlabel(r'c$_{\mathrm{SB}}$')
plt.ylabel('Fraction')
plt.xlim(0, 0.4)
plt.legend()
plt.savefig('concentrations.png')

