import matplotlib.pyplot as plt

x = [0, 0.26, 0.52, 1.37, 2.9]
y = [1.612801337466251, 3.1042748222709338, 4.205390052667621, 19.32244413280197, 8.737896360077075]
yerr = [0.7529735956243392, 0.9853799743769246, 0.6714865934918186, 3.3952459561756636, 0.0017937569367731953]


## w/o bkg
new_y = [2.7897026971179097[
new_yerr = [0.24212140829400225]

fig, ax = plt.subplots()
ax.errorbar(x, y, yerr, fmt='bo')
ax.errorbar(x, new_y, new_yerr, fmt='co')
ax.plot(x, x, 'r')
ax.plot([x[0],x[-1]], [y[0],y[-1]], 'g')
ax.set_xlabel('Nurgaliev values')
ax.set_ylabel('My values')
plt.savefig('compare_aphot.png')


