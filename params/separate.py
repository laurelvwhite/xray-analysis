archival = []

with open('archival_clusters.txt', 'r') as f:
    lines = f.readlines()

for line in lines:
    name = line.split(' ')[0]
    archival.append(name)

with open('all_clusters.txt', 'r') as f:
    lines = f.readlines()

with open('new_clusters.txt', 'w') as f:
    for line in lines:
        name = line.split(' ')[0]
        if name not in archival:
            f.write(line)
