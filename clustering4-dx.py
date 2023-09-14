import math, random
import numpy as np, pandas as pd, matplotlib.pyplot as plt
import scipy.spatial.distance as ssd
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import pdist, squareform

matfile = "pairwisedistancematrix_heuristicDX52_66.csv"
outfile = "clustersDX52_66.txt"
delta_level = 0.05
delta_threshold = 0.0 # will be set to average lowest 5% distance
cutoff = 0.995 #for proportion of intra-cluster distances less than delta_threshold
num_reps = 1000
min_clusters = 2
max_clusters = 64
random.seed(42)

#read distance matrix from .csv
#top left cell should be Seq_ID
#and first column should be sample names
df = pd.read_csv(matfile)
labels = df['Seq_ID'].tolist()
df = df.drop('Seq_ID', 1)

# convert square matrix form into a condensed nC2 array
distance_matrix = df.to_numpy()
num_strains = len(distance_matrix)
# Set the diagonals to zero
np.fill_diagonal(distance_matrix, 0)
distance_array = ssd.squareform(distance_matrix)

#heirarchical clustering
linkage_matrix = sch.linkage(distance_array, method='ward')

#plot, save as .png
plt.figure(figsize=(32, 16))
plt.title("test6")
dend = sch.dendrogram(linkage_matrix, labels=labels, show_leaf_counts=False)
ax = plt.gca()
ax.tick_params(axis='x', which='major', labelsize=24)
ax.tick_params(axis='y', which='major', labelsize=24)
plt.savefig('t.png')

def pool(clusters, nc):
    """
    input list of cluster assignments for the strains
    and number of clusters
    output a list of lists of members of each cluster
    """
    clust_pool = [ [] for _ in range(nc)]
    for i in range(len(clusters)):
        cnum = clusters[i][0]
        clust_pool[cnum].append(i)
    return clust_pool

#Determine max number of clusters that contain more than one strain each
for max_num_clust in range(min_clusters, max_clusters):
    clusters = sch.cut_tree(linkage_matrix, max_num_clust)
    counts = [0 for _ in range(max_num_clust)]
    for c in clusters:
        counts[c[0]] += 1
    minct = min(counts)
    if minct < 2:
        max_num_clust -= 1
        break
print("max number of clusters is", max_num_clust)
clusters = sch.cut_tree(linkage_matrix, max_num_clust)
clust_pool = pool(clusters, max_num_clust)

#assign species by first split 1= cayetanensis 2= ashfordi
#assuming cayetanensis is the larger population
species_arr = sch.cut_tree(linkage_matrix, 2)
species_arr += 1 #1 or 2
species = [species_arr[i][0] for i in range(len(species_arr))]

#Determine delta threshold
target = 3 #species 1, 2 or 3 for both together
n_choose = 1000 #find smallest cluster (probably 2)
for c in clust_pool:
    n_choose = min(n_choose, len(c))
num_sampled = 0
for cl in clust_pool:
    if (species[cl[0]] & target) > 0:
        num_sampled += n_choose
compares = num_sampled * (num_sampled - 1) // 2
rank = math.floor(delta_level * float(compares))
print("at", delta_level, ", rank is", rank, "of", compares, "distances.")

delta_threshold = np.zeros(compares)
for _ in range(num_reps):
    dlist = np.zeros(compares)
    samples = []
    n = 0
    for cl in clust_pool:
        if (species[cl[0]] & target) > 0:
            items = random.sample(cl, n_choose)
            samples.extend(items)
    k = 0
    for i in range(len(samples)):
        a = samples[i]
        for j in range(i+1, len(samples)):
            b = samples[j]
            dlist[k] = distance_matrix[a][b]
            k += 1
    sorted_array = np.sort(dlist)
    delta_threshold += sorted_array
delta_threshold /= float(num_reps)
print("delta threshold = ", delta_threshold[rank])

#Find optimal number of clusters
def check_cluster(delta_threshold):
    """
    input threshold distance
    output the minimal number of clusters
    where the number of intracluster distances under the delta_threshold
    is greater than "cutoff" proportion.
    Also returns cluster assigment list.
    """
    for num_clust in range(min_clusters, max_clusters):
        clusters2 = sch.cut_tree(linkage_matrix, num_clust)
        clust_pool2 = pool(clusters2, num_clust)
        count_all = 0
        count_high = 0
        for cl in clust_pool2:
            for i in range(len(cl)):
                a = cl[i]
                for j in range(i+1, len(cl)):
                    b = cl[j]
                    dist = distance_matrix[a][b]
                    count_all += 1
                    if dist > delta_threshold:
                        count_high += 1
        #print(num_clust, count_high, count_all)
        if 1.0 - float(count_high) / float(count_all) > cutoff:
            break
    return num_clust, clusters2

num_clust, clusters2 = check_cluster(delta_threshold[rank])
print("Optimal number of clusters at delta threshold is:",num_clust)

#output cluster assignment for each strain
out_file = open(outfile, 'w')
clust = []
for i in range(num_strains):
    clust.append((clusters2[i][0] + 1, labels[i]))
clust.sort()
for c in clust:
    out = c[1] + "\t" + str(c[0]) + "\n"
    out_file.write(out)
    #print(c[1], c[0])
out_file.close()

exit(0)
### extra plot stuff ###########################
#if you want to try a range of delta thresholds:
x = []
y = []
for i in range(math.ceil(0.4 * compares)):
    num_clust, _ = check_cluster(delta_threshold[i])
    delta_level1 = (i+1.0) / float(compares)
    print(i+1,  delta_level1, delta_threshold[i], num_clust)
    x.append(num_clust)
    y.append(delta_threshold[i])

#if you want to plot the distribution of distances
#in the original matrix to check normality
dlist = []
for i in range(num_strains):
    for j in range(i + 1, num_strains):
        d = round(distance_matrix[i][j], 4)
        dlist.append(d)

#or for distribution of the sampled distances
target = 3
dlist = np.zeros(compares)
samples = []
n = 0
for cl in clust_pool:
    if (species[cl[0]] & target) > 0:
        items = random.sample(cl, n_choose)
        samples.extend(items)
k = 0
for i in range(len(samples)):
    a = samples[i]
    for j in range(i+1, len(samples)):
        b = samples[j]
        dlist[k] = distance_matrix[a][b]
        k += 1
