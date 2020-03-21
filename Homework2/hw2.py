'''
@Description: 
@Author: Zhaoxi Chen
@Github: https://github.com/FrozenBurning
@Date: 2020-03-20 22:33:08
@LastEditors: Zhaoxi Chen
@LastEditTime: 2020-03-21 13:39:40
'''
import pandas as pd
from matplotlib import pyplot as plt
import scipy.cluster.hierarchy as hclst
from scipy.spatial.distance import pdist
import numpy as np
import seaborn as sns

visualized = True

# Read in Genematrix and preproc
with open("./GeneMatrix.txt", 'r') as f:
    mat = f.readlines()
raw_data = []
for data in mat:
    raw_data.append(data.split('\t')[1:])
genematrix_header = mat[0].split('\t')
genematrix_header[-1]=genematrix_header[-1].strip('\n')
genematrix = np.transpose(np.array(raw_data[1:])).astype('float')

# Read in clinical data, parsing ER_Status_nature2012 as criterion
df = pd.read_table('./clinical_data', sep='\t')
sampleid = df['sampleID'].tolist()
ER_measurement = df['ER_Status_nature2012']


def accuracy(gene2):
    true_pred = 0
    total = 0
    for i in range(0, len(gene2)):
        hot_id = genematrix_header[i]
        if hot_id in sampleid:
            idx = sampleid.index(hot_id)
            total = total+1
            if gene2[i] == 0 and ER_measurement[idx] == 'Positive':
                true_pred = true_pred + 1
            elif gene2[i] == 1 and ER_measurement[idx] == 'Negative':
                true_pred = true_pred + 1
    return true_pred/total

def myPCA(mat,pc_order = 1):
    cov = np.cov(mat,rowvar=False)
    eigen,eig_vec = np.linalg.eig(cov)
    print("top 5 PCs: ",eigen[:5])
    feature_t = np.transpose(eig_vec[:,:pc_order])
    dim_reduced = np.transpose(np.matmul(feature_t,np.transpose(genematrix)))
    return dim_reduced


# hierarchy cluster
cluster = hclst.linkage(genematrix, 'average')

# draw dendrogram and heatmap
if visualized:
    plt.figure(figsize=(50, 10))
    plt.title('Hierarchical Clustering')
    plt.xlabel('sample idx')
    plt.ylabel('distance')
    hclst.dendrogram(cluster, leaf_rotation=90., leaf_font_size=8.)
    # plt.show()
    plt.savefig('HClust.png')
    plt.figure()
    plt.title('Heatmap')
    plt.xlabel('gene')
    sns.clustermap(genematrix, figsize=(50, 50))
    plt.savefig('heatmap.png')

# classify into 2 classes
gene_in_2 = hclst.cut_tree(cluster, n_clusters=2)
print("Classification Accuracy: ",accuracy(gene_in_2))

# Using PCA and classify into 2 classes
reduced = myPCA(genematrix)
new_cluster = hclst.linkage(reduced, 'average')
if visualized:
    plt.figure(figsize=(50, 10))
    plt.title('Hierarchical Clustering Using top-1 PC')
    plt.xlabel('sample idx')
    plt.ylabel('distance')
    hclst.dendrogram(new_cluster, leaf_rotation=90., leaf_font_size=8.)
    # plt.show()
    plt.savefig('HClust-t1.png')
new_gene_in_2 = hclst.cut_tree(new_cluster,n_clusters=2)
print("Classification Accuracy with the first pinciple component:",accuracy(new_gene_in_2))