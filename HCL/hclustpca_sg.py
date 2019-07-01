
# Importing the libraries
from Grover.GroverOptimize import *
from time import time, sleep
import feather
from copy import deepcopy
import matplotlib.pyplot as plt
import pandas as pd
import scipy.cluster.hierarchy as sch
import re
import quantumtools.Qconfig as Qconfig
import quantumtools.qtools as qt
import quantumtools.accessories as acs
from qiskit import QuantumProgram, register, available_backends
from ete3 import PhyloTree
import re

print("imported libraries\n")
#filepath = "/Users/sgujja/OneDrive - NextCODE Health/Omar_Code/quantumMachineLearning - Documents/Cuts/pc1_top_44_genes_train_matrix.feather"
filepath = "lumAB_pc1_44.feather"

def unravel_matrix_uppertri(matrix):
    """unravels upper triangular elements of a matrix to a vector"""
    unravelled = []
    n = len(matrix)
    # insufficient to check just first row, should check other rows too
    assert len(matrix[0]) == n, "must be a square matrix"

    for i in range(n):
        for j in range(i + 1, n):
            unravelled.append(matrix[i][j])
    return unravelled


def unravel_index_uppertri(n, ind):
    """takes single index of upper triangular elements (row by row left to right counting) and unravels it to (row,col)
    n is dimension of square matrix, ind is the index to unravel
    i.e. takes the index in the unraveled upper triangular vector and returns coordinate (row,col) it came from in the
    nxn matrix """

    # get index of the row containing element ind. Also equal to the number of full rows all
    # previous elements occupy
    i = np.floor(n - 1 / 2 - np.sqrt((n - 1 / 2) ** 2 - 2 * ind))
    # column index must be measured in steps to the right from the diagonal,
    # term in brackets
    j = i + 1 + (ind - (i * n - i * (i + 1) / 2))
    return int(i), int(j)


def ravel_index_uppertri(n, i, j):
    """takes (row,col) and ravels it to single index of upper triangular elements (row by row left to right counting)
    n is dimension of square matrix, ind is the index to unravel"""
    return int(i * n - i * (i + 1) / 2) + (j - (i + 1))

def ward_dist(d_vs,d_vt,d_st,n_v,n_s,n_t):
    """calculate d(u,v) ward distance between new cluster u and existing cluster v,
    where cluster u results from merger of s and t.
    Based on: https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html
    d_vs = d(v,s), distance between u and v, etc.
    n_v = number of points in cluster v, etc.
    """
    T = n_v + n_s + n_t
    d2 = ((n_v + n_s) * d_vs ** 2 + (n_v + n_t) * d_vt ** 2 - n_v * d_st ** 2) / T

    return np.round(np.sqrt(d2),3)

class Cluster:
    """Cluster class used in the algorithm
        # from: https: // github.com / rflynn / python - examples / blob / master / src / stats / cluster / agglomerative.py
    """
    def __init__(self):
        pass

    def __repr__(self):
        return '(%s,%s)' % (self.left, self.right)

    def add(self, clusters, grid, cluster_sizes, lefti, righti):
        self.left = clusters[lefti]
        self.right = clusters[righti]
        d_st = grid[lefti][righti]
        n_s = cluster_sizes[lefti]
        n_t = cluster_sizes[righti]

        # merge columns grid[row][righti] and row grid[righti] into corresponding lefti
        # update righti column while popping lefti column
        # update right i row (from right i column), and pop lefti row
        for k in range(len(grid)): # loop rows
            r = grid[k]
            # print("k ", k)
            # print("r ", r)
            # print("i,j ", lefti, " ", righti)
            if k != lefti and k != righti: # don't update entry for the two merged rows
                d_vs=r[lefti]
                d_vt=r[righti]
                n_v=cluster_sizes[k]
                r[lefti]=ward_dist(d_vs,d_vt,d_st,n_v,n_s,n_t)
                # print("r[lefti] ",r[lefti])
            r.pop(righti)
        # sleep(2)
        grid.pop(righti)

        for k in range(len(grid)): #shorter by 1 than last for loop
            grid[lefti][k] = grid[k][lefti]
        clusters.pop(righti)
        return (clusters, grid)


def agglomerate(labels, grid, mode=0,backendcode=0):
    """
    given a list of labels and a 2-D grid of distances, iteratively agglomerate
    hierarchical Cluster
    mode = 0 for classical minimize, 1 for quantum, 2 for quantum with recursive split
    """
    clusters = labels
    cluster_sizes = [1] * len(labels)

    # properties of hclust object in R, that are eventually exported .. Note that R indices start from 1
    merge = []
    height = []
    # similar to clusters object above, but keeps track in the hclust method, with leaves negative previous clusters positive
    merge_clusters = list(reversed(range(-len(labels),0)))  # start (-1, -2, ..., -len(labels) )


    while len(clusters) > 1:
        n = len(grid)
        # find two closest clusters
        print(clusters)
        distances = [(1, 0, grid[1][0])]
        # print([len(r) for r in grid])

        for i, row in enumerate(grid[2:]):
            distances += [(i + 2, j, c) for j, c in enumerate(row[:i + 2])]

        if mode == 0: # classical
            j, i, _ = min(distances, key=lambda x: x[2])
            # print("distances ", distances)
        else: # quantum
            values = unravel_matrix_uppertri(grid)
            if len(values)==1:
                (i,j)=(0,1)
            else:
                if mode == 1:
                    ind = groverminimize(opt_oracle_lessthan, backendcode=backendcode, withancilla=True, fvalues=values,silent=True)
                elif mode ==2:
                    ind = groverminimize_split(opt_oracle_lessthan, backendcode=backendcode, withancilla=True, fvalues=values,
                                               silent=True)
                (i, j) = unravel_index_uppertri(n, ind)

        # merge i<-j
        # print("i,j ",i," ",j)
        c = Cluster()
        h=grid[i][j]
        clusters, grid = c.add(clusters, grid, cluster_sizes,i, j)
        clusters[i] = c
        cluster_sizes[i] += cluster_sizes.pop(j)
        # update R hclust objects, careful with indices starting at 0 or 1. Python indices (i,j,nc) keep at 0+, but the merge object
        # should be transferred to R without modification, its indices (left, right, nc+1)
        nc = len(merge) # number of clusters so far
        left = merge_clusters[i]
        right = merge_clusters.pop(j)
        merge.append([left,right])
        merge_clusters[i]=nc+1 # index of new cluster
        height.append(h)

        # print("cluster_sizes ",cluster_sizes)
    return clusters.pop(), merge, height


def hier_cluster2(doQuantum=True,mode=2, backendcode=0):
    """hierarchichal clustering on rows (samples) instead of columns"""
    # to do: hier_cluster2 should be merged into hier_cluster, since the two are almost identical, but for the data
    # pre processing, and calculation of the distance matrix. Also loading of the data should be more customizable

    if "{}" in filepath: # which dataset to use ... this part as well as filepath above are customizable
        dataset_train = feather.read_dataframe(filepath.format('lumAB', 'train'))#.values
        X_train = dataset_train.iloc[:, 2:].values
        y_train = dataset_train.iloc[:, 1].values
    else:
        dataset_train = feather.read_dataframe(filepath)#.values
        X_train = dataset_train.iloc[:, 2:].values
        y_train = dataset_train.iloc[:, 1].values

    nindices = len(X_train)
    indices=list(range(nindices))
    print("Distance matrix of ", nindices, "samples")
    dist_matrix = [[np.round(np.linalg.norm(X_train[i] - X_train[j]), 1) for i in indices] for j in indices]

    if doQuantum:
        indices_q = indices.copy()
        dist_matrix_q=deepcopy(dist_matrix)

    # need to map ntrial*(ntrial-1)/2 distances to the numbers 0 to 2**nqubits-1
    # nqubits = int(np.ceil(num_qubits_from_samples(nindices)))
    # use the same number of qubits to represent each index. Note for some ntrial, you can save several qubit with
    # more complex mappings, e.g. ntrial = 11, each index 4 qubits (2**4>11) or you can use 7 qubits total (2**7 > 11^2)
    # then get rid of the duplicate half distances, since d(i,j)=d(j,i)
    start_time = time()
    clusters_c, merge_c, height_c=agglomerate(indices, dist_matrix,0)
    order_c=list(map(lambda x:int(x)+1, re.findall(r'\d+', clusters_c.__str__()))) # add 1 so order matches R indices
    #print(order_c)
    tree_c = PhyloTree(str(clusters_c)+";")

    print(clusters_c)
    print("merge_c ",merge_c)
    print("height_c ", height_c)
    print("Classical HCl takes ", round(time() - start_time, 4), "s\n")

    feather.write_dataframe(pd.DataFrame(merge_c),"temp/merge_c.feather")
    feather.write_dataframe(pd.DataFrame(height_c),"temp/height_c.feather")
    feather.write_dataframe(pd.DataFrame(order_c),"temp/order_c.feather")

    if doQuantum:
        start_time = time()

        clusters_q,merge_q,height_q=agglomerate(indices_q, dist_matrix_q,mode=mode,backendcode=backendcode)
        order_q = list( map(lambda x: int(x) + 1, re.findall(r'\d+', clusters_q.__str__())))  # add 1 so order matches R indices

        tree_q = PhyloTree(str(clusters_q) + ";")

        print(clusters_q)
        print("merge_q ",merge_q)
        print("height_q ", height_q)
        print("Quantum HCl takes ", acs.timestring(time() - start_time),"\n")

        qctree_sim = tree_q.compare(tree_c)['source_edges_in_ref']
        print("Similarity score, quantum to classic:", round(qctree_sim,4))

        feather.write_dataframe(pd.DataFrame(merge_q), "temp/merge_q.feather")
        feather.write_dataframe(pd.DataFrame(height_q), "temp/height_q.feather")
        feather.write_dataframe(pd.DataFrame(order_q), "temp/order_q.feather")

        return [clusters_c,clusters_q]

hier_cluster2()
