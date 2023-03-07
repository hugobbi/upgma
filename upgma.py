import numpy as np
from newick import loads

def upgma(d, ids):
    while d.shape[0] > 1:
        n = d.shape[0]
        i, j = find_lowest_value(d)

        dist_node_to_u = d[i, j]/2
        new_node = f"({ids[i]}: {round(dist_node_to_u, 4)}, {ids[j]}: {round(dist_node_to_u, 4)})"
        new_ids = [new_node]
        for id in ids:
            if id not in {ids[i], ids[j]}:
                new_ids.append(id)

        dprime = np.zeros((n-1, n-1))
        ij_indexes = [i, j]
        dprime[1:, 1:] = np.delete(np.delete(d, ij_indexes, axis=1), 
                            ij_indexes, axis=0)
        
        dist_k_to_u = (d[i] + d[j])/2
        dist_k_to_u = np.delete(dist_k_to_u, ij_indexes) 
        dist_k_to_u = np.concatenate([[0], dist_k_to_u])
        
        for k in range(n-1):
           dprime[0, k] = dprime[k, 0] = dist_k_to_u[k]

        d = np.copy(dprime)
        ids = new_ids.copy()

    tree = ids[0]
    return tree

def find_lowest_value(d):
    lowest_value = np.Inf
    l, k = -1, -1

    for i in range(len(d)):
        for j in range(len(d[i])):
            if d[i, j] < lowest_value and i != j:
                lowest_value = d[i, j]
                l, k = i, j

    return (l, k)

d = np.array([[0.,  49.,  30.,  40.,  33.,  51., 124.,  33.,    39.,  23.],
              [ 49.,   0.,  58.,  43.,  38.,  49., 124.,  29.,  39.,  45.],
              [ 30.,  58.,   0., 29.,  25.,   36., 119.,  37.,  31.,  44.],
              [ 40.,  43.,  29.,   0.,  21.,  21., 134.,  29.,  37.,  42.],
              [ 33.,  38.,  25.,  21.,   0.,  12., 128.,  29.,  45.,  28.],
              [ 51.,  49.,  36.,  21.,  12.,   0., 144.,  29.,  47.,  52.],
              [124., 124., 119., 134., 128., 144.,   0., 121., 124., 121.],
              [ 33.,  29.,  37.,  29.,  29.,  29., 121.,   0.,  39.,  38.],
              [ 39.,  39.,  31.,  37.,  45.,  47., 124.,  39.,   0.,  37.],
              [ 23.,  45.,  44.,  42.,  28.,  52., 121.,  38.,  37.,   0.]])

ids = ["seq_1", "seq_2", "seq_3", "seq_4", "seq_5", "seq_6", "seq_7", "seq_8", "seq_9", "seq_10"]

tree = upgma(d, ids)
print(tree)
print(loads(tree)[0].ascii_art())