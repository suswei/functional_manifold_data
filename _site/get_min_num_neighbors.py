from sklearn import manifold
from scipy.sparse.csgraph import connected_components


def get_min_num_neighbors(data):

    number_connected_components = 100
    num_neigh = 1
    while number_connected_components > 1:
        embedding = manifold.Isomap(n_neighbors=num_neigh)
        embedding.fit_transform(data)
        nbrs = embedding.nbrs_
        sparse_graph = nbrs.kneighbors_graph(data).toarray()
        (number_connected_components, labels) = connected_components(sparse_graph)
        num_neigh += 1

    best_num_neigh = num_neigh - 1
    

    return best_num_neigh



