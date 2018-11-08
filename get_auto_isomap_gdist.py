from sklearn import manifold
from scipy.sparse.csgraph import connected_components


def get_auto_isomap_gdist(data):

    number_connected_components = 100
    num_neigh = 1
    while number_connected_components > 1:
        embedding = manifold.Isomap(n_neighbors=num_neigh)
        embedding.fit_transform(data)
        nbrs = embedding.nbrs_
        sparse_graph = nbrs.kneighbors_graph(data).toarray()
        (number_connected_components, labels) = connected_components(sparse_graph)
        print(number_connected_components)
        num_neigh += 1

    best_num_neigh = num_neigh - 1
    print(best_num_neigh)
    embedding = manifold.Isomap(n_neighbors=best_num_neigh)
    embedding.fit_transform(data)
    gdist = embedding.dist_matrix_

    return gdist



