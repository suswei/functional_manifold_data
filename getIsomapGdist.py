from sklearn import manifold

def getIsomapGdist(data, num_neigh):

  # need to do this because reticulate stupidly converts integers to floats, and then complains about it
  num_neigh = int(num_neigh)

  embedding = manifold.Isomap(n_neighbors = num_neigh, path_method = 'FW')
  embedding.fit_transform(data)
  gdist = embedding.dist_matrix_

  return gdist
