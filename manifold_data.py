# TODO: find out where this code came from, some scms python package?

import numpy as np
import pandas as pd

def manifold_data(samplesize):

  # data sampling
  n_manifold = samplesize
  n_bg = samplesize
  sampling_sigma = .2

  # sample the manifold
  manifold_x = np.linspace(0, 4*np.pi, n_manifold)
  manifold_y = np.sin(manifold_x)

  manifold_x += np.random.normal(0, sampling_sigma, n_manifold)
  manifold_y += np.random.normal(0, sampling_sigma, n_manifold)

  manifold = np.column_stack((manifold_x, manifold_y))

  # add some background noise
  bg_x = np.random.uniform(-1, 14, n_bg)
  bg_y = np.random.uniform(-2, 2, n_bg)
  bg = np.column_stack((bg_x, bg_y))

  # combine the noise and the signal
  # data = np.vstack((manifold, bg))
  data = pd.DataFrame(manifold = manifold, noise = bg)

  return data
