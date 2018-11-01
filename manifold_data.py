# TODO: find out where this code came from, some scms python package?

import numpy as np
import pandas as pd

def manifold_data(samplesize):

  # data sampling
  n_manifold = int(samplesize)
  n_bg = int(samplesize)
  sampling_sigma = .2

  # sample the manifold
  manifold_x = np.linspace(0, 4*np.pi, n_manifold)
  manifold_y = np.sin(manifold_x)
  true_manix = manifold_x
  true_maniy = manifold_y
  true_mani = np.column_stack((true_manix, true_maniy))

  # this does not appear to be noise normal to the manifold
  manifold_x += np.random.normal(0, sampling_sigma, n_manifold)
  manifold_y += np.random.normal(0, sampling_sigma, n_manifold)
  manifold = np.column_stack((manifold_x, manifold_y))

  # add some background noise
  bg_x = np.random.uniform(-1, 14, n_bg)
  bg_y = np.random.uniform(-2, 2, n_bg)
  bg = np.column_stack((bg_x, bg_y))

  # combine the noisy manifold and true manifold
  trueplusnoise_stacked = np.vstack((manifold, true_mani))

  # TODO: can't pick up the pandas dataframe structure!
  # data = pd.DataFrame(data = manifold, temp = true_mani)

  return trueplusnoise_stacked
