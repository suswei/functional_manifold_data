import numpy as np


def make_anisotropic_gaussian_kernel(Hdiag):
    """
    Makes an isotropic gaussian kernel function with bandwidth sigma.
    """

    # the kernel function is a function of ||x-y||^2, the squared Euclidean
    # distance between x and y

    # hardcode bandwidth matrix H according to Silverman's
    H = np.diag(Hdiag)

    def kernel(diff):
        # diff should be dim by sample size
        diff = np.transpose(diff)
        d = diff.shape[0]
        prod = np.matmul(np.transpose(diff), np.linalg.inv(H))
        prod = np.matmul(prod, diff)
        weight_matrix = (2 * np.pi) ** (-d / 2)*np.linalg.det(H)**(-1/2) * np.exp( (-1/2)* prod)
        weights = np.diag(weight_matrix)
        return weights

    return kernel


def kernel_density_estimate(x, data, kernel):
    """
    Estimates the probability density at x, given the data and the kernel.
    This obviously does not use data-dependent kernels.
    """
    # first, we calculate the interpoint distances
    d_squared = np.sum((data - x) ** 2, axis=1)

    # and evaluate the kernel at each distance
    weights = kernel(data-x)

    # return the kernel density estimate
    return np.mean(weights)


def gradient(x, data, kernel, Hdiag):
    """
    Computes the gradient at x of the kde of the data using isotropic Guassian
    kernels with bandwidth sigma.
    """
    # the number of data points and the dimensionality
    n, d = data.shape

    # first, we form the u_i
    u = np.matmul((data - x),np.diag(1./Hdiag))

    # now we form the c_i
    d_squared = np.sum((data - x) ** 2, axis=1)
    c = kernel(data-x)

    return -1. * np.mean(c * u.T, axis=1)


def hessian(x, data, kernel, Hdiag):
    """
    Computes the Hessian at x of the kde of the data.
    """
    # the number of data points and the dimensionality
    n, d = data.shape

    # first, we form the u_i
    u = np.matmul((data - x),np.diag(1./Hdiag))

    # now we form the c_i
    # d_squared = np.sum((data - x) ** 2, axis=1)
    c = kernel(data-x)

    Z = np.diag(1./Hdiag)

    return 1. / n * (c * u.T).dot(u) - 1. / n * Z * np.sum(c)


def local_inv_cov(x, data, kernel, Hdiag):
    """
    Computes the local inverse covariance from the gradient and hessian.
    """

    p = kernel_density_estimate(x, data, kernel)
    g = gradient(x, data, kernel, Hdiag)[:, np.newaxis]
    h = hessian(x, data, kernel, Hdiag)
    return -1. / p * h + 1. / p ** 2 * g.dot(g.T)


def mean_shift_update(x, data, kernel):
    """
    Applies a kernel mean shift update to the point x, given the data and a
    1-d kernel function of the squared distance.
    """
    # first, calculate the interpoint distance
    d_squared = np.sum((data - x) ** 2, axis=1)

    # eqvaluate the kernel at each distance
    weights = kernel(data-x)

    # now reweight each point
    shift = (data.T.dot(weights)) / np.sum(weights)

    # return the shift update
    return shift - x


def subspace_constrained_mean_shift_update(x, data, kernel, Hdiag):
    """
    Compute the constrained mean shift update of the point at x.
    """
    # first, we evaluate the mean shift update
    msu = mean_shift_update(x, data, kernel)

    # next we compute the local inverse covariance
    lic = local_inv_cov(x, data, kernel, Hdiag)

    # now we perform the eigendecomposition of lic
    vals, vects = np.linalg.eig(lic)

    # sort the eigenvectors by increasing eigenvalue
    idx = np.argsort(vals)
    vects = vects[:, idx]

    # project the mean shift update onto the orthogonal subspace spanned
    # by the largest eigenvectors of the inverse covariance
    V = vects[:, 1:]
    return V.dot(V.T).dot(msu)


def scms(data, sigma=None, n_iterations=5):
    """
    Performs subspace constrained mean shift on the entire data set using
    a Gaussian kernel with bandwidth of sigma.
    """
    denoised = data.copy()
    n, d = data.shape
    if sigma is None:
        variance_min = min(np.var(denoised,axis=0))
        Hdiag = (4/(d+2))**(1/(d+4)) * n **(-1/(d+4)) * np.sqrt(variance_min *np.ones(d))
    else:
        Hdiag = sigma * np.ones(d)

    Hdiag = Hdiag ** 2
    kernel = make_anisotropic_gaussian_kernel(Hdiag)

    for j in range(n_iterations):
        for i in range(data.shape[0]):
            # denoised[i] += subspace_constrained_mean_shift_update(
            #        denoised[i], data, kernel, sigma)
            np.add(denoised[i], subspace_constrained_mean_shift_update(denoised[i], data, kernel, Hdiag),
                   out=denoised[i], casting="unsafe")
    return denoised
