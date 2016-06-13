import numpy as np
from numpy import linalg


def freq_covariance(map1, map2, weight1, weight2, freq1, freq2, no_weight=False, 
        normalize=False):
    r"""Calculate the weighted nu nu' covariance"""
    input_map1 = map1[freq1, :]
    input_map2 = map2[freq2, :]
    if no_weight:
        input_weight1 = np.ones_like(input_map1)
        input_weight2 = np.ones_like(input_map2)
    else:
        input_weight1 = weight1[freq1, :]
        input_weight2 = weight2[freq2, :]

    map1shp = input_map1.shape
    map2shp = input_map2.shape

    map1_flat = input_map1
    map2_flat = input_map2
    weight1_flat = input_weight1
    weight2_flat = input_weight2

    wprod1 = map1_flat * weight1_flat
    wprod2 = map2_flat * weight2_flat

    # TODO: or should this be wprod2, wprod1^T?
    quad_wprod = np.ma.dot(wprod1, np.transpose(wprod2))
    quad_weight = np.ma.dot(weight1_flat, np.transpose(weight2_flat))

    if normalize:
        # For check C_{AB}/(C_A*C_B)**0.5
        quad_wprod1 = np.ma.dot(wprod1, np.transpose(wprod1))
        quad_wprod2 = np.ma.dot(wprod2, np.transpose(wprod2))
        quad_wprod /= (quad_wprod1*quad_wprod2)**0.5

    if not np.any(quad_wprod.mask):
        print "No mask in covariance matrix"
    else:
        print "%f%% covariance matrix are masked"%(
                float(np.ma.count_masked(quad_wprod))/quad_wprod.size * 100)

    mask = (quad_weight < 1e-20)

    quad_weight[mask] = 1.
    quad_wprod /= quad_weight
    quad_wprod[mask] = 0
    quad_weight[mask] = 0

    #return quad_wprod[..., np.newaxis], quad_weight[..., np.newaxis]
    return quad_wprod, quad_weight

def get_freq_svd_modes(corr, n_modes):
    r"""Same as get freq eigenmodes, but treats left and right maps
    separatly with an SVD.
    """
    u_matrix, singular_values, v_matrix = linalg.svd(corr)
    v_matrix = v_matrix.T
    sorted_singular_values = list(singular_values)
    sorted_singular_values.sort()
    left_vectors = []
    right_vectors = []
    for mode_index in range(n_modes):
        ind, = np.where(abs(singular_values) ==
                        sorted_singular_values[-mode_index - 1])
        #if len(ind) > 1:
        #    raise NotImplementedError('2 eigenvalues bitwise equal.')

        left_vectors.append(u_matrix[:, ind[0]])
        right_vectors.append(v_matrix[:, ind[0]])

    return singular_values, left_vectors, right_vectors
