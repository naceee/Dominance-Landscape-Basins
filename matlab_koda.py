import math
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.measurements import label


def find_paths_to_basins(M, labeled_basins, labeled_slopes, x, y, resolution):
    """
    funcion that recursively finds all the paths to basins
    :param M: num_objectives * resolution * resolution matrix, 
        holding objective vector for each mesh location
    :param labeled_basins: num_objectives * resolution * resolution matrix,
        with TODO
    :param labeled_slopes: num_objectives * resolution * resolution matrix,
        with TODO
    :param x: resolution by 1 array of ordered x locations of samples
        (e.g. x = linspace(-1, 1, resolution) )
    :param y: resolution by 1 array of ordered y locations of samples
        (e.g. y = linspace(-1, 1, resolution) )
    :param resolution: resolution = mesh size (number of cells on each dimension)
        to use when plotting the objective response
    :return:
    """
    if labeled_basins[x, y] != 0:
        labeled_slopes[x, y] = labeled_basins[x, y]
        return labeled_basins[x, y]
    elif labeled_slopes[x, y] != -1:
        return labeled_slopes[x, y]

    all_possible_basins = set()
    for i in [-1, 0, 1]:
        for j in [-1, 0, 1]:
            if 0 <= x + i < resolution and 0 <= y + j < resolution:
                if np.all(M[:, x, y] > M[:, x + i, y + j]):
                    b = find_paths_to_basins(M, labeled_basins, labeled_slopes, x + i, y + j, resolution)
                    all_possible_basins.add(b)
    if len(all_possible_basins) == 1:
        (b, ) = all_possible_basins
        labeled_slopes[x, y] = b
        return b
    else:
        labeled_slopes[x, y] = 0
        return 0


def enumerate_basins_and_slopes(M, resolution, neighbourhood="Moore"):
    basins = np.zeros((resolution, resolution))

    for i in range(resolution):
        for j in range(resolution):
            if is_local_min(M, i, j, neighbourhood, resolution):
                basins[i, j] = 1

    if neighbourhood == "Moore":
        structure = np.ones((3, 3), dtype=int)
    else:
        structure = np.matrix('0 1 0; 1 1 1; 0 1 0')

    labeled_basins, n_components = label(basins, structure)

    return labeled_basins, basins


def is_local_min(M, x, y, neighbourhood, resolution):
    if neighbourhood == "Moore":
        directions = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]
    else:
        directions = [(-1, 0), (0, -1), (0, 1), (1, 0)]

    for d in directions:
        if 0 <= x + d[0] < resolution and 0 <= y + d[1] < resolution:
            if np.all(M[:, x, y] > M[:, x + d[0], y + d[1]]):
                return False
    return True


def getDominanceLandscapeBasinsFromMatrix(M, x, y, neighbourhood):
    numObjectives, resolution, r = np.shape(M)

    if resolution != r:
        raise Exception('Second and third dimension of Y must be the same size')
    if len(x) != resolution:
        raise Exception('must be as many x grid labels as elements')
    if len(y) != resolution:
        raise Exception('must be as many y grid labels as elements')

    labeled_basins, basins = enumerate_basins_and_slopes(M, resolution, neighbourhood)
    labeled_slopes = np.zeros((resolution, resolution)) - 1

    for i in range(resolution):
        for j in range(resolution):
            labeled_slopes[i, j] = find_paths_to_basins(M, labeled_basins, labeled_slopes, i, j, resolution)

    basins = (labeled_slopes > 0) + basins

    return labeled_basins, basins, labeled_slopes


def create_matrix(f1, f2, lim_x, lim_y, resolution):

    M = np.zeros((2, resolution, resolution))
    x = np.linspace(lim_x[0], lim_x[1], resolution)
    y = np.linspace(lim_y[0], lim_y[1], resolution)
    neighbourhood = "Moore"

    for i, xx in enumerate(x):
        for j, yy in enumerate(y):
            M[0, j, i] = f1(xx, yy)
            M[1, j, i] = f2(xx, yy)

    labeled_basins, basins, labeled_slopes = getDominanceLandscapeBasinsFromMatrix(M, x, y, neighbourhood)

    slopes_only = np.ma.masked_where((0 >= labeled_slopes), labeled_slopes)
    basins_only = np.ma.masked_where((0 >= labeled_basins), labeled_basins)

    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(14, 6))
    axs[0].imshow(2 - basins, interpolation='none', cmap="gray")

    axs[1].imshow(slopes_only, alpha=0.3, interpolation='none')
    axs[1].imshow(basins_only, interpolation='none')

    plt.show()


def main():
    f0 = lambda x, y: (x - 0) ** 2 + (y - 1) ** 2
    f1 = lambda x, y: (x - 3) ** 2 + (y + 4) ** 2
    f2 = lambda x, y: ((x ** 2 - 5 * math.cos(2 * math.pi * x)) +
                       (y ** 2 - 5 * math.cos(2 * math.pi * y)))

    limits_x = (-5, 5)
    limits_y = (-5, 5)

    start = time.time()
    create_matrix(f1, f2, limits_x, limits_y, 501)
    end = time.time()
    print("elapsed time:", int(end - start), "s")


if __name__ == '__main__':
    main()


