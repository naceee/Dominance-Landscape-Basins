import math
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import label


def find_paths_to_basins(M, labeled_basins, labeled_slopes, x, y, resolution):
    """ function that recursively finds all the paths to basins
    """

    # if location (x, y) is part of a basin, return this basin index
    if labeled_basins[x, y] != 0:
        labeled_slopes[x, y] = labeled_basins[x, y]
        return labeled_basins[x, y]

    # if location (x, y) is part of a slope that leads to some basin, return that basin index
    elif labeled_slopes[x, y] != -1:
        return labeled_slopes[x, y]

    # check all neighbours that dominate this point (x, y)
    # and recursively find which basin they lead to
    leading_to_basins = set()
    for i in [-1, 0, 1]:
        for j in [-1, 0, 1]:
            if 0 <= x + i < resolution and 0 <= y + j < resolution:
                if np.all(M[:, x, y] > M[:, x + i, y + j]):
                    b = find_paths_to_basins(M, labeled_basins, labeled_slopes, x + i, y + j,
                                             resolution)
                    leading_to_basins.add(b)
    # if all dominating paths lead to the same basin b, then return b
    if len(leading_to_basins) == 1:
        (b, ) = leading_to_basins
        labeled_slopes[x, y] = b
        return b

    # if dominating paths lead to the more than one basin, then return 0
    else:
        labeled_slopes[x, y] = 0
        return 0


def enumerate_basins(M, resolution, neighbourhood="Moore"):
    """ find all the locally dominance-neutral regions, where all local moves estimated at the
    discretization used are mutually non-dominating and enumerate each region with numbers,
    starting from 1, ...
    """
    basins = np.zeros((resolution, resolution))

    # For each point, check the neighbourhood if there is any dominating neighbour.
    for i in range(resolution):
        for j in range(resolution):
            if is_local_min(M, i, j, neighbourhood, resolution):
                basins[i, j] = 1

    # for a given neighbourhood type, label all the connected components in the basins matrix
    if neighbourhood == "Moore":
        structure = np.ones((3, 3), dtype=int)
    else:
        structure = np.matrix('0 1 0; 1 1 1; 0 1 0')
    labeled_basins, n_components = label(basins, structure)

    return labeled_basins, basins


def enumerate_slopes(M, labeled_basins, resolution):
    """ find all the points in the discrete space that are locally dominated.
    With labeled_slopes[x, y] = r we denote points (x, y), for which all the dominating
    movement paths from neighbors lead to the same local optima region r.
    With labeled_slopes[x, y] = 0 we denote points (x, y), for which some dominating
    movement paths from neighbors lead to one local optima region and others to the other.
    """
    labeled_slopes = np.zeros((resolution, resolution)) - 1

    for i in range(resolution):
        for j in range(resolution):
            find_paths_to_basins(M, labeled_basins, labeled_slopes, i, j, resolution)

    return labeled_slopes


def is_local_min(M, x, y, neighbourhood, resolution):
    """ for a given neighbourhood type, check if point (x, y) is locally non-dominated """
    if neighbourhood == "Moore":
        directions = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]
    else:
        directions = [(-1, 0), (0, -1), (0, 1), (1, 0)]

    for d in directions:
        if 0 <= x + d[0] < resolution and 0 <= y + d[1] < resolution:
            if np.all(M[:, x, y] > M[:, x + d[0], y + d[1]]):
                return False
    return True


def get_dominance_landscape_basins_from_matrix(M, x, y, neighbourhood):
    num_objectives, resolution, r = np.shape(M)

    if resolution != r:
        raise Exception('Second and third dimension of Y must be the same size')
    if num_objectives < 2:
        raise Exception('Must have at least two objectives')
    if len(x) != resolution:
        raise Exception('must be as many x grid labels as elements')
    if len(y) != resolution:
        raise Exception('must be as many y grid labels as elements')

    #
    labeled_basins, basins = enumerate_basins(M, resolution, neighbourhood)
    labeled_slopes = enumerate_slopes(M, labeled_basins, resolution)
    return labeled_basins, basins, labeled_slopes


def create_plot(M, labeled_basins, basins, labeled_slopes, x, y):
    """ Create 2 plots:
    LEFT:
     - black: locally non-dominating regions
     - gray: regions that lead to the same local optima region
     - white: regions that lead to different optima regions

    RIGHT: ...
    """
    basins = (labeled_slopes > 0) + basins

    slopes_only = np.ma.masked_where((0 >= labeled_slopes), labeled_slopes)
    basins_only = np.ma.masked_where((0 >= labeled_basins), labeled_basins)

    # Set up plot
    fig1, ax = plt.subplots(nrows=1, ncols=2, subplot_kw=dict(projection='3d'))
    X, Y = np.meshgrid(x, y)
    ax[0].plot_surface(X, Y, M[0, :, :], rstride=1, cstride=1, linewidth=1, antialiased=False, shade=True)
    ax[1].plot_surface(X, Y, M[1, :, :], rstride=1, cstride=1, linewidth=1, antialiased=False, shade=True)
    plt.show()


    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(14, 14))
    axs[0].imshow(2 - basins, interpolation='none', cmap="gray", vmin=0, vmax=2,
                  extent=(x[0], x[-1], y[0], y[-1]))
    axs[1].imshow(slopes_only, alpha=0.3, interpolation='none', extent=(x[0], x[-1], y[0], y[-1]))
    axs[1].imshow(basins_only, interpolation='none', extent=(x[0], x[-1], y[0], y[-1]))

    plt.show()


def create_matrix(f1, f2, lim_x, lim_y, resolution):

    M = np.zeros((2, resolution, resolution))
    x = np.linspace(lim_x[0], lim_x[1], resolution)
    y = np.linspace(lim_y[0], lim_y[1], resolution)
    neighbourhood = "Moore"

    for i, xx in enumerate(x):
        for j, yy in enumerate(y):
            M[0, j, i] = f1(xx, yy)
            M[1, j, i] = f2(xx, yy)

    labeled_basins, basins, labeled_slopes = \
        get_dominance_landscape_basins_from_matrix(M, x, y, neighbourhood)

    create_plot(M, labeled_basins, basins, labeled_slopes, x, y)


def main():
    f0 = lambda x, y: (x + 2) ** 2 + (y + 3) ** 2
    f1 = lambda x, y: (x - 2) ** 2 + (y + 1) ** 2
    f2 = lambda x, y: ((x ** 2 - 5 * math.cos(2 * math.pi * x)) +
                       (y ** 2 - 5 * math.cos(2 * math.pi * y)))
    f3 = lambda x, y: math.sin(x) + math.cos(y)

    limits_x = (-4, 4)
    limits_y = (-4, 4)

    start = time.time()
    create_matrix(f1, f0, limits_x, limits_y, 401)
    end = time.time()
    print("elapsed time:", int(end - start), "s")


if __name__ == '__main__':
    main()


