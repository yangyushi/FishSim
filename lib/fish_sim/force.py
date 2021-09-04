import numpy as np


def force_lj(rij):
    """
    Calculate the force on particles with Lenneard-Jones potential

    Refernce:

        https://github.com/Allen-Tildesley/examples/blob/master/\
        python_examples/md_lj_module.py

    Args:
        rij (np.ndarray): the separation between particles, shape (dim, N, N)

    Return:
        tuple: the force (np.ndarray with shape (dim, N)) and
            other quantities of interest
    """
    r_cut = 2.5  # usual constant for LJ system
    pot_cut = r_cut ** -12 - r_cut ** -6
    sr2_overlap = 1.77  # overlap threshold, ~ r = 0.752 sigma

    rij_sq = np.sum(rij ** 2, axis=0)  # squared pairwise distance, shape (N, N)
    np.fill_diagonal(rij_sq, np.inf)
    in_range  = rij_sq < r_cut ** 2

    sr2 = np.where(in_range, 1.0 / rij_sq, 0.0)

    if (sr2 > sr2_overlap).any(): raise RuntimeError("Particle Overlap")

    sr6  = sr2 ** 3
    sr12 = sr6 ** 2
    cut  = sr12 - sr6
    vir  = cut + sr12
    pot  = np.where(in_range, cut - pot_cut, 0.0) # LJ pair potential (cut-and-shifted)
    f    = vir * sr2  # shape (N, N)
    f    = rij * f[np.newaxis, :, :]  # shape (dim, N, N)
    lap  = (22.0 * sr12 - 5.0 * sr6) * sr2  # shape (N, N)
    f    = np.sum(f, axis=2) * 24   # sum over columns of the displacement matrix

    interest = {
        "potential": np.sum(pot)  * 2.0,    # 4 / 2, because (ij - ji) double counted
        "virial":    np.sum(vir)  * 4.0,    # 24 / 3 / 2
        "force_sq":  np.sum(f**2),
        "laplacian": np.sum(lap)  * 24.0,
    }

    return f, interest  # sum over neighbours


def force_wca(rij):
    """
    Calculate the force on particles with Lenneard-Jones potential

    Refernce:

        https://github.com/Allen-Tildesley/examples/blob/master/\
        python_examples/md_lj_module.py

    Args:
        rij (np.ndarray): the separation between particles, shape (dim, N, N)

    Return:
        tuple: the force (np.ndarray with shape (N, dim)) and
            other quantities of interest
    """
    r_cut = 1.122462048  # 2^(1/6)
    pot_cut = r_cut ** -12 - r_cut ** -6
    sr2_overlap = 1.77  # overlap threshold, ~ r = 0.752 sigma

    rij_sq = np.sum(rij ** 2, axis=0)  # squared pairwise distance, shape (N, N)
    np.fill_diagonal(rij_sq, np.inf)
    in_range  = rij_sq < r_cut ** 2

    sr2 = np.where(in_range, 1.0 / rij_sq, 0.0)

    if (sr2 > sr2_overlap).any(): raise RuntimeError("Particle Overlap")

    sr6  = sr2 ** 3
    sr12 = sr6 ** 2
    cut  = sr12 - sr6
    vir  = cut + sr12
    pot  = np.where(in_range, cut - pot_cut - 1, 0.0) # WCA potential (cut-and-shifted)
    f    = vir * sr2  # shape (N, N)
    f    = rij * f[np.newaxis, :, :]  # shape (dim, N, N)
    lap  = (22.0 * sr12 - 5.0 * sr6) * sr2  # shape (N, N)

    interest = {
        "potential": np.sum(pot)  * 2.0,    # 4 / 2, because (ij - ji) double counted
        "virial":    np.sum(vir)  * 4.0,    # 24 / 3 / 2
        "force_sq":  np.sum(f**2) * 576.0, # 24^2
        "laplacian": np.sum(lap)  * 24.0,
    }

    return np.sum(f, axis=2) * 24, interest  # sum over neighbours
