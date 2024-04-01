import numpy as np
from numba import njit, prange


@njit(parallel=True)
def ndot1(A, B):

    """
    Calculate C with
        C[i, j, k] = sum_q( A[i, j, q] · B[i, q, k] )

    Args:
        A (np.ndarray): shape is (N, 3, 3)
        B (np.ndarray): shape is (N, 3, 3)

    Return:
        np.ndarray: shape is (N, 3, 3)
    """
    C = np.zeros(A.shape)
    for i in prange(A.shape[0]):
        for j in range(3):
            for k in range(3):
                for q in range(3):
                    C[i, j, k] += A[i, j, q] * B[i, q, k]
    return C


@njit(parallel=True)
def ndot2(A, B):
    """
    Calculate C with
        C[i, j] = sum_q( A[i, j, q] · B[i, q] )

    Args:
        A (np.ndarray): shape is (N, 3, 3)
        B (np.ndarray): shape is (N, 3)

    Return:
        np.ndarray: shape is (N, 3)
    """
    C = np.zeros(B.shape)
    for i in prange(A.shape[0]):
        for j in range(3):
            for q in range(3):
                C[i, j] += A[i, j, q] * B[i, q]
    return C


def add_vicsek_noise_3d(phi, theta, eta):
    """
    Generate 3D noise for vicsek model, the procedure is
        1. Generate N random vectores inside a cap on the north pole.
           This corresponds to noise for vector (0, 0, 1)
           The area of the cap is a fraction of the surface area
               of the unit sphere
           A_cap / A_sphere = eta
        2. Rotate the generated noise with matrix R = F^-1 @ G @ F
           R rotates (0, 0, 1) to direction of different velocity vectors

    Args:
        phi (np.ndarray): the azimuth of velocity vectors, shape (N, )
        theta (np.ndarray): the elevation of velocity vectors, shape (N, )
        eta (float): the noise for the simulation, ranges from 0 to 1

    Returns:
        (np.ndarray, np.ndarray): ( noisy_azimuth(phi), noisy_elevation(theta))

    Example:
        >>> phi = np.array((0, 0))
        >>> theta = np.array((np.pi/4, -np.pi/4))
        >>> phi_noisy, theta_noisy = add_vicsek_noise_3d(phi, theta, 0.001)
        >>> np.allclose(phi_noisy, phi, atol=1e-1)
        True
        >>> np.allclose(theta_noisy, theta, atol=1e-1)
        True
    """
    N = phi.shape[0]
    r, z = np.cos(theta), np.sin(theta)
    x, y = r * np.cos(phi), r * np.sin(phi)  # shape N

    B = np.array((x, y, z)).T  # shape (N, 3)
    A = np.array([(0, 0, 1)] * N)  # shape (N, 3)

    noise_phi = np.random.uniform(-np.pi, np.pi, N)
    noise_z = np.random.uniform(1 - 2 * eta, 1, N)
    noise_r = np.sqrt(1 - noise_z**2)
    noise_x = noise_r * np.cos(noise_phi)
    noise_y = noise_r * np.sin(noise_phi)
    noise = np.array((noise_x, noise_y, noise_z)).T  # shape (N, 3)
    # A x B -> (-y, x, 0), |A x B| -> r
    G = np.array((  # rotate (0, 0, 1) to (x, y, z)
        (z, -r, [0] * N),
        (r, z, [0] * N),
        ([0] * N, [0] * N, [1] * N)
    ))  # (3, 3, N)
    G = np.moveaxis(G, -1, 0)  # (N, 3, 3)
    u = A  # (N, 3)
    v = B - (A * B).sum(axis=1)[:, np.newaxis] * A  # (N, 3)
    v /= np.linalg.norm(v, axis=1)[:, np.newaxis]

    w = np.cross(B, A, axis=1)
    F = np.stack((u, v, w), axis=1)  # shape (N, 3, 3)

    F_inv = np.linalg.inv(F)  # same as F_inv[i] = np.linalg.inv(F[i])
    rot = ndot1(F_inv, ndot1(G, F))
    noise = ndot2(rot, noise)  # apply rotation

    noise_phi = np.arctan2(noise[:, 1], noise[:, 0])
    noise_theta = np.arcsin(noise[:, 2])

    return noise_phi, noise_theta


def to_xyz(velocities):
    v0, phi, theta = velocities
    x = v0 * np.cos(theta) * np.cos(phi)
    y = v0 * np.cos(theta) * np.sin(phi)
    z = v0 * np.sin(theta)
    return np.vstack((x, y, z))
