import os
from tqdm import tqdm
import pickle
import itertools
import numpy as np
from noise_3d import add_vicsek_noise_3d
import fish_track as ft
from scipy import spatial
import matplotlib.pyplot as plt
import csimulate


def vicsek_move(positions, velocities, radius, eta, box):
    """
    move particles according to Vecsek model in 3D, see Vicsek et al PRL 1995

    Args:
        positions (np.ndarray): shape, (3, N)
        velocities (np.ndarray): the velocity vector in the spherical coordinate system shape, (3, N)
        radius (float): agents align with neighbours within the interaction radius
        eta (float): the noise parameter in Vicsek model
        box(np.ndarray): simulation box, (Lx, Ly, Lz)

    Return:
        (np.ndarray, np.ndarray): updated positions and velocities, shape ((3, N), (3, N))
    """
    bx, by, bz = box.ravel()
    neighbours = list(itertools.product([0, bx, -bx], [0, by, -by], [0, bz, -bz]))
    distances_4d = [spatial.distance.cdist(positions.T, positions.T + n) for n in neighbours]
    distances = np.min(distances_4d, axis=0)
    is_interact = distances <= radius

    interact_num = np.sum(is_interact, axis=0)

    n = positions.shape[1]

    azimuth_table = velocities[1].reshape((n, 1)) * np.ones((1, n))
    elevate_table = velocities[2].reshape((n, 1)) * np.ones((1, n))
    r = np.cos(elevate_table)

    x_mean = np.sum(r * np.cos(azimuth_table) * is_interact, axis=0) / interact_num
    y_mean = np.sum(r * np.sin(azimuth_table) * is_interact, axis=0) / interact_num
    z_mean = np.sum(np.sin(elevate_table) * is_interact, axis=0) / interact_num

    azimuth = np.arctan2(y_mean, x_mean)  # -pi, pi
    elevation = np.arctan(z_mean / np.sqrt(y_mean ** 2 + x_mean ** 2))  # -pi/2, pi/2

    azimuth, elevation = add_vicsek_noise_3d(azimuth, elevation, eta=eta)

    azimuth[azimuth >= np.pi] -= 2 * np.pi
    azimuth[azimuth < -np.pi] += 2 * np.pi
    elevation[elevation >= np.pi / 2] -= np.pi
    elevation[elevation < -np.pi / 2] += np.pi

    velocities[1] = azimuth
    velocities[2] = elevation

    dr = np.cos(elevation)
    dz = velocities[0] * np.sin(elevation)
    dx = velocities[0] * dr * np.cos(azimuth)
    dy = velocities[0] * dr * np.sin(azimuth)

    positions += np.array((dx, dy, dz))
    positions -= (positions // box) * box


def simulate(number, noise, equilibrium=3000, sample=15000, box=10, radius=1.0, speed=0.03):
    """
    Run a 3D vicsek model simulation

    Args:
        number (int): Number of particles/agents
        noise (float): The noise parameter in vicsek model, range from 0 to 1
        equilibrium (int): Number of frames running before samping so that the system go to steady state
        sample (int): Number of frames being sampled
        box (float): the cubic box size of the simulation with periodic boundary condition
        radius (float): the interaction radius of vicsek model
        speed (float): the speed of vicsek model, dt is 1 in simulation

    Returns:
        None, write simulation result in a pickled file to hard drive
    """
    box = np.ones([3, 1]) * box

    positions = np.random.rand(3, number) * box  # 3, n
    azimuth   = np.random.uniform(-np.pi, np.pi, size=number)
    elevation = np.random.uniform(-np.pi/2, np.pi/2, size=number)
    velocities = np.array([speed * np.ones(number), azimuth, elevation])

    for _ in tqdm(range(equilibrium)):
        vicsek_move(positions, velocities, radius, noise, box)

    positions_movie = np.empty((sample, number, 3))
    velocities_movie = np.empty((sample, number, 3))

    for i in tqdm(range(sample)):
        vicsek_move(positions, velocities, radius, noise, box)
        positions_movie[i] = positions.T
        velocities_movie[i] = velocities.T

    vr = np.cos(velocities_movie[:, :, 2])
    vx = (1 * (np.cos(velocities_movie[:, :, 1]) * vr))
    vy = (1 * (np.sin(velocities_movie[:, :, 1]) * vr))
    vz = (1 * np.sin(velocities_movie[:, :, 2]))
    order = np.sqrt(vx.mean(1)**2 + vy.mean(1)**2 + vz.mean(1)**2)
    return order.mean()


if __name__ == "__main__":
    import time
    N = 50
    density = 1.0
    box = (N / density) ** (1/3)
    spd = 0.1
    noises = np.arange(0.0, 1.2, 0.2)
    orders = []
    t0 = time.time()
    for reduced_noise in noises:
        noise = reduced_noise * np.sqrt(density)
        order = simulate(number=N, noise=noise, equilibrium=1000, sample=1000, box=box, radius=1, speed=spd )
        orders.append(order)
    print(f"Python version: {time.time() - t0:.2f}s")

    plt.scatter(noises, orders, color='tomato', facecolor='none')

    orders = []
    t0 = time.time()
    for reduced_noise in noises:
        noise = reduced_noise * np.sqrt(density)
        result = csimulate.vicsek_3d_pbc(N, box, 1.0, noise, spd, 1000, 1000, 1)  # (T, n, 6)
        order = np.sqrt((result[:, :, 3:].mean(1) ** 2).sum(-1)).mean()
        orders.append(order / spd)
    print(f"C++ version: {time.time() - t0:.2f}s")

    plt.scatter(noises, orders, color='tomato', marker='+')
    plt.show()
