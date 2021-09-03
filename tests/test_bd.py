import sys
sys.path.append('../lib')
import numpy as np
from fish_sim.model import BD
from numba import njit
import matplotlib.pyplot as plt

@njit
def get_msd(trajectories, size, step=1):
    msd = np.empty((len(trajectories), size))
    for i, traj in enumerate(trajectories):
        length = len(traj)
        for tau in np.arange(0, size):
            if tau < length:
                msd[i, tau] = np.sum(
                    (traj[tau : length] - traj[: length - tau]) ** 2,
                    axis=1
                ).mean()
            else:
                msd[i, tau] = np.nan
    return np.arange(0, size), msd

def test_bd(plot=False):
    N, dim = 100, 4
    T, dt = 1000, 0.01
    D, kT, m = 1, 1, 1
    x = BD(N=N, dim=dim, dt=dt, D=D, kT=kT, m=m)
    trajs = np.empty((T, N, dim))
    for t in range(T):
        x.move_overdamp()
        trajs[t] = x.r
    trajs = np.moveaxis(trajs, 1, 0)
    tau, msd = get_msd(trajs, 50)
    msd = np.mean(msd, axis=0)
    theory = 2 * D * dim * tau * dt
    assert np.allclose(msd, theory, rtol=0.1), "Simulation and theory do not match"
    if plot:
        plt.scatter(
            tau, msd, color='tomato', fc='w',
            label=f'BD Simulation ({dim}D)'
        )
        plt.plot(
            tau, 2 * D * dim * tau * dt, color='teal',
            label=f'Theory ({dim}D)'
        )
        plt.title('Ideal Gas')
        plt.xlabel('$\\tau$')
        plt.ylabel('msd')
        plt.legend()
        plt.tight_layout()
        plt.show()
