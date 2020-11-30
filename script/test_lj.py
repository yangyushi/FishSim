import numpy as np
import sys
sys.path.append('../src')
sys.path.append('src')
import simulate as sim


def test_lj():
    N, dim = 256, 3
    density = 0.75
    nblock, block = 10, 1000
    box = (N / density) ** (1 / dim)

    system_lj = sim.BDLJPBC(
        N, dim, box=box, dt=0.005,
        gamma=1, kT=1, m=1
    )

    system_lj.r = np.loadtxt('cnf.inp', skiprows=2, usecols=[0, 1, 2])
    system_lj.v = np.loadtxt('cnf.inp', skiprows=2, usecols=[3, 4, 5])

    obs = sim.Thermodynamic(block=block)
    system_lj.attach(obs)

    sim.animate(system_lj, r=10, jump=10, show=True)


if __name__ == "__main__":
    test_lj()
