import numpy as np
import sys
sys.path.append('../src')
sys.path.append('src')
import bd


def test_lj():
    N, dim = 512, 2
    density = 0.5
    nblock, block = 10, 10000
    box = (N / density) ** (1 / dim)

    system_lj = bd.BDLJPBC(
        N, dim, box=box, dt=0.005,
        gamma=1, kT=1, m=1
    )
    #system_lj.r = np.loadtxt('cnf.inp', skiprows=2, usecols=[0, 1, 2])
    #system_lj.v = np.loadtxt('cnf.inp', skiprows=2, usecols=[3, 4, 5])

    #obs = bd.Thermodynamic(block=block)
    #system_lj.attach(obs)

    bd.animate(system_lj, r=10, jump=10)


if __name__ == "__main__":
    test_lj()
