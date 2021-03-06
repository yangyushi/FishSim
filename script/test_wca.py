import numpy as np
import sys
sys.path.append('../src')
sys.path.append('src')
import fish_sim as sim
import force


def test_lj():
    N, dim = 512, 2
    density = 0.5
    nblock, block = 10, 1000
    box = (N / density) ** (1 / dim)

    system_wca = sim.BDLJPBC(
        N, dim, box=box, dt=0.005,
        gamma=1, kT=1, m=1
    )
    system_wca.force_func = force.force_wca

    sim.animate(system_wca, r=10, jump=10, show=True)


if __name__ == "__main__":
    test_lj()
