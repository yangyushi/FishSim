import sys
import numpy as np
sys.path.insert(0, '../lib')
import fish_sim as fs


def test_lj():
    N, dim = 256, 3
    density = 0.75
    block = 100
    box = (N / density) ** (1 / dim)

    system_lj = fs.model.BDLJPBC(
        N, dim, box=box, dt=0.005,
        gamma=1, kT=1, m=1
    )

    system_lj.r = np.loadtxt('cnf.inp', skiprows=2, usecols=[0, 1, 2])
    system_lj.v = np.loadtxt('cnf.inp', skiprows=2, usecols=[3, 4, 5])

    obs = fs.utility.Thermodynamic(block=block)
    system_lj.attach(obs)

    fs.utility.animate(
        system_lj, r=10, jump=2, show=True, frames=100, repeat=False,
        title='LJ Fluid'
    )

    tk = np.mean(obs.result['T_kinetic'])
    tc = np.mean(obs.result['T_configuration'])

    assert np.isclose(tk, tc, rtol=0.1)

if __name__ == "__main__":
    test_lj()
