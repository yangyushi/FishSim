import numpy as np
import sys
sys.path.append('../src')
sys.path.append('src')
import fish_sim as sim


@sim.Boundary("pbc")
class M1(sim.Vicsek3D): pass


@sim.Boundary("align_sphere")
class M2(sim.Vicsek3D): pass


@sim.Boundary("align_fish_bowl")
class M3(sim.Vicsek3D): pass


def test_vicsek():
    N, dim = 50, 3
    density = 1
    eta = 0.2
    v0 = 0.05
    r0 = 1
    Pe = 100
    R = np.power(3 * N / (density * np.pi * 4), 1/3)
    nblock, block = 10, 1000
    box = (N / density) ** (1 / dim)


    system = M1(N, eta=eta, r0=r0, v0=v0, box=box, D=1, kT=1, m=1)
    sim.animate( system, r=10, jump=1, box=(0, box), show=True)

    system = M2(N, eta=eta, r0=r0, v0=v0, R=R, D=1, kT=1, m=1)
    sim.animate(system, r=10, jump=1, box=(-R, R), show=True)

    system = M3(N, eta=eta, r0=r0, v0=v0, z_max=1, c=0.25, D=1, kT=1, m=1)
    sim.animate(system, r=10, jump=1, box=(-2, 2), show=True)

if __name__ == "__main__":
    test_vicsek()
