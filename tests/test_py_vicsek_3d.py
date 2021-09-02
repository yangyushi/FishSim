import numpy as np
import sys
sys.path.insert(0, '../lib')
import fish_sim as fs


@fs.model.Boundary("pbc")
class M1(fs.model.Vicsek3D): pass


@fs.model.Boundary("align_sphere")
class M2(fs.model.Vicsek3D): pass


@fs.model.Boundary("align_fish_bowl")
class M3(fs.model.Vicsek3D): pass


def test_vicsek():
    N, dim = 50, 3
    density = 1
    eta = 0.1
    v0 = 0.05
    r0 = 1
    Pe = 100
    R = np.power(3 * N / (density * np.pi * 4), 1/3)
    nblock, block = 3, 10
    frames = nblock * block
    box = (N / density) ** (1 / dim)


    system = M1(N, eta=eta, r0=r0, v0=v0, box=box, D=1, kT=1, m=1)
    obs_dyn = fs.utility.Dynamic(block=block, report=True)
    system.attach(obs_dyn)
    fs.utility.animate(
        system, r=10, jump=1, box=(0, box), show=True, frames=frames,
        title='Vicsek 3D PBC'
    )

    system = M2(N, eta=eta, r0=r0, v0=v0, R=R, D=1, kT=1, m=1)
    obs_dyn = fs.utility.Dynamic(block=block, report=True)
    system.attach(obs_dyn)
    fs.utility.animate(
        system, r=10, jump=1, box=(-R, R), show=True, frames=frames,
        title='Vicsek 3D Sphere'
    )

    system = M3(N, eta=eta, r0=r0, v0=v0, z_max=1, c=0.25, D=1, kT=1, m=1)
    obs_dyn = fs.utility.Dynamic(block=block, report=True)
    system.attach(obs_dyn)
    fs.utility.animate(
        system, r=10, jump=1, box=(-2, 2), show=True, frames=frames,
        title='Vicsek 3D Fish Bowl'
    )

if __name__ == "__main__":
    test_vicsek()
