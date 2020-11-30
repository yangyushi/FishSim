import numpy as np
import sys
sys.path.append('../src')
sys.path.append('src')
import simulate as sim
import fish_corr as fc

@sim.Boundary("align_half_sphere")
class M1(sim.Vicsek3D): pass

@sim.Boundary("pbc")
class M2(sim.Vicsek3D): pass


def test_vicsek():
    N, dim = 50, 3
    density = 0.5
    eta = 0.2
    v0 = 0.05
    r0 = 1
    Pe = 100
    R = np.power(3 * N / (density * np.pi * 4), 1/3)
    nblock, block = 10, 1000
    box = (N / density) ** (1 / dim)

    system = M1(
        N, eta=eta, r0=r0, v0=v0, box=box,
        D=1, kT=1, m=1, R=R
    )

    dump_xyz = sim.DumpXYZ(frequency=10, filename='vicsek-3d-ahs')
    dump_mod = sim.DumpModel(frequency=10, filename='vicsek-3d-ahs')

    system.attach(dump_xyz)
    system.attach(dump_mod)

    sim.animate(
        system, r=10,
        jump=1, box=(-R*1, R*1),
        show=True,
       # save='vicsek_3d_aligned_sphere.gif',
    )
    dump_mod.dump(fc.SimMovie)

    """
    system = M2(
        N, eta=eta, r0=r0, v0=v0, box=box,
        D=1, kT=1, m=1, R=R
    )
    sim.animate(
        system, r=10,
        jump=1, box=(0, box),
        show=False, save='vicsek_3d_pbc.gif'
    )
    """


if __name__ == "__main__":
    test_vicsek()
