import numpy as np
import sys
sys.path.append('../src')
sys.path.append('src')
import simulate as sim

@sim.Boundary("align_sphere")
class M(sim.ABP3D): pass


def test_abp():
    model = M
    N, dim = 100, 3
    density = 1
    dt = 0.002
    Pe = 100
    R = 40
    nblock, block = 10, 1000
    box = (N / density) ** (1 / dim)

    system_abp = model(
        N, dt=dt, Pe=Pe, box=box,
        D=1, kT=1, m=1, R=R
    )

    sim.animate(
        system_abp, r=10,
        jump=10, box=(-R*2, R*2),
        show=True
    )


if __name__ == "__main__":
    test_abp()
