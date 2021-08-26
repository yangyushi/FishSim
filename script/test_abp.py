import numpy as np
import sys
sys.path.append('../src')
sys.path.append('src')
import fish_sim as sim


@sim.Boundary("PBC")
class M1(sim.ABP2D): pass

def test_abp():
    model = sim.ABP2DWCAPBC
    #model = sim.ABP2D_AS
    #model = sim.ABP2D_AS_WCA
    #model = M1
    N, dim = 100, 2
    density = 0.5
    dt = 0.0002
    Pe = 50
    R = 50
    nblock, block = 10, 1000
    box = (N / density) ** (1 / dim)

    system_abp = model(
        N, dt=dt, Pe=Pe, box=box,
        D=1, kT=1, m=1, R=R
    )
    sim.animate_active_2d(
        system_abp, r=10,
        jump=1, box=(0, box),
        show=True
    )


if __name__ == "__main__":
    test_abp()
