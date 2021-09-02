import numpy as np
import sys
sys.path.insert(0, '../lib')
import fish_sim as fs


@fs.model.Boundary("PBC")
class M1(fs.model.ABP2D): pass

def test_abp():
    model = fs.model.ABP2DWCAPBC
    N, dim = 200, 2
    density = 0.5
    dt = 0.0002
    Pe = 10
    R = 50
    nblock, block = 3, 100
    box = (N / density) ** (1 / dim)

    system_abp = model(
        N, dt=dt, Pe=Pe, box=box,
        D=1, kT=1, m=1, R=R
    )
    obs_dyn = fs.utility.Dynamic(block=block)
    system_abp.attach(obs_dyn)

    fs.utility.animate_active_2d(
        system_abp, r=10, jump=5, box=(0, box),
        show=True, frames=100, title="2D ABP PBC"
    )


if __name__ == "__main__":
    test_abp()
