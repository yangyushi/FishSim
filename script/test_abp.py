import numpy as np
import sys
sys.path.append('../src')
sys.path.append('src')
import bd


def test_abp():
    model = bd.ABP2DWCAPBC
    model = bd.ABP2D_AS
    model = bd.ABP2D_AS_WCA
    N, dim = 100, 2
    density = 0.5
    dt = 0.0002
    Pe = 100
    R = 40
    nblock, block = 10, 1000
    box = (N / density) ** (1 / dim)

    system_abp = model(
        N, dt=dt, Pe=Pe, box=box,
        D=1, kT=1, m=1, R=R
    )
    bd.animate_active_2d(
        system_abp, r=10,
        jump=100, box=(-R*2, R*2)
    )


if __name__ == "__main__":
    test_abp()
