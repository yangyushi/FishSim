import numpy as np
import sys
sys.path.append('../src')
sys.path.append('src')
import bd


def test_abp():
    model = bd.ABP2DWCAPBC
    N, dim = 200, 2
    density = 0.6
    dt = 0.00005
    Pe = 100
    nblock, block = 10, 1000
    box = (N / density) ** (1 / dim)

    system_abp = model(
        N, dt=dt, Pe=Pe, box=box,
        D=1, kT=1, m=1
    )
    bd.animate_active_2d(
        system_abp, r=10,
        jump=20, box=box
    )


if __name__ == "__main__":
    test_abp()
