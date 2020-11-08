import numpy as np
import sys
sys.path.append('../src')
sys.path.append('src')
import bd


def test_abp():
    model = bd.ABP2D
    model = bd.ABP2DWCAPBC
    N, dim = 512, 2
    density = 0.5
    nblock, block = 10, 1000
    box = (N / density) ** (1 / dim)

    system_abp = model(
        N, dim, box=box, dt=0.0001,
        D=1, Pe=50, kT=1, m=1
    )


    bd.animate(system_abp, r=10, jump=10)


if __name__ == "__main__":
    test_abp()
