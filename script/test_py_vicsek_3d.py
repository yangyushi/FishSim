import numpy as np
import sys
sys.path.append('../src')
sys.path.append('src')
import bd

@bd.Boundary("align_sphere")
class M(bd.Vicsek3D): pass


def test_vicsek():
    model = M
    N, dim = 100, 3
    density = 1
    eta = 0.1
    v0 = 0.05
    r0 = 1
    Pe = 100
    R = np.power(3 * N / (density * np.pi * 4), 1/3)
    nblock, block = 10, 1000
    box = (N / density) ** (1 / dim)

    system_abp = model(
        N, eta=eta, r0=r0, v0=v0, box=box,
        D=1, kT=1, m=1, R=R
    )

    bd.animate(
        system_abp, r=10,
        jump=1, box=(-R*2, R*2),
        show=True
    )


if __name__ == "__main__":
    test_vicsek()
