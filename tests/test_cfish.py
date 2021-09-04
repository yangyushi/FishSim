import sys
sys.path.insert(0, "../lib")
import numpy as np
from time import time

import fish_sim as fs


def test_vicsek_3d():
    traj = fs.cmodel.vicsek_3d_pbc(
        n=20, r=1.0, box=10, eta=0.5, v0=0.01, pre_steps=100, run_steps=100
    )
    assert traj.shape == (100, 20, 6)
    assert not np.isnan(traj).any()

def test_vicsek_2d():
    traj = fs.cmodel.vicsek_2d_pbc(
        n=20, r=1.0, box=10, eta=0.5, v0=0.01, pre_steps=100, run_steps=100
    )
    assert traj.shape == (100, 20, 4)
    assert not np.isnan(traj).any()

def test_cpp_optimise():
    N = 100
    v0 = 0.05
    r0 = 1
    density = 1
    eta = 0.5
    box = np.power(N, 1.0/3.0)


    @fs.model.Boundary("pbc")
    class V3PBC(fs.model.Vicsek3D): pass
    t_py = time()
    system = V3PBC(N, eta, v0, r0, box=box)
    for _ in range(100):
        system.move()
    t_py = time() - t_py
    print(f"Python: {t_py:.4f}s")

    t_cpp = time()
    system = fs.cmodel.Vicsek3DPBC(n=N, r=r0, eta=eta, v0=v0, box=box)
    for _ in range(100):
        system.move()
    t_cpp = time() - t_cpp
    print(f"C++: {t_cpp:.4f}s")

    assert (t_cpp * 10) < t_py, "10x speed up not achieved, check the optimisation setup"



if __name__ == "__main__":
    test_vicsek_2d()
    test_vicsek_3d()
    test_cpp_optimise()
