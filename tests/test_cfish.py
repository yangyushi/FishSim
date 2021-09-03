import sys
sys.path.append("../lib")
import numpy as np
import fish_sim as fs


def test_vicsek_3d():
    traj = fs.cfish_sim.vicsek_3d_pbc(
        n=20, r=1.0, box=10, eta=0.5, v0=0.01, pre_steps=100, run_steps=100
    )
    assert traj.shape == (100, 20, 6)
    assert not np.isnan(traj).any()

def test_vicsek_2d():
    traj = fs.cfish_sim.vicsek_2d_pbc(
        n=20, r=1.0, box=10, eta=0.5, v0=0.01, pre_steps=100, run_steps=100
    )
    assert traj.shape == (100, 20, 4)
    assert not np.isnan(traj).any()
