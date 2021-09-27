import numpy as np
import sys
sys.path.insert(0, '../lib')
import fish_sim as fs


def test_couzin(repeat=False):
    args = {
        'n' : 50,
        'rr' : 1,
        'ro' : 2,
        'ra' : 20,
        'perception' : 90 / 180 * np.pi,
        'noise':0.1,
        'speed': 3,
        'turn_rate':40 / 180 * np.pi,
        'dt':0.01
    }
    system = fs.cmodel.Couzin3D(**args)
    fs.utility.animate(
        system, show=True, repeat=repeat, r=4,
        frames=30, box=(-25, 25), fps=30,
        figsize=(5, 5), arrow=1
    )


def test_couzin_tank(repeat=False):
    args = {
        'n' : 50,
        'rr' : 0 / 1e3,
        'ro' : 0 / 1e3,
        'ra' : 0 / 1e3,
        'perception' : 90 / 180 * np.pi,
        'noise':0.5,
        'speed': 0.1,
        'turn_rate':40 / 180 * np.pi,
        'dt':0.005,
        'c': 0.743,
        'h': 0.5,
        'kw': 0.01,
        'align_cap': False,
        'align_base': True,
        'g': 1,
    }
    system = fs.cmodel.CouzinTank3D(**args)
    fs.utility.animate(
        system, show=True, repeat=repeat, r=4,
        frames=30, fps=30,
        figsize=(5, 5), arrow=0.02
    )



if __name__ == "__main__":
    test_couzin_tank(repeat=True)
    test_couzin(repeat=True)
