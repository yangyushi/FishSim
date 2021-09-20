import numpy as np
import sys
sys.path.insert(0, '../lib')
import fish_sim as fs


def test_couzin():
    args = {
        'n' : 100,
        'rr' : 1,
        'ro' : 2,
        'ra' : 20,
        'perception' : 90 / 180 * np.pi,
        'noise':0.1,
        'speed': 3,
        'turning_rate':40 / 180 * np.pi,
        'dt':0.005
    }
    system = fs.cmodel.Couzin3D(**args)
    fs.utility.animate(
        system, show=True, repeat=False, r=4,
        frames=30, box=(-25, 25), fps=30,
        figsize=(5, 5)
    )

if __name__ == "__main__":
    test_couzin()
