import numpy as np
import sys
sys.path.insert(0, '../lib')
import fish_sim as fs


def test_network():
    n = 50
    k = 5
    eta_vals = np.linspace(0, 1, 20)

    pol_vals = []
    for eta in eta_vals:
        system = fs.cmodel.Network3D(n, k, eta)
        for _ in range(100):
            system.move()

        pol = []
        for _ in range(100):
            system.move()
            pol.append(system.get_polarisation())

        pol_vals.append(np.mean(pol))
    assert np.any(np.diff(pol_vals[::-1]) > 0)
    return eta_vals, pol_vals

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    mpl.style.use('yushi')

    result = test_network()
    plt.plot(*result, marker='o')
    plt.savefig('network-fast.pdf')
