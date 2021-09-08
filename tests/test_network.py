import numpy as np
import sys
sys.path.insert(0, '../lib')
import fish_sim as fs


def test_network():
    n = 10
    k = 5
    eta_vals = np.linspace(0, 0.5, 10)

    pol_vals = []
    for eta in eta_vals:
        system = fs.cmodel.Network3D(n, k, eta)
        for _ in range(1000):
            system.move()

        pol = []
        for _ in range(1000):
            system.move()
            pol.append(system.get_polarisation())

        pol_vals.append(np.mean(pol))
    assert np.any(np.diff(pol_vals[::-1]) > 0)

if __name__ == "__main__":
    test_network()
