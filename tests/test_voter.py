import pytest
import fish_sim as fs
import numpy as np
import matplotlib.pyplot as plt


@pytest.mark.parametrize("n", [10, 50])
def test_voter(n, steps=100):
    k = 3
    eta_vals = np.linspace(0, 1, 10)
    magnitisation = np.empty(10)
    for i, eta in enumerate(eta_vals):
        system = fs.cmodel.Voter(n, k, eta)
        system.evolve(steps)
        mag_vals = np.empty(steps)
        for j in range(steps):
            system.move(True)
            mag_vals[j] = np.abs(system.get_polarisation())
        magnitisation[i] = mag_vals.mean()
    return eta_vals, magnitisation


if __name__ == "__main__":
    eta, mag = test_voter(n=10000, steps=100)
    eta_theory = np.linspace(0, 1, 100)
    mag_theory = (3 - 2 / (1 - eta_theory)) ** 0.5
    mag_theory[eta_theory > 1/3] = 0
    plt.scatter(eta, mag)
    plt.plot(eta_theory, mag_theory)
    plt.show()
