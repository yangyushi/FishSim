import pytest
import fish_sim as fs
import numpy as np


@pytest.mark.parametrize("n", [10, 50])
def test_voter(n, steps=100):
    k = 3
    eta_vals = np.linspace(0.01, 0.99, 10)
    magnitisation = np.empty(10)
    for i, eta in enumerate(eta_vals):
        system = fs.cmodel.Voter(n, k, eta)
        system.evolve(steps)
        mag_vals = np.empty(steps)
        for j in range(steps):
            system.move(True)
            mag_vals[j] = np.abs(system.get_polarisation())
        magnitisation[i] = mag_vals.mean()

    tmp = 3 - 2 / (1 - eta_vals)
    tmp[tmp < 0] = 0
    theory = (tmp) ** 0.5
    assert (theory - magnitisation).mean() < 1e-2


if __name__ == "__main__":
    test_voter(n=2000, steps=100)
