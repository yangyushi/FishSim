import numpy as np
import csimulate as cs
from matplotlib import cm
import matplotlib.pyplot as plt

n = 100
dim = 2
box = 1000
eta = 0.2
v0 = 0.17
r = 1.0
run_steps = 1000

positions = np.ones((dim, n)) * box / 2
velocities = np.zeros((dim, n))
velocities[1, :] = -v0

sim = cs.continue_vicsek_2d_pbc(
    n=n, box=box, eta=eta, v0=v0, r=r, run_steps=run_steps,
    positions=positions, velocities=velocities
)

fig = plt.figure(figsize=(2.6, 2.4))
ax = fig.add_subplot()
for i in range(n):
    ax.plot(*sim[:, i, :2].T, color='k', lw=1.6)
    ax.plot(*sim[:, i, :2].T, color=cm.Set3(i % 10), lw=0.8)
ax.set_xticks([])
ax.set_yticks([])
plt.tight_layout()
plt.savefig('flow.svg')
plt.show()
