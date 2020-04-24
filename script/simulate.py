import os
import sys
sys.path.append('../lib')
import time
import numpy as np
import matplotlib.pyplot as plt
import csimulate


N = 200
density = 1.0
box = (N / density) ** (1/3)
spd = 0.1
r = 1.0
noises = np.arange(0.0, 1.2, 0.2)

orders = []

t0 = time.time()
for reduced_noise in noises:
    noise = reduced_noise * np.sqrt(density)
    result = csimulate.vicsek_3d_pbc(N, box, noise, spd, r, 2000, 2000)  # (T, n, 6)
    order = np.sqrt((result[:, :, 3:].mean(1) ** 2).sum(-1)).mean()
    orders.append(order / spd)
t_cpp = time.time() - t0

# plotting results from a Python code
result_py = np.load('result_py.npy')
plt.scatter(*result_py, color='tomato', s=64, facecolor='none', label=f'Python: ~280s')

# plotting C++ results
plt.scatter(noises, orders, color='teal', marker='+', label=f'C++: {t_cpp:.2f}s', s=72)
plt.xlabel('Noise', fontsize=14)
plt.ylabel('Dynamical Order', fontsize=14)
plt.legend(fontsize=14)
plt.tight_layout()
plt.savefig('result.pdf')
plt.show()
