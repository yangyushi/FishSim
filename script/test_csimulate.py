#!/usr/bin/env python3
import os
import sys
sys.path.append('../lib')
import time
import numpy as np
import matplotlib.pyplot as plt
import fish_sim as fs


N = 500
density = 5.0
box = (N / density) ** (1/3)
spd = 0.1
r = 1.0
noises = np.arange(0.0, 1.0, 0.05)

orders = []

t_cpp = 0
for noise in noises:
    t0 = time.time()
    result = fs.csimulate.vicsek_3d_pbc(N, box, noise, spd, r, 1000, 1000)  # (T, n, 6)
    t_cpp += time.time() - t0
    order = np.sqrt((result[:, :, 3:].mean(1) ** 2).sum(-1)).mean()
    orders.append(order / spd)

# plotting results from a Python code
result_py = np.load('result_py.npy')
plt.scatter(*result_py, color='tomato', s=64, facecolor='none', label=f'Python: ~5000s')

# plotting C++ results
plt.scatter(noises, orders, color='teal', marker='+', label=f'C++: {t_cpp:.2f}s', s=72)
plt.xlabel('Noise', fontsize=14)
plt.ylabel('Dynamical Order', fontsize=14)
plt.legend(fontsize=14)
plt.tight_layout()
plt.savefig('result.pdf')
plt.show()
