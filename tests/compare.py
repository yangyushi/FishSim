import numpy as np
import fish_sim as fs
from time import time


N = 100
v0 = 0.05
r0 = 1
density = 1
eta = 0.5
box = np.power(N, 1.0/3.0)


@fs.model.Boundary("pbc")
class V3PBC(fs.model.Vicsek3D): pass


t_py = time()
system = V3PBC(N, eta, v0, r0, box=box)
for _ in range(100):
    system.move()
t_py = time() - t_py
print(f"Python: {t_py:.4f}s")

t_cpp = time()
result = fs.cfish_sim.vicsek_3d_pbc(
    n=N, box=box, eta=eta, v0=v0, r=r0,
    pre_steps=50, run_steps=50, jump=1, use_nl=True
)
t_cpp = time() - t_cpp
print(f"C++: {t_cpp:.4f}s")
