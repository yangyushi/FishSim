import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist
import fish_sim as fs


def mc(pos, r, sigma, dx):
    dim, n = pos.shape
    while True:
        idx = np.random.randint(0, n)
        print(idx)
        trial = pos[:, idx] + np.random.uniform(-dx, dx, 2)
        if np.linalg.norm(trial) >= r:
            continue
        mask = np.ones(n, dtype=bool)
        mask[idx] = 0
        dist = np.linalg.norm(pos - trial[:, None], axis=0)
        md = np.min(dist[mask])
        print(md)
        if md > sigma:
            pos[:, idx] = trial
            return pos



def no_overlap(n, r, sigma):
    res = [np.array((0, 0))]
    while len(res) < n:
        pos = np.random.uniform(-r, r, 2)
        if np.linalg.norm(pos) < 2 * r:
            trial = np.array([pos] + res)
            if pdist(trial).min() > sigma:
                print(len(res))
                res = [pos] + res
    return np.array(res).T



@fs.model.Boundary("align_sphere")
class M(fs.model.Vicsek2D):
    def __init__(self, *arg, **kwarg):
        fs.model.Vicsek2D.__init__(self, *arg, **kwarg)
        self.positions = no_overlap(self.n, kwarg['R'] * 1.4, sigma=1.000)
        #self.positions = mc(self.positions, kwarg['R'] * 1.0, sigma=1.000, dx=0.1)
        #self.force_func = fs.force.force_wca
        #self.force_func = fs.force.force_wca_no_overlap
        self.force_func = fs.force.force_lj_no_overlap


N, dim = 61, 2
density = 1.0
alpha = 0.6
eta = 0.1
v0 = 0.02
m = 1e4
r0 = 1.5
R = np.power(3 * N / (density * np.pi * 4), 1/2)
nblock, block = 5, 30
frames = nblock * block
box = (N / density) ** (1 / dim)


system = M(N, eta=eta, r=r0, v0=v0, m=m, alpha=alpha, box=box, R=R)

obs_dyn = fs.utility.Dynamic(block=block, report=True)

system.attach(obs_dyn)


for _ in range(1000):
    system.move()

fs.utility.animate_active_2d(
    system, r=500, jump=1, box=(-R * 1.4, R*1.4), show=True, frames=frames,
    repeat=False, arrow=True, save='disk.gif', fps=30,
)
