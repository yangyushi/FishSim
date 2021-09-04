import numpy as np
import sys
sys.path.insert(0, '../lib')
import fish_sim as fs


model = fs.model.Lavergne_2019_Science

N, dim = 100, 2
nblock, block = 10, 1000
alpha=np.pi / 4
p_act = 0.5

system = model(
    N, dt=0.05, alpha=alpha, p_act=p_act,
    D=1, Pe=10, kT=1, m=1,
    R0=100
)
system.move()

fs.utility.animate_active_2d(
    system, r=5, jump=1,
    box=(-120, 120), show=True, frames=100,
    title='vision'
)
