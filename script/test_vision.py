import numpy as np
import sys
sys.path.append('../src')
sys.path.append('src')
import simulate as sim


model = sim.Lavergne_2019_Science

N, dim = 50, 2
nblock, block = 10, 1000
alpha=np.pi / 4
p_act = 0.5

system_abp = model(
    N, dt=0.005, alpha=alpha, p_act=p_act,
    D=1, Pe=10, kT=1, m=1,
    R0=100
)


system_abp.move()
sim.animate_active_2d(
    system_abp, r=5, jump=100,
    box=(-120, 120), save='vision.gif', show=True, frames=100
)

