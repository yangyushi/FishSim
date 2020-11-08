import numpy as np
import sys
sys.path.append('../src')
sys.path.append('src')
import bd


model = bd.Lavergne_2019_Science
model = bd.L2SWCA
N, dim = 50, 2
nblock, block = 10, 1000
alpha=np.pi / 2
p_act = 0.5

system_abp = model(
    N, dt=0.001, alpha=alpha, p_act=p_act,
    D=1, Pe=10, kT=1, m=1,
    R0=100
)


system_abp.move()
bd.animate_active_2d(system_abp, r=5, jump=50, box=(-500, 500))
