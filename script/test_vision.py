import numpy as np
import sys
sys.path.append('../src')
sys.path.append('src')
import bd


model = bd.Lavergne_2019_Science
N, dim = 100, 2
nblock, block = 10, 1000
alpha=np.pi / 4
p_act = 0.5

system_abp = model(
    N, dt=0.01, alpha=alpha, p_act=p_act,
    D=1, Pe=2, kT=1, m=1,
    R0=24
)


system_abp.move()
bd.animate(system_abp, r=10, jump=100, box=(-30, 30))
