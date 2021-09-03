import fish_sim as fs


system = fs.cmodel.Vicsek3D(
    n=100, r=1.0, eta=0.1, v0=0.05
)
fs.utility.animate(
    system, show=True, frames=25, repeat=False, r=1, box=(-50, 50)
)

n = 100
system = fs.cmodel.Vicsek3DPBC(
    n=n, r=1, eta=0.2, box=5, v0=0.005
)
fs.utility.animate(system, show=True, frames=25, repeat=False, r=5)
