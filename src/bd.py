from abc import ABCMeta, abstractmethod
import matplotlib.pyplot as plt
from itertools import product
from scipy.spatial.distance import pdist, cdist, squareform
import matplotlib.animation as animation
import numpy as np
from numba import njit
from force import force_lj, force_wca
from scipy.special import gamma as gamma_func
from noise_3d import add_vicsek_noise_3d


def animate(system, r=100, jump=100, box=None,
            save='', fps=60, show=False, frames=100):
    fig = plt.figure(figsize=(5, 5), tight_layout=True)
    if system.dim == 2:
        ax = fig.add_subplot()
        ax.set_xticks([])
        ax.set_yticks([])
        if not isinstance(box, type(None)):
            if type(box) == tuple:
                ax.set_xlim(*box)
                ax.set_ylim(*box)
            else:
                ax.set_xlim(0, box)
                ax.set_ylim(0, box)
        elif hasattr(system, 'box'):
            ax.set_xlim(0, system.box)
            ax.set_ylim(0, system.box)
    elif system.dim == 3:
        ax = fig.add_subplot(projection='3d')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([])
        if not isinstance(box, type(None)):
            if type(box) == tuple:
                ax.set_xlim(*box)
                ax.set_ylim(*box)
                ax.set_zlim(*box)
            else:
                ax.set_xlim(0, box)
                ax.set_ylim(0, box)
                ax.set_zlim(0, box)
        elif hasattr(system, 'box'):
            ax.set_xlim(0, system.box)
            ax.set_ylim(0, system.box)
            ax.set_zlim(0, system.box)
    else:
        return NotImplementedError("Only 2D and 3D systems are Supported")

    def update(num, system, scatter):
        for _ in range(jump):
            system.move()
        if system.dim == 3:
            scatter.set_data(system.r.T[:2])
            scatter.set_3d_properties(system.r.T[2])
        else:
            scatter.set_data(system.r.T)
        return scatter

    scatter = ax.plot(
        *system.r.T, color='teal', mfc='w', ls='None', marker='o',
        markersize=r
    )[0]
    ani = animation.FuncAnimation(
        fig, update, frames=frames, fargs=(system, scatter), interval=1
    )
    if show:
        plt.show()
    if save:
        ani.save(save, writer='imagemagick', fps=fps)


def animate_active_2d(
        system, r=100, jump=100, box=None,
        save='', fps=60, show=True, frames=100):
    fig = plt.figure(figsize=(5, 5), tight_layout=True)
    ax = fig.add_subplot()
    ax.set_xticks([])
    ax.set_yticks([])
    if not isinstance(box, type(None)):
        if type(box) == tuple:
            ax.set_xlim(*box)
            ax.set_ylim(*box)
        else:
            ax.set_xlim(0, box)
            ax.set_ylim(0, box)
    elif hasattr(system, 'box'):
        ax.set_xlim(0, system.box)
        ax.set_ylim(0, system.box)
    if system.dim != 2:
        return NotImplementedError("Only 2D system is Supported")

    def update(num, system, scatter):
        for _ in range(jump):
            system.move()
            scatter.set_data(system.r.T)
            quiver.set_offsets(system.r)
            quiver.set_UVC(
                np.cos(system.phi),
                np.sin(system.phi),
            )
        return scatter

    scatter = ax.plot(
        *system.r.T, color='teal', mfc='w', ls='None', marker='o',
        markersize=r
    )[0]
    quiver = ax.quiver(
        *system.r.T, np.cos(system.phi), np.sin(system.phi),
        pivot='mid', units='width', color='teal', zorder=5, scale=250/r,
    )
    ani = animation.FuncAnimation(
        fig, update, frames, fargs=(system, scatter), interval=1
    )
    if show:
        plt.show()
    if save:
        ani.save(save, writer='imagemagick', fps=fps)


class Observer():
    def __init__(self, block):
        self.block = block
        self.count = 0

    def update(self, system):
        if self.count == self.block - 1:
            self.collect(system)
            self.aggregate()
            self.count = 0
            system.interest = {}
        else:
            self.collect(system)
            self.count += 1

    def aggregate(self, system):
        """
        called at the end of each block
        """
        pass

    def collect(self, system):
        """
        called at the end of each running step
        """
        pass


class Thermodynamic(Observer):
    """
    Calculate the block average of thermodynamic quantities.
        The required quantities of interest are,

         - kinetic        : the total kinetic energy
         - potential      : the total potential energy (cut & shifted)
         - potential_full : the total potential energy
         - virial         : the total virial
         - laplacian      : the total Laplacian
         - force_sq       : the total squared force
    """
    def __init__(self, block):
        Observer.__init__(self, block)
        self.Ek = 0
        names = ('E/N cut&shifted', 'P cut&shifted', 'T Kinetic', 'T Config')
        print('|'.join([f'{n.split(" ")[0]: ^20}' for n in names]))
        print('|'.join([f'{n.split(" ")[1]: ^20}' for n in names]))
        print('-' * (21 * len(names)))

    def aggregate(self, system):
        e_kin = np.mean(system.interest['kinetic']) / system.N
        e_tot = e_kin + np.mean(system.interest['potential']) / system.N
        press = system.kT * system.density + np.mean(system.interest['virial']) / system.volume
        t_kin = 2 * e_kin / system.dim
        t_con = np.mean(
            np.array(system.interest['force_sq']) / np.array(system.interest['laplacian'])
        )
        quantities = (e_tot, press, t_kin, t_con)
        print('|'.join([f'{val: ^20.4f}' for val in quantities]))


class BD():
    """
    The basic Brownian Dynamic simulation of ideal gas
    """
    def __init__(self, N, dim, dt, gamma=None, kT=None, m=None, D=None, **kwargs):
        self.N, self.dim, self.dt = N, dim, dt
        self.gamma = gamma
        self.kT = kT
        self.m = m
        self.D = D
        self.update_parameters()
        self.r = np.zeros((N, dim))  # positions
        self.f = np.zeros_like(self.r)
        self.interest = {}  # collect quantity of interest
        self.__observers = []  # ovserver design pattern
        sigma = np.sqrt(self.kT / (self.m * self.N))
        self.v = np.random.normal(0, sigma, (self.N, self.dim))

    def attach(self, observer):
        self.interest = {}
        if observer not in self.__observers:
            self.__observers.append(observer)

    def detach(self, observer):
        self.interest = {}
        try:
            self.__observers.remove(observer)
        except:
            pass

    def notify(self):
        for observer in self.__observers:
            observer.update(self)

    def update_parameters(self):
        par = [self.gamma, self.kT, self.m, self.D]
        n_unknown = sum([isinstance(x, type(None)) for x in par])
        if n_unknown == 1:
            if isinstance(self.D, type(None)):
                self.D = self.kT / self.m / self.gamma
            elif isinstance(self.gamma, type(None)):
                self.gamma = self.kT / self.m / self.D
            elif isinstance(self.kT, type(None)):
                self.kT = self.m * self.gamma * self.D
            elif isinstance(self.m, type(None)):
                self.m = self.kT / self.D / self.gamma

    def get_pair_shift(self):
        """
        The content in the result is rij = ri - rj, shape (dim, N, N)

        ..code-block::

                 r11      | r12 = r1 - r2 | ... | r13 = r1 - r3 |
            r21 = r2 - r1 |      r22      | ... | r23 = r2 - r3 |
                  .       |       .       |  .  |       .       |
                  .       |       .       |  .  |       .       |
                  .       |       .       |  .  |       .       |
            rn1 = rn - r1 | rn2 = r2 - rn | ... |      rnn      |

        To get force, sum over the row. (axis = 0)

        To get sfhits, take columns. (vectors from i to all others = rij[:, i])
        """
        return self.r.T[:, :, np.newaxis] - self.r.T[:, np.newaxis, :]

    def fix_boundary(self): pass

    def force_func(self, shift): return np.zeros((self.N, self.dim)), {}

    def get_force(self):
        rij = self.get_pair_shift()
        self.f, interest = self.force_func(rij)
        if self.m:
            interest['kinetic'] = 0.5 * self.m * np.sum(self.v ** 2)
        else:
            interest['kinetic'] = 0.5 * np.sum(self.v ** 2)
        for key in interest:
            if key in self.interest:
                self.interest[key].append(interest[key])
            else:
                self.interest.update({key : [ interest[key] ] })

    def move(self):
        """
        use baoab algorithm to integrate the SDE

        Args:
            gamma (float): the damping constant
            kT (float): the temperature times Boltzmann constant
            m (float): the mass
        """
        self.v += 0.5 * self.dt * self.f / self.m  # b propagator

        self.r += 0.5 * self.dt * self.v  # a propagator
        self.fix_boundary()

        # --- o propagator ---
        x = self.gamma * self.dt
        c = np.sqrt(
            1 - np.exp(- 2 * x) if x > 0.0001 else \
                np.polyval([-2/3,4/3,-2.0,2.0,0.0], x)  # Tayler Expansion
        )
        self.v = self.v * np.exp(-x) + c * np.sqrt(self.m * self.kT) * np.random.randn(self.N, self.dim)
        # -------

        self.r += 0.5 * self.dt * self.v  # a propagator
        self.fix_boundary()

        self.get_force()

        self.v += 0.5 * self.dt * self.f / self.m  # b propagator

        self.notify()

    def move_overdamp(self):
        self.get_force()
        self.r += self.D / self.kT * self.f * self.dt + \
                  np.sqrt(2 * self.D * self.dt) * np.random.randn(self.N, self.dim)
        self.fix_boundary()
        self.notify()


def Boundary(condition):
    """
    A class decorator with parameter to enforce periodic boundary for a
        MD/MC simulation. This decorator can be extended to cover other
        geometries.

    Args:
        condition (str): the type of the boundary condition
            corrently pbc is the only option. If it is not the pbc then
            the simulation is running in a free space
    """
    def pbc_decorator(system):
        """
        A class decorator to enforce periodic boundary for a MD/MC simulation

        Args:
            system (class): a class inhereting the abstract class [Dynamic] for
                a simulation.
        """
        class SystemPBC(system):
            def __init__(self, *args, **kwargs):
                self.box = kwargs['box']
                system.__init__(self, *args, **kwargs)
                try:
                    self.box = float(self.box)
                except (ValueError, TypeError) as err:
                    raise ValueError("the box should be a numerical value")
                self.volume = self.box ** self.dim
                self.density = self.N / self.volume
                repeat = np.ceil(np.power(self.N, 1 / self.dim)).astype(int)
                lattice = np.array(list(product(*[range(repeat)] * self.dim)))
                self.r = lattice[:self.N] / repeat * kwargs['box']
                self.get_force()

            def fix_boundary(self):
                self.r %= self.box

            def get_pair_shift(self):
                pos_in_box = self.r.T / self.box
                pair_shift = pos_in_box[:, :, np.newaxis] - pos_in_box[:, np.newaxis, :]
                pair_shift = (pair_shift - np.rint(pair_shift)) * self.box
                return pair_shift

        return SystemPBC

    def align_sphere(system):
        class SystemAS(system):
            def __init__(self, *args, **kwargs):
                self.R = kwargs['R']  # radius of the bounding sphere
                system.__init__(self, *args, **kwargs)
                try:
                    self.R = float(self.R)
                except (ValueError, TypeError) as err:
                    raise ValueError("the radius should be a numerical value")
                self.volume = np.pi**(self.dim / 2) /\
                              gamma_func(self.dim/2 + 1) * self.R**self.dim
                self.density = self.N / self.volume
                p = np.random.randn(self.dim, self.N)  # positions
                r = np.linalg.norm(p, axis=0)  # radii
                l = np.random.uniform(0, self.R**2, self.N) ** (1 / self.dim)
                self.r = (p / r * l).T  # uniform inside n-sphere
                self.get_force()

            def fix_boundary(self):
                """
                fish that are outside of the circle should not move further
                """
                v0 = self.v.copy()
                is_outside = np.linalg.norm(self.r, axis=1) >= self.R
                if not np.any(is_outside):
                    return
                v_orient = self.v[is_outside]  # n, 2
                r_orient = self.r[is_outside]  # n, 2
                v_mag = np.linalg.norm(v_orient, axis=1)  # n,
                r_mag = np.linalg.norm(r_orient, axis=1)  # n,
                v_orient /= v_mag[:, np.newaxis]
                r_orient /= r_mag[:, np.newaxis]
                angles = np.arccos(np.einsum('ij,ij->i', v_orient, r_orient))  # n,
                signs = np.sign(np.cross(v_orient, r_orient))
                # if greater than pi/2, no change, otherwise change to pi/2
                should_fix = np.abs(angles) < np.pi/2  # n
                if not np.any(should_fix):
                    return
                if self.dim == 2:
                    r_angle = np.arctan2(r_orient[:, 1], r_orient[:, 0])  # n,
                    v_orient[should_fix, 0] = np.cos(
                        r_angle[should_fix] - signs[should_fix] * np.pi/2
                    ) * v_mag[should_fix]
                    v_orient[should_fix, 1] = np.sin(
                        r_angle[should_fix] - signs[should_fix] * np.pi/2
                    ) * v_mag[should_fix]
                    self.v[is_outside, :] = v_orient
                    self.phi[is_outside] = np.arctan2(v_orient[:, 1], v_orient[:, 0])
                elif self.dim == 3:
                    e_align = np.cross(
                        r_orient[should_fix],
                        np.cross(v_orient[should_fix], r_orient[should_fix]),
                    )
                    e_align = e_align / np.linalg.norm(e_align, axis=1)[:, np.newaxis]
                    to_fix = is_outside.copy()
                    to_fix[to_fix] *= should_fix   # to fix = outside & angle < np.pi/2
                    self.o[to_fix] = e_align
                    self.v[is_outside] = self.o[is_outside] * v_mag[:, np.newaxis]
                else:
                    raise ValueError("Only 2D and 3D systems are supported")

        return SystemAS

    if condition == "pbc":
        return pbc_decorator
    elif condition == "align_sphere":
        return align_sphere
    else:
        return lambda x: x


@Boundary("pbc")
class BDLJPBC(BD):
    """
    Integrating the overdamped Langevin equation, performing
        the Brownian dynamic simulation.

    Attributes:
        r (np.ndarray): the positions of all the particles in the system.
        m (float): the mass of the particles.
        D (float): the diffusion coefficient.
        kT (float): the (initial) temperature of the system.
        T (float): the temperature of the system.
        kB (float): the boltzman constant.
    """
    def __init__(self, N, dim, dt, **kwargs):
        """
        Set initial positions on a simple cubic lattice, and
            the initial velocities were randomly sampled according
            to the Maxwell-Boltzmann distribution.

        Args:
            N (int): the number of particles.
            dim (int): the dimension of the system.
            box (float): the size of the cubic periodic box.
        """
        BD.__init__(self, N, dim, dt, **kwargs)

        r_cut = 2.5  # usual constant for LJ system
        if self.box / 2 <= r_cut:  # naïve pbc calculation breaks otherwise
            raise ValueError("The box size must be grater than 5")

        self.force_func = force_lj


class Vicsek3D(BD):
    """
    Vicsek Model Simulation in 3D
    """
    def __init__(self, N, eta, v0, r0, **kwargs):
        BD.__init__(self, N, dim=3, dt=1, **kwargs)
        self.v0 = v0  # speed
        self.r0 = r0  # interaction range
        self.eta = eta
        self.r = np.random.randn(self.N, self.dim)
        self.o = np.random.randn(self.N, self.dim)  # orientation
        self.o = self.o / np.linalg.norm(self.o, axis=1)[:, np.newaxis]
        self.v = self.o * self.v0

    def move(self):
        self.get_force()
        # active
        self.r += self.v * self.dt

        rij = self.get_pair_shift()
        dij = np.linalg.norm(rij, axis=0)  # distances, shape (N, N)
        aij = dij < self.r0  # adjacency matrix (N, N)
        nni = np.sum(aij, axis=0)  # number of neighbours for each particle

        vii = self.v.T[:, np.newaxis, :] * np.ones((3, self.N, 1))  # (3, N, N)
        v_mean = np.sum(vii * aij[np.newaxis, :, :], axis=2) / nni  # (3, N)
        v_xy = np.linalg.norm(v_mean[:2], axis=0)  # (N,)
        azi = np.arctan2(v_mean[1], v_mean[0]) # (N,)
        ele = np.arctan(v_mean[2] / v_xy)  # (N,)
        azi, ele = add_vicsek_noise_3d(azi, ele, eta=self.eta)
        e_xy = np.cos(ele)
        self.o[:, 2] = np.sin(ele)
        self.o[:, 1] = e_xy * np.sin(azi)
        self.o[:, 0] = e_xy * np.cos(azi)
        self.v = self.o * self.v0

        self.fix_boundary()
        self.notify()

    def move_overdamp(self): self.move()


class ABP2D(BD):
    """
    Pe = 3 * v0 * τr / σ = 3 v0 / Dr = v0 / D  (σ = 1)
    Dr = 3 * Dt / σ^2    = 3 D                 (σ = 1)
    D  = kT / (m * γ)                          (σ = 1)
    v0 = Pe * D

    see wysocki-2014-EPL (3D simulation)
    """
    def __init__(self, N, dt, Pe, **kwargs):
        BD.__init__(self, N, 2, dt, **kwargs)
        if not self.D:
            raise KeyError("Missing the Diffusion constant of the system")
        self.Pe = Pe
        self.v0 = self.Pe * self.D
        self.Dr = 3 * self.D
        self.r = np.random.randn(self.N, self.dim)
        self.o = np.random.randn(self.N, self.dim)  # orientation
        self.o = self.o / np.linalg.norm(self.o, axis=1)[:, np.newaxis]
        self.v = self.o * self.v0
        self.phi = np.random.uniform(0, 2 * np.pi, self.N)

    def move_overdamp(self):
        """
        dx = v0 cosφ + √2D ξ_x
        dy = v0 sinφ + √2D ξ_x
        dφ = √2Dr ξ_φ
        """
        self.get_force()
        # active
        self.r += self.v * self.dt
        # Brownian
        self.r += self.D / self.kT * self.f * self.dt + \
                  np.sqrt(2 * self.D * self.dt) * np.random.randn(self.N, self.dim)
        # Re-orient
        self.phi += np.sqrt(2 * self.Dr * self.dt) * np.random.randn(self.N)
        self.phi = self.phi % (2 * np.pi)
        self.o = np.array((np.cos(self.phi), np.sin(self.phi))).T  # (n, 2)
        self.v = self.o * self.v0
        self.fix_boundary()
        self.notify()

    def move(self):
        self.move_overdamp()


class ABP3D(BD):
    """
    Pe = 3 * v0 * τr / σ = 3 v0 / Dr = v0 / D  (σ = 1)
    Dr = 3 * Dt / σ^2    = 3 D                 (σ = 1)
    D  = kT / (m * γ)                          (σ = 1)
    v0 = Pe * D

    see wysocki-2014-EPL (3D simulation)
    """
    def __init__(self, N, dt, Pe, **kwargs):
        BD.__init__(self, N, dim=3, dt=dt, **kwargs)
        if not self.D:
            raise KeyError("Missing the Diffusion constant of the system")
        self.Pe = Pe
        self.v0 = self.Pe * self.D
        self.Dr = 3 * self.D
        self.r = np.random.randn(self.N, self.dim)
        self.o = np.random.randn(self.N, self.dim)  # orientation
        self.o /= np.linalg.norm(self.o, axis=1)[:, np.newaxis]
        self.v = self.o * self.v0

    def rot_diffuse(self):
        """
        Perform the rotational diffusion
        (https://github.com/FTurci/active-lammps, I copied Francesco's code)
        (Winkler et al. Soft Matter, 2015, 11, 6680)
        """
        ox, oy, oz = self.o.T
        cos_theta = oz
        sin_theta = np.sqrt(1 - oz**2)
        phi = np.arctan2(oy, ox)
        cos_phi = np.cos(phi)
        sin_phi = np.sin(phi)

        e_theta = np.empty((3, self.N))
        e_phi = np.empty((3, self.N))
        e_theta[0] =  cos_phi * cos_theta
        e_theta[1] =  sin_phi * cos_theta
        e_theta[2] =           -sin_theta
        e_phi[0]   = -sin_phi
        e_phi[1]   =  cos_phi
        e_phi[2]   =  0.0

        r1, r2 = np.random.randn(2, self.N) * np.sqrt(2 * self.Dr)
        dt_wiener = np.sqrt(self.dt)
        for i in range(3):
            self.o[:, i] += e_theta[i] * dt_wiener * r1 +\
                            e_phi[i]   * dt_wiener * r2 -\
                        2 * self.Dr * self.dt * self.o[:, i]
        self.o /= np.linalg.norm(self.o, axis=1)[:, np.newaxis]

    def move_overdamp(self):
        """
        dx = v0 cos cosφ + √2D ξ_x
        dy = v0 sinφ + √2D ξ_x
        dz = v0 sinφ + √2D ξ_x
        dφ = √2Dr ξ_φ
        """
        self.get_force()
        # active
        self.r += self.v * self.dt
        # Brownian
        self.r += self.D / self.kT * self.f * self.dt + \
                  np.sqrt(2 * self.D * self.dt) * np.random.randn(self.N, self.dim)
        # Re-orient
        self.rot_diffuse()
        self.v = self.o * self.v0
        self.fix_boundary()
        self.notify()

    def move(self):
        self.move_overdamp()


@Boundary("pbc")
class ABP2DWCAPBC(ABP2D):
    def __init__(self, N, dt, Pe, **kwargs):
        ABP2D.__init__(self, N, dt, Pe, **kwargs)
        self.force_func = force_wca


class Lavergne_2019_Science(ABP2D):
    def __init__(self, N, dt, Pe, alpha, p_act, R0=24, **kwargs):
        """
        Args:
            N (int): the number of particles
            Pe (float): The Péclet number
            alpha (float): the vision cone range, from 0 to π
            p_act (float): the perception activation threshold,
                scaled by the centre perception. p_act = P* / Pc
                (see Fig.4 where p_act range from 0 to 1.2)
        """
        ABP2D.__init__(self, N, dt, Pe, **kwargs)
        self.alpha = alpha
        density = self.N / (np.pi * R0**2)
        self.p_centre = self.alpha / np.pi * R0 * density
        self.p_act = p_act * self.p_centre
        xy = np.random.randn(self.dim, self.N)
        r = np.linalg.norm(xy, axis=0)
        l = np.random.uniform(0, R0**2, self.N) ** (1 / self.dim)
        self.r = (xy / r * l).T
        self.is_active = np.ones(self.N, dtype=bool)

    def get_perception(self):
        rij = self.get_pair_shift()  # 2, N, N
        """
        `rij` is something like
        r11 | r12 | r13
        r21 | r22 | r23
        r31 | r32 | r33
        where θij is the angle of vector from i to j
        """
        dist = np.linalg.norm(rij, axis=0)  # N, N
        np.fill_diagonal(dist, np.inf)
        angle = np.arctan2(rij[1], rij[0]) % (2 * np.pi)  # N, N
        """
        `angle` is something like
        θ11 | θ21 | θ31
        θ12 | θ22 | θ32
        θ13 | θ23 | θ33
        where θij is the angle of vector from i to j
        """
        angle_diff = angle - self.phi[np.newaxis, :]
        angle_diff = angle_diff - np.rint(angle_diff / 2 / np.pi) * 2 * np.pi
        """
        `angle_diff` is something like
        θ11 - φ1 | θ21 - φ2 | θ31 - φ3
        θ12 - φ1 | θ22 - φ2 | θ32 - φ3
        θ13 - φ1 | θ23 - φ2 | θ33 - φ3
        sum over rows for angle_diff per particle
        """
        adj_mat = np.abs(angle_diff) < self.alpha
        return np.sum(adj_mat / dist, axis=0) / (2 * np.pi)

    def move_overdamp(self):
        """
        perception = get_perception
        dx = v0(perception) cosφ + √2D ξ_x
        dx = v0(perception) sinφ + √2D ξ_x
        dφ = √2Dr ξ_φ
        """
        self.is_active = self.get_perception() > self.p_act
        # active
        self.r += self.v * self.dt * self.is_active[:, np.newaxis]
        # Brownian
        self.get_force()
        self.r += self.D / self.kT * self.f * self.dt  # force
        self.r += np.sqrt(2 * self.D * self.dt) * np.random.randn(self.N, self.dim)
        # Re-orient
        self.phi += np.sqrt(2 * self.Dr * self.dt) * np.random.randn(self.N)
        self.phi = self.phi % (2 * np.pi)
        self.o = np.array((np.cos(self.phi), np.sin(self.phi))).T  # (n, 2)
        self.v = self.o * self.v0
        self.fix_boundary()
        self.notify()


class Vision(Lavergne_2019_Science):
    def move_overdamp(self):
        """
        perception = get_perception
        dx = v0(perception) cosφ + √2D ξ_x
        dx = v0(perception) sinφ + √2D ξ_x
        dφ = √2Dr ξ_φ
        """
        perception, target = self.get_perception()
        self.is_active =  perception > self.p_act
        # active
        self.r += self.v * self.dt
        # Brownian (No translational diffusion)
        self.get_force()
        self.r += self.D / self.kT * self.f * self.dt  # force
        # Re-orient (Orientational diffusion)
        self.phi -= target * self.is_active
        self.phi += np.sqrt(2 * self.Dr * self.dt) * np.random.randn(self.N)
        self.phi = self.phi % (2 * np.pi)
        self.o = np.array((np.cos(self.phi), np.sin(self.phi))).T  # (n, 2)
        self.v = self.o * self.v0
        self.fix_boundary()
        self.notify()

    def get_perception(self):
        rij = self.get_pair_shift()  # 2, N, N
        dist = np.linalg.norm(rij, axis=0)  # N, N
        np.fill_diagonal(dist, np.inf)
        angle = np.arctan2(rij[1], rij[0]) % (2 * np.pi)  # N, N
        angle_diff = self.phi[np.newaxis, :] - angle  # N, N, sum shold be along axis 0
        angle_diff = angle_diff - np.rint(angle_diff / 2 / np.pi) * 2 * np.pi
        adj_mat = np.abs(angle_diff) < self.alpha
        neighbour_num = adj_mat.sum(axis=0)
        target = (angle_diff * adj_mat).sum(axis=0) / neighbour_num
        target[neighbour_num==0]=0
        close = np.argmin(np.abs(angle_diff), axis=0)[None, :]
        target = np.take_along_axis(angle_diff, close, axis=0).ravel()
        return np.sum(adj_mat / dist, axis=0) / (2 * np.pi), target


@Boundary("align_sphere")
class Vision_AS(Vision): pass


class L2SWCA(Lavergne_2019_Science): force_func = lambda self, rij: force_wca(rij)


@Boundary("align_sphere")
class ABP2D_AS(ABP2D): pass


@Boundary("align_sphere")
class ABP2D_AS_WCA(ABP2D): force_func=lambda self, x: force_wca(x)
