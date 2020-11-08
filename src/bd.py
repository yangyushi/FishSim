from abc import ABCMeta, abstractmethod
import matplotlib.pyplot as plt
from itertools import product
from scipy.spatial.distance import pdist, cdist, squareform
import matplotlib.animation as animation
import numpy as np
from numba import njit
from force import force_lj, force_wca


def animate(system, r=100, jump=100, box=None):
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
        fig, update, 100, fargs=(system, scatter), interval=30
    )
    plt.show()


def animate_active_2d(system, r=100, jump=100, box=None):
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
        fig, update, 100, fargs=(system, scatter), interval=30
    )
    plt.show()


class Observer(metaclass=ABCMeta):
    @abstractmethod
    def update(self, subject): pass


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
        self.block = block
        self.count = 0
        self.Ek = 0
        names = ('E/N cut&shifted', 'P cut&shifted', 'T Kinetic', 'T Config')
        print('|'.join([f'{n.split(" ")[0]: ^20}' for n in names]))
        print('|'.join([f'{n.split(" ")[1]: ^20}' for n in names]))
        print('-' * (21 * len(names)))

    def update(self, sys):
        if self.count == self.block - 1:
            e_kin = np.mean(sys.interest['kinetic']) / sys.N
            e_tot = e_kin + np.mean(sys.interest['potential']) / sys.N
            press = sys.kT * sys.density + np.mean(sys.interest['virial']) / sys.volume
            t_kin = 2 * e_kin / sys.dim
            t_con = np.mean(
                np.array(sys.interest['force_sq']) / np.array(sys.interest['laplacian'])
            )
            quantities = (e_tot, press, t_kin, t_con)
            print('|'.join([f'{val: ^20.4f}' for val in quantities]))
            self.count = 0
            sys.interest = {}
        else:
            self.count += 1

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
        The content in the result is rij = ri - rj

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

            def fix_boundary(self):
                self.r %= self.box

            def get_pair_shift(self):
                pos_in_box = self.r.T / self.box
                pair_shift = pos_in_box[:, :, np.newaxis] - pos_in_box[:, np.newaxis, :]
                pair_shift = (pair_shift - np.rint(pair_shift)) * self.box
                return pair_shift

        return SystemPBC

    if condition == "pbc":
        return pbc_decorator
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
        self.get_force()


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
        sigma = np.sqrt(self.kT / (self.m * self.N))
        self.v = np.random.normal(0, sigma, (self.N, self.dim))
        self.phi = np.random.uniform(0, 2 * np.pi, self.N)

    def move_overdamp(self):
        """
        dx = v0 cosφ + √2D ξ_x
        dx = v0 sinφ + √2D ξ_x
        dφ = √2Dr ξ_φ
        """
        self.get_force()
        # active
        self.r[:, 0] += self.v0 * np.cos(self.phi) * self.dt
        self.r[:, 1] += self.v0 * np.sin(self.phi) * self.dt
        # Brownian
        self.r += self.D / self.kT * self.f * self.dt + \
                  np.sqrt(2 * self.D * self.dt) * np.random.randn(self.N, self.dim)
        # Re-orient
        self.phi += np.sqrt(2 * self.Dr * self.dt) * np.random.randn(self.N)
        self.phi = self.phi % (2 * np.pi)
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
        self.r[:, 0] += self.v0 * np.cos(self.phi) * self.dt * self.is_active
        self.r[:, 1] += self.v0 * np.sin(self.phi) * self.dt * self.is_active
        # Brownian
        self.get_force()
        self.r += self.D / self.kT * self.f * self.dt  # force
        self.r += np.sqrt(2 * self.D * self.dt) * np.random.randn(self.N, self.dim)
        # Re-orient
        self.phi += np.sqrt(2 * self.Dr * self.dt) * np.random.randn(self.N)
        self.phi = self.phi % (2 * np.pi)
        self.fix_boundary()
        self.notify()


class L2SWCA(Lavergne_2019_Science):
    force_func = lambda self, rij: force_wca(rij)
