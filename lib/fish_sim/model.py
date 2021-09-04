from itertools import product

import numpy as np
from scipy.special import gamma as gamma_func

from .force import force_lj, force_wca
from .noise_3d import add_vicsek_noise_3d


class BD():
    """
    The basic Brownian Dynamic simulation of ideal gas
    """
    def __init__(self, n, dim, dt, gamma=None, kT=None, m=None, D=None, **kwargs):
        self.n, self.dim, self.dt = n, dim, dt
        self.positions = np.zeros((self.dim, self.n))  # positions
        self.gamma = gamma
        self.kT = kT
        self.m = m
        self.D = D
        self.update_parameters()
        self.f = np.zeros_like(self.positions)
        self.interest = {}  # collect quantity of interest
        self.__observers = []  # ovserver design pattern
        sigma = np.sqrt(self.kT / (self.m * self.n))
        self.velocities = np.random.normal(0, sigma, (self.dim, self.n))

    def get_positions(self):
        return self.positions

    def get_velocities(self):
        return self.velocities

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
        return self.positions[:, :, np.newaxis] - self.positions[:, np.newaxis, :]

    def fix_boundary(self): pass

    def force_func(self, shift): return np.zeros((self.dim, self.n)), {}

    def get_force(self):
        rij = self.get_pair_shift()
        self.f, interest = self.force_func(rij)
        if self.m:
            interest['kinetic'] = 0.5 * self.m * np.sum(self.velocities ** 2)
        else:
            interest['kinetic'] = 0.5 * np.sum(self.velocities ** 2)
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
        self.velocities += 0.5 * self.dt * self.f / self.m  # b propagator

        self.positions += 0.5 * self.dt * self.velocities  # a propagator
        self.fix_boundary()

        # --- o propagator ---
        x = self.gamma * self.dt
        c = np.sqrt(
            1 - np.exp(- 2 * x) if x > 0.0001 else \
                np.polyval([-2/3,4/3,-2.0,2.0,0.0], x)  # Tayler Expansion
        )
        self.velocities = self.velocities * np.exp(-x) +\
            c * np.sqrt(self.m * self.kT) * np.random.randn(self.dim, self.n)
        # -------

        self.positions += 0.5 * self.dt * self.velocities  # a propagator
        self.fix_boundary()

        self.get_force()

        self.velocities += 0.5 * self.dt * self.f / self.m  # b propagator

        self.notify()

    def move_overdamp(self):
        self.get_force()
        self.positions += self.D / self.kT * self.f * self.dt + \
                  np.sqrt(2 * self.D * self.dt) * np.random.randn(self.dim, self.n)
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
                self.density = self.n / self.volume
                repeat = np.ceil(np.power(self.n, 1 / self.dim)).astype(int)
                lattice = np.array(list(product(*[range(repeat)] * self.dim)))
                self.positions = lattice[:self.n].T / repeat * kwargs['box']
                self.get_force()

            def fix_boundary(self):
                self.positions %= self.box

            def get_pair_shift(self):
                pos_in_box = self.positions / self.box  # (dim, n)
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
                self.density = self.n / self.volume
                p = np.random.randn(self.dim, self.n)  # positions
                r = np.linalg.norm(p, axis=0)  # radii
                l = np.random.uniform(0, self.R**self.dim, self.n) ** (1 / self.dim)
                self.positions = p / r * l  # uniform inside n-sphere
                self.get_force()

            def fix_boundary(self):
                """
                fish that are outside of the circle should not move further
                """
                v0 = self.velocities.copy()
                is_outside = np.linalg.norm(self.positions, axis=0) >= self.R
                if not np.any(is_outside):
                    return
                v_orient = self.velocities[:, is_outside]  # dim, n
                r_orient = self.positions[:, is_outside]  # dim, n
                v_mag = np.linalg.norm(v_orient, axis=0)  # n,
                r_mag = np.linalg.norm(r_orient, axis=0)  # n,
                v_orient /= v_mag[np.newaxis, :]  # dim, n
                r_orient /= r_mag[np.newaxis, :]  # dim, n
                angles = np.arccos(np.einsum('ij,ij->j', v_orient, r_orient))  # n,
                signs = np.sign(np.cross(v_orient.T, r_orient.T).T)  # dim, n
                # if greater than pi/2, no change, otherwise change to pi/2
                should_fix = np.abs(angles) < np.pi/2  # n
                if not np.any(should_fix):
                    return
                if self.dim == 2:
                    r_angle = np.arctan2(r_orient[1], r_orient[0])  # n,
                    v_orient[0, should_fix] = np.cos(
                        r_angle[:, should_fix] - signs[:, should_fix] * np.pi/2
                    ) * v_mag[should_fix]
                    v_orient[1, should_fix] = np.sin(
                        r_angle[:, should_fix] - signs[:, should_fix] * np.pi/2
                    ) * v_mag[should_fix]
                    self.velocities[:, is_outside] = v_orient
                    self.phi[is_outside] = np.arctan2(v_orient[1], v_orient[0])
                elif self.dim == 3:
                    e_align = np.cross(
                        r_orient[:, should_fix].T,
                        np.cross(
                            v_orient[:, should_fix].T,
                            r_orient[:, should_fix].T
                        ),
                    ).T
                    e_align = e_align / np.linalg.norm(e_align, axis=0)[None, :]
                    to_fix = is_outside.copy()
                    to_fix[to_fix] *= should_fix   # to fix = outside & angle < np.pi/2
                    self.o[:, to_fix] = e_align
                    self.velocities[:, is_outside] =\
                        self.o[:, is_outside] * v_mag[np.newaxis, :]
                else:
                    raise ValueError("Only 2D and 3D systems are supported")

        return SystemAS

    def align_half_sphere(system):
        class SystemAHS(system):
            """
            Particles aligh with the boundary formed by the cut of
                x^2 + y^2 + ... = R  (sphere) and
                y = 0 or z = 0  (cap)
            """
            def __init__(self, *args, **kwargs):
                self.R = kwargs['R']  # radius of the bounding sphere
                system.__init__(self, *args, **kwargs)
                try:
                    self.R = float(self.R)
                except (ValueError, TypeError) as err:
                    raise ValueError("the radius should be a numerical value")
                self.volume = np.pi**(self.dim / 2) /\
                              gamma_func(self.dim/2 + 1) * self.R**self.dim
                self.density = self.n / self.volume
                p = np.random.randn(self.dim, self.n)  # positions
                r = np.linalg.norm(p, axis=0)  # radii
                l = np.random.uniform(0, self.R**self.dim, self.n) ** (1 / self.dim)
                self.positions = p / r * l  # uniform inside n-sphere
                self.positions[self.positions[:, self.dim-1] > 0, self.dim-1] *= -1
                self.get_force()

            def fix_boundary(self):
                """
                fish that are outside of the circle should not move further
                """
                v0 = self.velocities.copy()
                out_sphere = np.linalg.norm(self.positions, axis=0) >= self.R
                out_cap = self.positions[:, self.dim - 1] > 0
                is_outside = np.logical_or(out_sphere, out_cap)
                if not np.any(is_outside):
                    return
                excess_cap = self.positions[is_outside, self.dim - 1]  # shape (N', ) prime means N outside
                excess_sph = np.linalg.norm(self.positions[is_outside], axis=1) - self.R  # shape (N', )
                v_orient = self.velocities[is_outside]  # n, 2
                r_orient = self.positions[is_outside]  # n, 2
                v_mag = np.linalg.norm(v_orient, axis=1)  # n,
                r_mag = np.linalg.norm(r_orient, axis=1)  # n,
                v_orient /= v_mag[:, np.newaxis]
                r_orient /= r_mag[:, np.newaxis]
                angles = np.arccos(np.einsum('ij,ij->i', v_orient, r_orient))  # n,
                # if greater than pi/2, no change, otherwise change to pi/2
                should_fix_cap = v_orient[:, self.dim - 1] > 0
                should_fix_sph = np.abs(angles) < np.pi/2  # n
                should_fix_cap[excess_cap < excess_sph] = False
                should_fix_sph[excess_cap >= excess_sph] = False
                should_fix = np.logical_or(should_fix_sph, should_fix_cap)
                if not np.any(should_fix):
                    return
                if self.dim == 2:
                    # fix the particles that exceeds the sphere
                    if np.sum(should_fix_sph) > 0:
                        signs = np.sign(np.cross(v_orient, r_orient))
                        r_angle = np.arctan2(r_orient[:, 1], r_orient[:, 0])  # n,
                        v_orient[should_fix_sph, 0] = np.cos(
                            r_angle[should_fix_sph] - signs[should_fix_sph] * np.pi/2
                        ) * v_mag[should_fix_sph]
                        v_orient[should_fix_sph, 1] = np.sin(
                            r_angle[should_fix_sph] - signs[should_fix_sph] * np.pi/2
                        ) * v_mag[should_fix_sph]
                    # fix the particles that exceeds the cap
                    fix_cap_num = np.sum(should_fix_cap)
                    if fix_cap_num > 0:
                        signs = np.sign(v_orient[should_fix_cap, 0])
                        cap_fixed = np.repeat(np.array((1, 0))[np.newaxis, :], fix_cap_num, axis=0)
                        v_orient[should_fix_cap] = cap_fixed
                        v_orient[should_fix_cap, 0] *= signs
                        v_orient[should_fix_cap] *= v_mag[should_fix_cap][:, np.newaxis]
                    # update velocity
                    self.velocities[is_outside, :] = v_orient
                    self.phi[is_outside] = np.arctan2(v_orient[:, 1], v_orient[:, 0])
                elif self.dim == 3:
                    # fix the particles that exceeds the sphere
                    if np.sum(should_fix_sph) > 0:
                        e_align = np.cross(
                            r_orient[should_fix_sph],
                            np.cross(v_orient[should_fix_sph], r_orient[should_fix_sph]),
                        )
                        e_align = e_align / np.linalg.norm(e_align, axis=1)[:, np.newaxis]
                        to_fix = is_outside.copy()
                        to_fix[to_fix] *= should_fix_sph   # to fix = outside & angle < np.pi/2
                        self.o[to_fix] = e_align
                    # fix the particles that exceeds the cap
                    if np.sum(should_fix_cap) > 0:
                        e_align = v_orient[should_fix_cap].copy()
                        e_align[:, -1] = 0
                        e_align /= np.linalg.norm(e_align, axis=1)[:, np.newaxis]
                        to_fix = is_outside.copy()
                        to_fix[to_fix] *= should_fix_cap   # to fix = outside & angle < np.pi/2
                        self.o[to_fix] = e_align
                    # update the velocity
                    self.velocities[is_outside] = self.o[is_outside] * v_mag[:, np.newaxis]
                else:
                    raise ValueError("Only 2D and 3D systems are supported")
        return SystemAHS

    def align_fish_bowl(system):
        class SystemAFB(system):
            """
            Particles aligh with the boundary formed by the cut of
                y = c * x^2 or z = c * (x^2 + y^2)  (hyperbolic) and
                y = 0 or z = 0  (cap)

            Useful:
                z = c * r^2
                r = sqrt(z / c)
            """
            def __init__(self, *args, **kwargs):
                self.c = kwargs['c']
                self.z_max = kwargs['z_max']
                self.r_max = np.sqrt(self.z_max / self.c)
                system.__init__(self, *args, **kwargs)
                try:
                    self.z_max = float(self.z_max)
                except (ValueError, TypeError) as err:
                    raise ValueError("z_max should be a numerical value")
                try:
                    self.c = float(self.c)
                except (ValueError, TypeError) as err:
                    raise ValueError("c should be a numerical value")
                rand_z = np.random.uniform(0, self.z_max, self.n)
                rand_r = np.random.uniform(-1, 1, self.n) * np.sqrt(rand_z / self.c)
                if self.dim == 2:
                    self.volume = 2 * (
                       self.z_max * self.r_max -\
                       2 * self.c / 3 * self.r_max**3
                    )
                    self.positions = np.array((rand_z, rand_r))
                elif self.dim == 3:
                    self.volume = np.pi / self.c * self.z_max**2
                    rand_phi = np.random.uniform(0, np.pi * 2, self.n)
                    self.positions = np.array((
                        rand_r * np.cos(rand_phi),
                        rand_r * np.sin(rand_phi),
                        rand_z
                    ))
                else:
                    raise ValueError("Invalid Dimension")
                self.density = self.n / self.volume
                self.get_force()

            def __project(self, r, z):
                r"""
                Calculate the projection of 2D point in cylindar coordinate system to a hyprobolic function

                .. math::

                    y = \textsf{self.c} \cdot x^2

                Args:
                    r (:obj:`float` or :obj:`numpy.ndarray`): the radii of a point in cylinder coordinate system
                    z (:obj:`float` or :obj:`numpy.ndarray`): the z (height) of a point

                Return:
                    (:obj:`float` or :obj:`numpy.ndarray`): The projected coordinates in 2D
                """
                p = 2**(1/3)
                term = 108 * self.c ** 4 * r + np.sqrt(
                    11664 * self.c**8 * r**2 +
                    864 * self.c**6 * (1 - 2 * self.c * z)**3
                )
                term = np.power(term, 1/3)
                radii_proj = -(p * (1 - 2 * self.c * z)) / term +\
                    term / (6 * p * self.c**2)
                z_proj = self.c * radii_proj ** 2
                return radii_proj, z_proj

            def fix_boundary(self):
                """
                fish that are outside of the boundary should not move further
                """
                radii = np.linalg.norm(  # shape (n, )
                    self.positions[:self.dim - 1, :], axis=0
                )
                z = self.positions[self.dim - 1, :]  # it is actually y in 2D, shape n
                out_hyp = self.c * radii**2 > z  # n,
                out_cap = z > self.z_max  # n,
                is_outside = np.logical_or(out_hyp, out_cap)  # n,
                if not np.any(is_outside):
                    return
                radii_out = radii[is_outside]  # n,
                z_out = z[is_outside]  # n,
                radii_out_proj, z_out_proj = self.__project(radii_out, z_out)  # n,
                excess_cap = z_out - self.z_max
                excess_hyp = np.sqrt(
                    (radii_out - radii_out_proj)**2 + (z_out - z_out_proj)**2
                )
                excess_hyp[self.c * radii_out**2 < z_out] = 0  # ignore the inside cases

                v_orient = self.velocities[:, is_outside]  # dim, n
                v_mag = np.linalg.norm(v_orient, axis=0)  # n,
                v_orient /= v_mag[np.newaxis, :]  # dim, n
                if self.dim == 2:
                    r_proj = np.array((  # dim, n
                        np.sign(self.positions[0, is_outside]) * radii_out_proj,  # projected x
                        z_out_proj  # projected y
                    ))
                elif self.dim == 3:
                    phi = np.arctan2(
                        self.positions[1, is_outside], self.positions[0, is_outside]
                    )
                    r_proj = np.array((  # dim, n
                        radii_out_proj * np.cos(phi),
                        radii_out_proj * np.sin(phi),
                        z_out_proj
                    ))
                else:
                    raise ValueError("Only 2D and 3D systems are supported")
                r_orient = self.positions[:, is_outside] - r_proj # dim, n
                r_mag = np.linalg.norm(r_orient, axis=0)  # n,
                r_orient /= r_mag[np.newaxis, :]
                angles = np.arccos(np.einsum('ij,ij->j', v_orient, r_orient))  # n,
                # if greater than pi/2, no change, otherwise change to pi/2
                should_fix_cap = v_orient[self.dim-1, :] > 0  # n,
                should_fix_hyp = np.abs(angles) < np.pi / 2  # n
                should_fix_cap[excess_cap < excess_hyp] = False
                should_fix_hyp[excess_cap >= excess_hyp] = False
                should_fix = np.logical_or(should_fix_hyp, should_fix_cap)
                if not np.any(should_fix):
                    return
                if self.dim == 2:
                    # fix the particles that exceeds the hyperbolic
                    if np.sum(should_fix_hyp) > 0:
                        signs = np.sign(np.cross(v_orient.T, r_orient.T).T)
                        r_angle = np.arctan2(r_orient[1], r_orient[0])  # n,
                        v_orient[0, should_fix_hyp] = np.cos(
                            r_angle[should_fix_hyp] - signs[should_fix_hyp] * np.pi/2
                        ) * v_mag[should_fix_hyp]
                        v_orient[1, should_fix_hyp] = np.sin(
                            r_angle[should_fix_hyp] - signs[should_fix_hyp] * np.pi/2
                        ) * v_mag[should_fix_hyp]
                    # fix the particles that exceeds the cap
                    fix_cap_num = np.sum(should_fix_cap)
                    if fix_cap_num > 0:
                        signs = np.sign(v_orient[0, should_fix_cap])
                        cap_fixed = np.repeat(
                            np.array((1, 0))[:, np.newaxis], fix_cap_num, axis=1
                        )  # dim, n
                        v_orient[should_fix_cap] = cap_fixed
                        v_orient[0, should_fix_cap] *= signs
                        v_orient[should_fix_cap] *= v_mag[should_fix_cap][np.newaxis, :]
                    # update velocity
                    self.velocities[:, is_outside] = v_orient
                    self.phi[is_outside] = np.arctan2(v_orient[1], v_orient[0])
                elif self.dim == 3:
                    # fix the particles that exceeds the hyperbolic
                    if np.sum(should_fix_hyp) > 0:
                        e_align = np.cross(
                            r_orient[:, should_fix_hyp].T,
                            np.cross(
                                v_orient[:, should_fix_hyp].T,
                                r_orient[:, should_fix_hyp].T
                            ),
                        ).T  # dim, n
                        e_align = e_align / np.linalg.norm(e_align, axis=0)[np.newaxis, :]
                        to_fix = is_outside.copy()
                        to_fix[to_fix] *= should_fix_hyp   # to fix = outside & angle < np.pi/2
                        self.o[:, to_fix] = e_align
                    # fix the particles that exceeds the cap
                    if np.sum(should_fix_cap) > 0:
                        e_align = v_orient[:, should_fix_cap].copy()
                        e_align[-1] = 0
                        e_align /= np.linalg.norm(e_align, axis=0)[np.newaxis, :]
                        to_fix = is_outside.copy()
                        to_fix[to_fix] *= should_fix_cap   # to fix = outside & angle < np.pi/2
                        self.o[:, to_fix] = e_align
                    # update the velocity
                    self.velocities[:, is_outside] = self.o[:, is_outside] *\
                        v_mag[np.newaxis, :]
                else:
                    raise ValueError("Only 2D and 3D systems are supported")
        return SystemAFB

    if condition == "pbc":
        return pbc_decorator
    if condition == "align_sphere":
        return align_sphere
    if condition == "align_half_sphere":
        return align_half_sphere
    if condition == "align_fish_bowl":
        return align_fish_bowl
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
    def __init__(self, n, r, eta, v0, **kwargs):
        kwargs.update({'m': 1, 'D': 1, 'kT': 1})  # meaningless in Vicsek Model
        BD.__init__(self, n, dim=3, dt=1, **kwargs)
        self.v0 = v0  # speed
        self.r0 = r   # interaction range
        self.eta = eta
        self.positions = np.random.randn(self.dim, self.n)
        self.o = np.random.randn(self.dim, self.n)  # orientation
        self.o = self.o / np.linalg.norm(self.o, axis=0)[np.newaxis, :]
        self.velocities = self.o * self.v0

    def move(self):
        self.get_force()
        # active
        self.positions += self.velocities * self.dt

        rij = self.get_pair_shift()
        dij = np.linalg.norm(rij, axis=0)  # distances, shape (N, N)
        aij = dij < self.r0  # adjacency matrix (N, N)
        nni = np.sum(aij, axis=0)  # number of neighbours for each particle

        vii = self.velocities[:, np.newaxis, :] * np.ones((3, self.n, 1))  # (3, N, N)
        v_mean = np.sum(vii * aij[np.newaxis, :, :], axis=2) / nni  # (3, N)
        v_xy = np.linalg.norm(v_mean[:2], axis=0)  # (N,)
        azi = np.arctan2(v_mean[1], v_mean[0]) # (N,)
        ele = np.arctan(v_mean[2] / v_xy)  # (N,)
        azi, ele = add_vicsek_noise_3d(azi, ele, eta=self.eta)
        e_xy = np.cos(ele)
        self.o[2] = np.sin(ele)
        self.o[1] = e_xy * np.sin(azi)
        self.o[0] = e_xy * np.cos(azi)
        self.velocities = self.o * self.v0

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
    def __init__(self, n, dt, Pe, **kwargs):
        BD.__init__(self, n, 2, dt, **kwargs)
        if not self.D:
            raise KeyError("Missing the Diffusion constant of the system")
        self.Pe = Pe
        self.v0 = self.Pe * self.D
        self.Dr = 3 * self.D
        self.positions = np.random.randn(self.dim, self.n)
        self.o = np.random.randn(self.dim, self.n)  # orientation
        self.o = self.o / np.linalg.norm(self.o, axis=0)[np.newaxis, :]
        self.velocities = self.o * self.v0
        self.phi = np.random.uniform(0, 2 * np.pi, self.n)

    def move_overdamp(self):
        """
        dx = v0 cosφ + √2D ξ_x
        dy = v0 sinφ + √2D ξ_x
        dφ = √2Dr ξ_φ
        """
        self.get_force()
        # active
        self.positions += self.velocities * self.dt
        # Brownian
        self.positions += self.D / self.kT * self.f * self.dt + \
                  np.sqrt(2 * self.D * self.dt) * np.random.randn(self.dim, self.n)
        # Re-orient
        self.phi += np.sqrt(2 * self.Dr * self.dt) * np.random.randn(self.n)
        self.phi = self.phi % (2 * np.pi)
        self.o = np.array((np.cos(self.phi), np.sin(self.phi)))  # (2, n)
        self.velocities = self.o * self.v0
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
    def __init__(self, n, dt, Pe, **kwargs):
        BD.__init__(self, n, dim=3, dt=dt, **kwargs)
        if not self.D:
            raise KeyError("Missing the Diffusion constant of the system")
        self.Pe = Pe
        self.v0 = self.Pe * self.D
        self.Dr = 3 * self.D
        self.positions = np.random.randn(self.dim, self.n)
        self.o = np.random.randn(self.dim, self.n)  # orientation
        self.o /= np.linalg.norm(self.o, axis=0)[np.newaxis, :]
        self.velocities = self.o * self.v0

    def rot_diffuse(self):
        """
        Perform the rotational diffusion
        (https://github.com/FTurci/active-lammps, I copied Francesco's code)
        (Winkler et al. Soft Matter, 2015, 11, 6680)
        """
        ox, oy, oz = self.o
        cos_theta = oz
        sin_theta = np.sqrt(1 - oz**2)
        phi = np.arctan2(oy, ox)
        cos_phi = np.cos(phi)
        sin_phi = np.sin(phi)

        e_theta = np.empty((3, self.n))
        e_phi = np.empty((3, self.n))
        e_theta[0] =  cos_phi * cos_theta
        e_theta[1] =  sin_phi * cos_theta
        e_theta[2] =           -sin_theta
        e_phi[0]   = -sin_phi
        e_phi[1]   =  cos_phi
        e_phi[2]   =  0.0

        r1, r2 = np.random.randn(2, self.n) * np.sqrt(2 * self.Dr)
        dt_wiener = np.sqrt(self.dt)
        for i in range(3):
            self.o[i, :] += e_theta[i] * dt_wiener * r1 +\
                            e_phi[i]   * dt_wiener * r2 -\
                            2 * self.Dr * self.dt * self.o[i, :]
        self.o /= np.linalg.norm(self.o, axis=0)[np.newaxis, :]

    def move_overdamp(self):
        """
        dx = v0 cos cosφ + √2D ξ_x
        dy = v0 sinφ + √2D ξ_x
        dz = v0 sinφ + √2D ξ_x
        dφ = √2Dr ξ_φ
        """
        self.get_force()
        # active
        self.positions += self.velocities * self.dt
        # Brownian
        self.positions += self.D / self.kT * self.f * self.dt + \
                  np.sqrt(2 * self.D * self.dt) * np.random.randn(self.dim, self.n)
        # Re-orient
        self.rot_diffuse()
        self.velocities = self.o * self.v0
        self.fix_boundary()
        self.notify()

    def move(self):
        self.move_overdamp()


@Boundary("pbc")
class ABP2DWCAPBC(ABP2D):
    def __init__(self, n, dt, Pe, **kwargs):
        ABP2D.__init__(self, n, dt, Pe, **kwargs)
        self.force_func = force_wca


class Lavergne_2019_Science(ABP2D):
    def __init__(self, n, dt, Pe, alpha, p_act, R0=24, **kwargs):
        """
        Args:
            n (int): the number of particles
            Pe (float): The Péclet number
            alpha (float): the vision cone range, from 0 to π
            p_act (float): the perception activation threshold,
                scaled by the centre perception. p_act = P* / Pc
                (see Fig.4 where p_act range from 0 to 1.2)
        """
        ABP2D.__init__(self, n, dt, Pe, **kwargs)
        self.alpha = alpha
        density = self.n / (np.pi * R0**2)
        self.p_centre = self.alpha / np.pi * R0 * density
        self.p_act = p_act * self.p_centre
        xy = np.random.randn(self.dim, self.n)
        r = np.linalg.norm(xy, axis=0)
        l = np.random.uniform(0, R0**2, self.n) ** (1 / self.dim)
        self.positions = xy / r * l
        self.is_active = np.ones(self.n, dtype=bool)

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
        self.positions += self.velocities * self.dt * self.is_active[np.newaxis, :]
        # Brownian
        self.get_force()
        self.positions += self.D / self.kT * self.f * self.dt  # force
        self.positions += np.sqrt(2 * self.D * self.dt) * np.random.randn(self.dim, self.n)
        # Re-orient
        #self.phi -= target * self.is_active
        self.phi += np.sqrt(2 * self.Dr * self.dt) * np.random.randn(self.n)
        self.phi = self.phi % (2 * np.pi)
        self.o = np.array((np.cos(self.phi), np.sin(self.phi)))  # (n, 2)
        self.velocities = self.o * self.v0
        self.fix_boundary()
        self.notify()


class Vision3D(ABP3D):
    def __init__(self, n, dt, Pe, alpha, R0=24, **kwargs):
        """
        Args:
            n (int): the number of particles
            Pe (float): The Péclet number
            alpha (float): the vision cone range, from 0 to π
        """
        ABP3D.__init__(self, n, dt, Pe, **kwargs)
        self.alpha = alpha
        xy = np.random.randn(self.dim, self.n)
        r = np.linalg.norm(xy, axis=0)
        l = np.random.uniform(0, R0**2, self.n) ** (1 / self.dim)
        self.positions = xy / r * l
        self.is_active = np.ones(self.n, dtype=bool)
        self.p_act = 0


    def move_overdamp(self):
        perception, target = self.get_perception()
        self.is_active =  perception > self.p_act
        # active
        self.positions += self.velocities * self.dt
        # Brownian (No translational diffusion)
        self.get_force()
        self.positions += self.D / self.kT * self.f * self.dt  # force
        # Re-orient (Orientational diffusion)
        self.o[self.is_active] = target[self.is_active]
        self.rot_diffuse()
        self.velocities = self.o * self.v0
        self.fix_boundary()
        self.notify()

    def get_perception(self):
        rij = self.get_pair_shift()  # 3, N, N
        dist = np.linalg.norm(rij, axis=0)  # N, N
        np.fill_diagonal(dist, np.inf)
        eij = rij / dist[np.newaxis, :, :]  # 3, N, N
        # self.o boardcast to 3, N, N; should sum over rows
        oi = self.o[:, np.newaxis, :] * np.ones((1, self.n, 1))  # 3, N, N
        angle_diff = np.arccos(np.einsum('ijk,ijk->jk', eij, oi))  # N, N
        angle_diff_ij = angle_diff.copy()
        np.fill_diagonal(angle_diff_ij, np.inf)
        min_diff = np.argmin(angle_diff_ij, axis=0)
        target = eij[:, min_diff, np.arange(self.n)].T  # N, 3
        adj_mat = np.abs(angle_diff) < self.alpha
        np.fill_diagonal(adj_mat, 0)
        return np.sum(adj_mat / dist, axis=0) / (2 * np.pi), target


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
        self.positions += self.velocities * self.dt
        # Brownian (No translational diffusion)
        self.get_force()
        self.positions += self.D / self.kT * self.f * self.dt  # force
        # Re-orient (Orientational diffusion)
        self.phi -= target * self.is_active
        self.phi += np.sqrt(2 * self.Dr * self.dt) * np.random.randn(self.n)
        self.phi = self.phi % (2 * np.pi)
        self.o = np.array((np.cos(self.phi), np.sin(self.phi)))  # (dim, n)
        self.velocities = self.o * self.v0
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
        close = np.argmin(np.abs(angle_diff), axis=0)
        target = angle_diff[close, np.arange(self.n)]
        return np.sum(adj_mat / dist, axis=0) / (2 * np.pi), target


@Boundary("align_sphere")
class Vision_AS(Vision): pass


class L2SWCA(Lavergne_2019_Science): force_func = lambda self, rij: force_wca(rij)


@Boundary("align_sphere")
class ABP2D_AS(ABP2D): pass


@Boundary("align_sphere")
class ABP2D_AS_WCA(ABP2D): force_func=lambda self, x: force_wca(x)
