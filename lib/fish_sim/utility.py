import pickle
import numpy as np
from itertools import product
import matplotlib.pyplot as plt
from force import force_lj, force_wca
import matplotlib.animation as animation
from noise_3d import add_vicsek_noise_3d
from scipy.special import gamma as gamma_func
from scipy.spatial.distance import pdist, cdist, squareform


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
        self.__block = block
        self.__count = 0

    def update(self, system):
        if self.__count == self.__block - 1:
            self.collect(system)
            self.aggregate(system)
            self.__count = 0
            system.interest = {}
        else:
            self.collect(system)
            self.__count += 1

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


class DumpXYZ(Observer):
    """
    Dump the configurations (positions & velocities) for each block
        into a xyz file
    """
    def __init__(self, frequency, filename):
        Observer.__init__(self, frequency)
        if '.xyz' == filename[-4:]:
            fname = filename
        else:
            fname = filename + '.xyz'
        self.f = open(fname, 'w')
        self.count = 0
        self.active = True

    def aggregate(self, system):
        """
        Append many frames to an xyz file
        """
        if self.active:
            configuration = np.concatenate((system.r, system.v), axis=1)
            np.savetxt(
                self.f, configuration, delimiter='\t',
                fmt=['A\t%.8e'] + ['%.8e' for i in range(2 * system.dim - 1)],
                comments='',
                header='%s\nframe %s' % (system.N, self.count)
            )
            self.count += 1

    def stop(self):
        self.active = False

    def start(self):
        self.active = True


class DumpModel(Observer):
    def __init__(self, frequency, filename):
        Observer.__init__(self, frequency)
        if '.pkl' == filename[-4:]:
            fname = filename
        else:
            fname = filename + '.pkl'
        self.f = open(fname, 'wb')
        self.positions, self.velocities = [], []
        self.active = True
        self.not_dumped = True

    def aggregate(self, system):
        if self.active and self.not_dumped:
            self.positions.append(system.r)
            self.velocities.append(system.v)

    def stop(self):
        self.active = False

    def start(self):
        self.active = True

    def dump(self, Model):
        """
        The constructor of class Model should be

        ..code-block::

            Model(positions, velocities)
        """
        p = np.array(self.positions)
        v = np.array(self.velocities)
        model = Model(p, v)
        pickle.dump(model, self.f)
        self.f.close()
        self.not_dumped = False
