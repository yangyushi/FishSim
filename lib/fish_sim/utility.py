import pickle
import numpy as np
from itertools import product
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.special import gamma as gamma_func
from scipy.spatial.distance import pdist, cdist, squareform

from .noise_3d import add_vicsek_noise_3d
from .force import force_lj, force_wca


class FuncAnimationDisposable(animation.FuncAnimation):
    def __init__(self, fig, func, **kwargs):
        super().__init__(fig, func, **kwargs)

    def _step(self, *args):
        still_going = animation.Animation._step(self, *args)
        if not still_going and self.repeat:
            super()._init_draw()
            self.frame_seq = self.new_frame_seq()
            self.event_source.interval = self._repeat_delay
            return True
        elif (not still_going) and (not self.repeat):
            plt.close()
            return False
        else:
            self.event_source.interval = self._interval
            return still_going

    def _stop(self, *args):
        # On stop we disconnect all of our events.
        if self._blit:
            self._fig.canvas.mpl_disconnect(self._resize_id)
        self._fig.canvas.mpl_disconnect(self._close_id)


def animate(
    system, r=100, jump=100, box=None, save='',
        fps=60, show=False, frames=100, repeat=False, title=""
):
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
    ax.set_title(title)

    end_frame_num = frames - 1

    def update(num, system, scatter):
        for _ in range(jump):
            system.move()
        if system.dim == 3:
            scatter.set_data(system.r.T[:2])
            scatter.set_3d_properties(system.r.T[2])
        else:
            scatter.set_data(system.r.T)
            if (num == end_frame_num) and not repeat:
                raise StopIteration
        return scatter

    scatter = ax.plot(
        *system.r.T, color='teal', mfc='w', ls='None', marker='o',
        markersize=r
    )[0]
    ani = FuncAnimationDisposable(
        fig, update, frames=frames, fargs=(system, scatter), interval=1,
        blit=False, repeat=repeat
    )
    if show:
        plt.show()
    if save:
        ani.save(save, writer='imagemagick', fps=fps)


def animate_active_2d(
        system, r=100, jump=100, box=None,
        save='', fps=60, show=True, frames=100,
        repeat=False, title=''
):
    fig = plt.figure(figsize=(5, 5), tight_layout=True)
    ax = fig.add_subplot()
    ax.set_title(title)
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
        if (num == frames - 1) and (not repeat):
            raise StopIteration
        return scatter

    scatter = ax.plot(
        *system.r.T, color='teal', mfc='w', ls='None', marker='o',
        markersize=r
    )[0]
    quiver = ax.quiver(
        *system.r.T, np.cos(system.phi), np.sin(system.phi),
        pivot='mid', units='width', color='teal', zorder=5, scale=250/r,
    )
    ani = FuncAnimationDisposable(
        fig, update, frames=frames, fargs=(system, scatter), interval=1,
        blit=False, repeat=repeat
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

    Attributes:
        report (bool): if true, the information will be printed on-the-fly
        result (dict): a collection of all calculation resuls. The elements\
            are, `E/N`, `P`, `T_kinetic`, and\
            `T_configuration`.
    """
    def __init__(self, block, report=True):
        Observer.__init__(self, block)
        self.Ek = 0
        names = ('E/N cut&shifted', 'P cut&shifted', 'T Kinetic', 'T Config')
        self.report = report
        self.result = {
            'E/N': [],
            'P': [],
            'T_kinetic': [],
            'T_configuration': [],
        }
        if self.report:
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
        self.result['E/N'].append(e_tot)
        self.result['P'].append(press)
        self.result['T_kinetic'].append(t_kin)
        self.result['T_configuration'].append(t_con)
        if self.report:
            quantities = (e_tot, press, t_kin, t_con)
            print('|'.join([f'{val: ^20.4f}' for val in quantities]))


class Dynamic(Observer):
    def __init__(self, block, report=True):
        Observer.__init__(self, block)
        self.names = ['polarisation']
        self.result_frames = {
            'polarisation': np.empty(block),
        }
        self.result = {
            'polarisation':[],
        }
        self.report = report
        if self.report:
            print('|'.join([f'{name: ^20}' for name in self.names]))
            print('-' * (21 * len(self.names)))

    def aggregate(self, system):
        for key in self.result:
            self.result[key].append(self.result_frames[key].mean())
        if self.report:
            print('|'.join([
                f'{self.result[name][-1]: ^20.4f}' for name in self.names
            ]))

    def collect(self, system):
        orient = system.v / np.linalg.norm(system.v, axis=1)[:, np.newaxis]
        pol = np.linalg.norm(np.mean(orient, axis=0))
        self.result_frames['polarisation'] = pol


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
