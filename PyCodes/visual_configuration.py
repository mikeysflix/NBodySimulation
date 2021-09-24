import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
from units_configuration import *

class AnimationConfiguration(Quantity):

    def __init__(self, savedir=None):
        super().__init__()
        self.savedir = savedir
        self._update_frame_ndim = None

    @property
    def update_frame_ndim(self):
        return self._update_frame_ndim

    @staticmethod
    def update_legend_counter(i, leg, subtitle, time_coordinates, time_units):
        title = '$t = {:.2f}$ {}s'.format(
            time_coordinates[i],
            time_units)
        if subtitle is not None:
            title = '{} {}'.format(
                subtitle,
                title)
        leg.set_title(title)
        return leg

    @staticmethod
    def update_frame_1d(i, scats, lines, bodies, coordinates, dims):
        xdim = dims[0]
        yi = 0
        yj = np.zeros(5)
        for scat, line, body, coords in zip(scats, lines, bodies, coordinates):
            xi = coords[xdim, i]
            # xi = coords[i, xdim]
            scat.set_offsets((xi, yi))
            if i > 5:
                xj = coords[xdim, i-5 : i]
                line.set_data((xj, yj))
        return scats, lines

    @staticmethod
    def update_frame_2d(i, scats, lines, bodies, coordinates, dims):
        xdim, ydim = dims[0], dims[1]
        for scat, line, body, coords in zip(scats, lines, bodies, coordinates):
            xi = coords[xdim, i]
            yi = coords[ydim, i]
            # xi = coords[i, xdim]
            # yi = coords[i, ydim]
            scat.set_offsets((xi, yi))
            if i > 5:
                xj = coords[xdim, i-5 : i]
                yj = coords[ydim, i-5 : i]
                line.set_data((xj, yj))
        return scats, lines

    @staticmethod
    def update_frame_3d(i, scats, lines, bodies, coordinates, dims):
        xdim, ydim, zdim = dims[0], dims[1], dims[2]
        for scat, line, body, coords in zip(scats, lines, bodies, coordinates):
            xi = coords[xdim, i]
            yi = coords[ydim, i]
            zi = coords[zdim, i]
            # xi = coords[i, xdim]
            # yi = coords[i, ydim]
            # zi = coords[i, zdim]
            scat.set_offsets((xi, yi))
            scat.set_3d_properties(zi, zdir='z')
            if i > 5:
                xj = coords[xdim, i-5 : i]
                yj = coords[ydim, i-5 : i]
                zj = coords[zdim, i-5 : i]
                line.set_data((xj, yj))
                line.set_3d_properties(zj)
        return scats, lines

    def configure_frame_update(self, dims):
        if dims.size == 1:
            self._update_frame_ndim = lambda *args, **kwargs : self.update_frame_1d(*args, **kwargs)
        elif dims.size == 2:
            self._update_frame_ndim = lambda *args, **kwargs : self.update_frame_2d(*args, **kwargs)
        elif dims.size == 3:
            self._update_frame_ndim = lambda *args, **kwargs : self.update_frame_3d(*args, **kwargs)
        else:
            raise ValueError("invalid dims.size: {}".format(dims.size))

    def update_frame(self, i, scats, lines, bodies, coordinates, dims, leg, subtitle, time_coordinates, time_units):
        scats, lines = self.update_frame_ndim(
            i=i,
            scats=scats,
            lines=lines,
            bodies=bodies,
            coordinates=coordinates,
            dims=dims)
        leg = self.update_legend_counter(
            i=i,
            leg=leg,
            subtitle=subtitle,
            time_coordinates=time_coordinates,
            time_units=time_units)

    def display_animation(self, fig, anim, fps=None, savename=None, extension=None):
        if savename is None:
            plt.show()
        elif isinstance(savename, str):
            if extension is None:
                extension = '.mp4'
            savepath = '{}{}{}'.format(self.savedir, savename, extension)
            anim.save(
                savepath,
                fps=fps)
        else:
            raise ValueError("invalid type(savename): {}".format(type(savename)))
        plt.close(fig)

class VisualConfiguration(AnimationConfiguration):

    def __init__(self, savedir=None, ticksize=7, labelsize=8, textsize=5, titlesize=9):
        super().__init__(savedir=savedir)
        self.ticksize = ticksize
        self.labelsize = labelsize
        self.textsize = textsize
        self.titlesize = titlesize

    @staticmethod
    def get_empty_handle(ax):
        empty_handle = ax.scatter([np.nan], [np.nan], color='none', alpha=0)
        return empty_handle

    @staticmethod
    def get_number_of_legend_columns(labels):
        if isinstance(labels, int):
            n = labels
        else:
            n = len(labels)
        if n > 2:
            if n % 3 == 0:
                ncol = 3
            else:
                ncol = n // 2
        else:
            ncol = n
        return ncol

    def update_legend_design(self, leg, title=None, textcolor=None, facecolor=None, edgecolor=None, borderaxespad=None):
        if title:
            leg.set_title(title, prop={'size': self.labelsize, 'weight' : 'semibold'})
            if textcolor:
                leg.get_title().set_color(textcolor)
            # leg.get_title().set_ha("center")
        leg._legend_box.align = "center"
        frame = leg.get_frame()
        if facecolor:
            frame.set_facecolor(facecolor)
        if edgecolor:
            frame.set_edgecolor(edgecolor)
        if textcolor:
            for text in leg.get_texts():
                text.set_color(textcolor)
        return leg

    def subview_legend(self, fig, ax, handles, labels, title=None, bottom=0.2, textcolor='darkorange', facecolor='k', edgecolor='steelblue'):
        if len(labels) == 1:
            ncol = 3
            empty_handle = self.get_empty_handle(ax)
            handles = [empty_handle, handles[0], empty_handle]
            labels = ['  ', labels[0], '  ']
        else:
            ncol = self.get_number_of_legend_columns(labels)
        fig.subplots_adjust(bottom=bottom)
        leg = fig.legend(
            handles=handles,
            labels=labels,
            ncol=ncol,
            loc='lower center',
            mode='expand',
            borderaxespad=0.1,
            fontsize=self.labelsize)
        leg = self.update_legend_design(
            leg=leg,
            title=title,
            textcolor=textcolor,
            facecolor=facecolor,
            edgecolor=edgecolor)
        return leg

    def display_image(self, fig, savename=None, dpi=800, bbox_inches='tight', pad_inches=0.1, extension='.png', **kwargs):
        if savename is None:
            plt.show()
        elif isinstance(savename, str):
            if self.savedir is None:
                raise ValueError("cannot save plot; self.savedir is None")
            savepath = '{}{}{}'.format(self.savedir, savename, extension)
            fig.savefig(savepath, dpi=dpi, bbox_inches=bbox_inches, pad_inches=pad_inches, **kwargs)
        else:
            raise ValueError("invalid type(savename): {}".format(type(savename)))
        plt.close(fig)












#
