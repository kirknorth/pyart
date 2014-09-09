"""
pyart.graph.griddisplay
=======================

A class for plotting grid objects onto model (Cartesian) axes.
"""


import numpy as np
import matplotlib.pyplot as plt


class GridDisplay:
    """
    """

    def __init__(self, grid, debug=False):
        """ Initialize graphing object. """

        # Populate fields attribute
        self.fields = grid.fields

        # Set source or instrument attribute
        if 'radar_0_instrument_name' in grid.metadata:
            self.source = str(grid.metadata['radar_0_instrument_name'])
        else:
            self.source = 'Model'

        # Populate axes attributes
        self.start_time = datetime_utils.datetime_from_grid(grid)
        self.x_disp = grid.axes['x_disp']['data']
        self.y_disp = grid.axes['y_disp']['data']
        self.z_disp = grid.axes['z_disp']['data']
        self.lat = grid.axes['lat']['data'][0]
        self.lon = grid.axes['lon']['data'][0]

        # Populate graph attributes
        self.plots = []
        self.plot_vars = []

    def plot_xy_slice(self, field, level=0, vmin=None, vmax=None, cmap='jet',
                      title=None, title_flag=True, colorbar_label=None,
                      colorbar_flag=True, orientation='vertical',
                      axislabels=(None, None), axislabels_flag=True, fig=None,
                      ax=None):
        """
        Plot a horizontal cross section of the grid.
        """

        # Parse plot parameters
        ax = self._parse_ax(ax)
        vmin, vmax = self._parse_vmin_vmax(field, vmin, vmax)

        # Create pcolormesh and add it to graph object
        qm = ax.pcolormesh(self.x_disp / 1000.0, self.y_disp / 1000.0,
                           self.fields[field]['data'][level],
                           vmin=vmin, vmax=vmax, cmap=cmap)
        self.plots.append(qm)
        self.plot_vars.append(field)

        # Set axis labels
        if axislabels_flag:
            self._label_axes_xy(axislabels, ax)

        # Plot color bar
        if colorbar_flag:
            self.plot_colorbar(field, mappable=qm, orientation=orientation,
                               label=colorbar_label, fig=fig, ax=ax)

        # Set title
        if title_flag:
            self._set_xy_title(field, level, title, ax)

        return

    def plot_colorbar(self, field, mappable=None, orientation='vertical',
                      label=None, fig=None, ax=None, cax=None):
        """
        Plot a color bar to a specified axis.

        Parameters
        ----------
        mappable : ScalarMappable
            None will use the last mappable plotted.
        orientation : 'vertical' or 'horizontal'
            Color bar orientation on plot.
        label : str
            Color bar label. None will use 'long_name' ['units'] for the
            last field plotted or an empty string if the field does not have
            these keys.
        fig : Figure
            A figure instance.
        ax : Axes
            An axes instance.
        """

        # Parse plot parameters
        mp = self._parse_mappable(mappable)

        # Create color bar instance
        cb = plt.colorbar(mappable=mp, orientation=orientation, ax=ax, cax=cax)

        # Set color bar label
        self._set_colorbar_label(field, label, cb)

        return

    def generate_xy_title(self, field, level):
        """
        Generate a title for a horizontal cross section plot.

        Parameters
        ----------
        field : str
            The field to be plotted.
        level : int
            The height level to be plotted.

        Return
        ------
        title : str
            The title for the plot.
        """

        date = self.start_time.isoformat() + 'Z'
        height = '%.2f km' % self.z_disp[level] / 1000.0
        name = self.fields[field]['long_name']

        return '%s %.2f km %s \n %s' % (self.source, height, date, name)

    def generate_colorbar_label(self, field):
        """
        Generate the color bar label for the plot.
        """

        name = self.fields[feild]['long_name']
        units = self.fields[field]['units']

        return '%s [%s]' % (name, units)

    def _label_xy_axes(self, labels, ax):
        """
        Label the axes for a horizontal cross section plot.
        """

        ax = self._parse_ax(ax)
        xlabel, ylabel = labels

        if xlabel is None:
            ax.set_xlabel('Eastward Distance from Origin [km]')
        else:
            ax.set_xlabel(xlabel)
        if ylabel is None:
            ax.set_ylabel('Northward Distance from Origin [km]')
        else:
            ax.set_ylabel(ylabel)

    def _set_xy_title(self, field, level, title, ax):
        """
        Label the title for a horizontal cross section plot.
        """

        ax = self._parse_ax(ax)

        if title is None:
            ax.set_title(self.generate_xy_title(field, level))
        else:
            ax.set_title(title)

    def _set_colorbar_label(self, field, label, cb):
        """
        """

        if label is None:
            cb.set_label(self.generate_colorbar_label(field))
        else:
            cb.set_label(label)

    def _parse_vmin_vmax(self, field, vmin, vmax):
        """
        Parse the vmin and vmax plot parameters. If these parameters are not
        supplied and the field does not have 'valid_min' and 'valid_max'
        metadata, then vmin and vmax will default to reflectivity values.
        """

        if vmin is None:
            if 'valid_min' in self.fields[field]:
                vmin = self.fields[field]['valid_min']
            else:
                vmin = -8
        if vmax is None:
            if 'valid_max' in self.fields[field]:
                vmax = self.fields[field]['valid_max']
            else:
                vmax = 64

        return vmin, vmax

    def _parse_ax(self, ax):
        """
        Parse the plot axis parameter.
        """

        if ax is None:
            ax = plt.gca()

        return ax

    def _parse_mappable(self, mappable):
        """
        Parse the plot ScalarMappable parameter.
        """

        if mappable is None:
            if self.plots:
                raise ValueError('No mappables (plots) found!')
            else:
                mp = self.plots[-1]

        return mp
