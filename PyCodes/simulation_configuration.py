from collections import OrderedDict
from diffeq_integrators import *
from visual_configuration import *
from digital_signal_processing import *

class SimulationConfiguration(VisualConfiguration):

    def __init__(self, bodies, gravitational_constant=6.67e-11, name='', savedir=None, epsilon=1e-6):
        """

        """
        super().__init__(savedir)
        self.vector_quantities = (
            'position',
            'velocity',
            'acceleration',
            'momentum',
            'angular momentum',
            'force',
            # 'eccentricity',
            )
        self.scalar_quantities = (
            'potential energy',
            'kinetic energy',
            'total energy',
            'lagrangian')
        self.epsilon = epsilon
        self._ode_system_solver = {
            'velocity-verlet': VerletLeapFrogIntegrator(bodies, gravitational_constant, name),
            'yoshida' : YoshidaLeapFrogIntegrator(bodies, gravitational_constant, name),
            'euler-cromer' : ModifiedEulerIntegrator(bodies, gravitational_constant, name),
            'rk-4' : RungeKuttaFourIntegrator(bodies, gravitational_constant, name),
            'adams-bashforth-moulton' : AdamsBashforthMoultonIntegrator(bodies, gravitational_constant, name)}
        self._bodies = OrderedDict()
        self._body_names = None

    @property
    def ode_system_solver(self):
        return self._ode_system_solver

    @property
    def bodies(self):
        return self._bodies

    @property
    def body_names(self):
        return self._body_names

    def update_bodies(self, cls):
        names = []
        for i in range(cls.nbodies):
            body = cls.bodies[i]
            r = cls.solution['position'][i, :, :]
            v = cls.solution['velocity'][i, :, :]
            a = cls.solution['acceleration'][i, :, :]
            p = body.mass * v
            h = np.cross(r, v, axis=0)
            rmag = cls.get_scalar(r, axis=0)
            # ecc = np.cross(v, h, axis=0) / cls.reduced_mass - r / rmag
            e_k = np.sum(p * v, axis=0) / 2
            e_p = 0
            for j in range(cls.nbodies):
                if i != j:
                    displacement = r - cls.solution['position'][j, :, :]
                    distance = cls.get_scalar(displacement, axis=0)
                    # distance[np.isnan(distance)] = self.epsilon
                    distance[distance == 0] = self.epsilon
                    e_p -= cls.standard_mu[i] * cls.bodies[j].mass / distance
            e_h = e_k + e_p
            self._bodies[body.name] = {
                'name' : body.name,
                'facecolor' : body.facecolor,
                'mass' : body.mass,
                'standard gravitational parameter' : cls.standard_mu[i],
                'position' : r,
                'velocity' : v,
                'acceleration' : a,
                'momentum' : p,
                'angular momentum' : body.mass * h,
                'force' : body.mass * a,
                # 'eccentricity' : ecc,
                'lagrangian' : e_k - e_p,
                'potential energy' : e_p,
                'kinetic energy' : e_k,
                'total energy' : e_h}
            names.append(body.name)
        self._body_names = names

    def update_system_center_of_mass(self, cls):
        r = np.zeros((cls.ndim, cls.t.size))
        v = np.zeros((cls.ndim, cls.t.size))
        a = np.zeros((cls.ndim, cls.t.size))
        e_p, e_k = 0, 0
        for body, d in self.bodies.items():
            r += d['position'] * d['mass'] / cls.total_mass
            v += d['velocity'] * d['mass'] / cls.total_mass
            a += d['acceleration'] * d['mass'] / cls.total_mass
            e_p += d['potential energy']
            e_k += d['kinetic energy']
        p = v * cls.reduced_mass # ... cls.total_mass
        h = np.cross(r, v, axis=0)
        rmag = cls.get_scalar(r, axis=0)
        # ecc = np.cross(v, h, axis=0) / cls.reduced_mass - r / rmag
        name = 'Center-of-Mass'
        e_h = e_k + e_p
        com_system = {
            'name' : name,
            'facecolor' : 'k',
            'mass' : cls.total_mass,
            'position' : r,
            'velocity' : v,
            'acceleration' : a,
            'momentum' : p,
            'angular momentum' : cls.reduced_mass * h,
            'force' : cls.reduced_mass * a, # ...cls.total_mass
            # 'eccentricity' : ecc,
            'lagrangian' : e_k - e_p,
            'potential energy' : e_p,
            'kinetic energy' : e_k,
            'total energy' : e_h}
        self._bodies[name] = com_system
        self._body_names.append(name)
        self._body_names = np.array(self._body_names)

    def update_orbital_periods(self, method, f_window=None):
        for name, body in self.bodies.items():
            if body['name'] != 'Center-of-Mass':
                scalars = self.get_scalar_coordinates(
                    quantity='position',
                    body=body,
                    ref_frame='Center-of-Mass',
                    scale=1)
                dt = np.mean(np.diff(self.ode_system_solver.t))
                dsp = DigitalSignalProcessing(
                    scalars=scalars,
                    t=self.ode_system_solver.t,
                    dt=dt,
                    sampling_rate=1/dt,
                    f_window=f_window)
                period = dsp.get_period(
                    method=method)
                body['orbital period'] = {
                    'method' : method,
                    'value' : period}

    def select_bodies(self, names=None, selection_criteria='inclusive'):
        ## get bodies by name
        if names is None:
            result = [dict(body) for name, body in self.bodies.items() if name != 'Center-of-Mass']
        else:
            if isinstance(names, str):
                names = [names]
            elif not isinstance(names, (tuple, list, np.ndarray)):
                raise ValueError("invalid type(names): {}".format(type(names)))
            ## apply selection criteria
            if selection_criteria == 'exclusive':
                result = [dict(body) for name, body in self.bodies.items() if name not in names]
            elif selection_criteria == 'inclusive':
                result = [dict(body) for name, body in self.bodies.items() if name in names]
            else:
                raise ValueError("invalid selection_criteria: {}".format(selection_criteria))
            ## verify selected bodies
            if len(result) == 0:
                raise ValueError("there are no bodies to return")
        return result

    def select_vector_dimensions(self, dims):
        if dims is None:
            dims = np.arange(self.ode_system_solver.ndim, dtype=int)
        elif isinstance(dims, int):
            dims = np.array([dims], dtype=int)
        elif isinstance(dims, (tuple, list, np.ndarray)):
            if not isinstance(dims, np.ndarray):
                dims = np.array(dims)
        else:
            raise ValueError("invalid type(dims): {}".format(type(dims)))
        if dims.size > 3: # dim-0, dim-1, dim-2
            raise ValueError("number of dimensions cannot exceed 3")
        return dims

    def get_vector_coordinates(self, quantity, body, ref_frame=None, scale=1):
        if quantity not in self.vector_quantities:
            raise ValueError("{} is not a valid vector quantity: {}".format(quantity))
        coordinates = np.array(body[quantity])
        if ref_frame:
            ref_body = self.bodies[ref_frame]
            coordinates -= ref_body[quantity]
        return coordinates * scale

    def get_scalar_coordinates(self, quantity, body, ref_frame=None, scale=1):
        if quantity in self.vector_quantities:
            vector = self.get_vector_coordinates(quantity, body, ref_frame, scale)
            scalar = self.ode_system_solver.get_scalar(vector, axis=0)
        elif quantity in self.scalar_quantities:
            # scalar = body[quantity]
            # if ref_frame:
            #     ref_body = self.bodies[ref_frame]
            #     scalar -= ... # ref_body[quantity]
            # scalar *= scale
            if 'energy' in quantity:
                if quantity == 'total energy':
                    e_kin = self.get_scalar_coordinates(
                        quantity='kinetic energy',
                        body=body,
                        ref_frame=ref_frame,
                        scale=scale)
                    e_pot = self.get_scalar_coordinates(
                        quantity='potential energy',
                        body=body,
                        ref_frame=ref_frame,
                        scale=scale)
                    scalar = e_kin + e_pot
                else:
                    if quantity == 'kinetic energy':
                        scalar = self.get_kinetic_energy_by_frame(
                            body=body,
                            ref_frame=ref_frame)
                    else: # quantity == 'potential energy'
                        names = ['Center-of-Mass']
                        if body['name'] != 'Center-of-Mass':
                            names.append(body['name'])
                        other_bodies = self.select_bodies(
                            names=names,
                            selection_criteria='exclusive')
                        scalar = self.get_potential_energy_by_frame(
                            body=body,
                            other_bodies=other_bodies,
                            ref_frame=ref_frame)
                    scalar *= scale
            elif quantity == 'lagrangian':
                e_kin = self.get_scalar_coordinates(
                    quantity='kinetic energy',
                    body=body,
                    ref_frame=ref_frame,
                    scale=scale)
                e_pot = self.get_scalar_coordinates(
                    quantity='potential energy',
                    body=body,
                    ref_frame=ref_frame,
                    scale=scale)
                scalar = e_kin - e_pot
            else:
                raise ValueError("not yet implemented; quantity={}".format(quantity))
        else:
            raise ValueError("{} is not a vector or scalar quantity: {}".format(quantity))
        return scalar

    def get_kinetic_energy_by_frame(self, body, ref_frame=None):
        vel = self.get_vector_coordinates(
            quantity='velocity',
            body=body,
            ref_frame=ref_frame,
            scale=1)
        # print(body.name, vel, vel.shape)
        # scalar = body.mass * np.dot(vel, vel)
        scalar = body['mass'] * self.ode_system_solver.get_scalar(
            vector=np.square(vel),
            axis=0)
        return scalar

    def get_potential_energy_by_frame(self, body, other_bodies, ref_frame=None):
        scalar = 0
        for other_body in other_bodies:
            pos = self.get_vector_coordinates(
                quantity='position',
                body=body,
                ref_frame=ref_frame,
                scale=1)
            other_pos = self.get_vector_coordinates(
                quantity='position',
                body=other_body,
                ref_frame=ref_frame,
                scale=1)
            displacement = pos - other_pos
            distance = self.ode_system_solver.get_scalar(displacement, axis=0)
            distance[distance == 0] = self.epsilon # distance[np.isnan(distance)] = self.epsilon
            scalar -= self.ode_system_solver.gravitational_constant * body['mass'] * other_body['mass'] / distance
        return scalar

    def view_vector_quantity(self, quantity, position_units=None, time_units=None, mass_units=None, energy_units=None, ref_frame='Center-of-Mass', names=None, selection_criteria='inclusive', subtitle=None, dims=None, elev=None, azim=None, save=False, **kwargs):
        ## verify quantity
        if quantity not in self.vector_quantities:
            raise ValueError("invalid quantity: {}".format(quantity))
        ## get time-scale
        if time_units is None:
            time_units = self.ode_system_solver.units['time']
        time_scale = self.get_scale('time', itime_units=self.ode_system_solver.units['time'], jtime_units=time_units)
        ## get duration
        duration = self.ode_system_solver.duration * time_scale
        if int(duration) == float(duration):
            duration = int(duration)
        ## get quantity_scale
        quantity_scale = self.get_scale(quantity,
            iposition_units=self.ode_system_solver.units['position'], jposition_units=position_units,
            itime_units=self.ode_system_solver.units['time'], jtime_units=time_units,
            imass_units=self.ode_system_solver.units['mass'], jmass_units=mass_units,
            ienergy_units=self.ode_system_solver.units['energy'], jenergy_units=energy_units)
        ## get bodies
        bodies = self.select_bodies(names, selection_criteria)
        ## get viewing dimensions
        dims = self.select_vector_dimensions(dims)
        ## initialize plot
        if dims.size == 3:
            fig = plt.figure(**kwargs)
            ax = fig.add_subplot(1, 1, 1, projection='3d')
        else:
            fig, ax = plt.subplots(**kwargs)
        ## view per body / center-of-mass
        for body in bodies:
            vector_coordinates = self.get_vector_coordinates(quantity, body, ref_frame, quantity_scale)
            if dims.size == 3:
                x, y, z = vector_coordinates[dims[0], :], vector_coordinates[dims[1], :], vector_coordinates[dims[2], :]
                ax.scatter(x, y, z,
                    marker='.',
                    label=body['name'],
                    color=body['facecolor'],
                    s=2)
            else:
                x = vector_coordinates[dims[0], :]
                if dims.size == 1:
                    y = np.zeros(x.shape)
                else:
                    y = vector_coordinates[dims[1], :]
                ax.scatter(x, y,
                    marker='.',
                    label=body['name'],
                    color=body['facecolor'],
                    s=2)
        ## modify angular viewing position
        if elev is not None:
            ax.elev = elev
        if azim is not None:
            ax.azim = azim
        ## axis labels
        axis_label = self.get_label(quantity, position_units, time_units, mass_units, energy_units, include_quantity=True)
        ax.set_xlabel(axis_label, fontsize=self.labelsize)
        ax.set_ylabel(axis_label, fontsize=self.labelsize)
        ## axis title
        time_label = self.get_label('time', time_units=time_units)
        title = '${}$-Body Simulation (${}$ {})'.format(self.ode_system_solver.nbodies, duration, time_label)
        if ref_frame:
            title += '\nvia {} Frame'.format(ref_frame)
        ax.set_title(title, fontsize=self.titlesize)
        ## axis ticks
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
        if dims.size == 3:
            ax.set_zlabel(axis_label, fontsize=self.labelsize)
            ax.zaxis.set_minor_locator(ticker.AutoMinorLocator())
            for axis in ('x', 'y', 'z'):
                ax.tick_params(axis=axis, labelsize=self.ticksize)
        else:
            ax.tick_params(axis='both', labelsize=self.ticksize)
            ax.grid(color='k', linestyle=':', alpha=0.3)
        ## equalize axis-scaling
        try:
            ax.set_aspect('equal')
        except NotImplementedError:
            pass
        ## legend
        handles, labels = ax.get_legend_handles_labels()
        if subtitle is not None:
            if not isinstance(subtitle, str):
                raise ValueError("invalid type(subtitle): {}".format(type(subtitle)))
            # subtitle = 'Solution via {} Method'.format(self.ode_system_solver.method)
        self.subview_legend(fig, ax, handles, labels, title=subtitle)
        ## save or show figure
        if save:

            savename = '' if subtitle is None else '{}-'.format(subtitle.replace(' ', '_'))
            savename += '{}-'.format(self.ode_system_solver.method.replace('-', '_').replace(' ', '_'))
            savename += 'vector_{}-dims_{}-'.format(quantity.replace(' ', '_'), '_'.join(dims.astype(str)))
            savename += 'ref_{}'.format(ref_frame.replace('-', '_').replace(' ', '_'))
            savename = savename.replace(' ', '-')

            # savename = 'vector_{}_dims-{}'.format(quantity.replace(' ', '_'), '_'.join(dims.astype(str)))
            # if ref_frame:
            #     savename += '_ref-{}'.format(ref_frame.replace(' ', '_'))
            # if subtitle:
            #     savename += '_{}'.format(subtitle.replace(' ', '_'))
            # savename += '_{}'.format(self.ode_system_solver.method.replace(' ', '_'))

        else:
            savename = None
        self.display_image(fig, savename)

    def view_scalar_quantity(self, quantity, position_units=None, time_units=None, mass_units=None, energy_units=None, ref_frame='Center-of-Mass', names=None, selection_criteria='inclusive', subtitle=None, save=False, **kwargs):
        if quantity not in list(self.vector_quantities) + list(self.scalar_quantities) + ['energy']:
            raise ValueError("invalid quantity: {}".format(quantity))
        ## get time-scale
        if time_units is None:
            time_units = self.ode_system_solver.units['time']
        time_scale = self.get_scale('time', itime_units=self.ode_system_solver.units['time'], jtime_units=time_units)
        ## get duration
        duration = self.ode_system_solver.duration * time_scale
        if int(duration) == float(duration):
            duration = int(duration)
        ## get quantity_scale
        quantity_scale = self.get_scale(quantity,
            iposition_units=self.ode_system_solver.units['position'], jposition_units=position_units,
            itime_units=self.ode_system_solver.units['time'], jtime_units=time_units,
            imass_units=self.ode_system_solver.units['mass'], jmass_units=mass_units,
            ienergy_units=self.ode_system_solver.units['energy'], jenergy_units=energy_units)
        ## get bodies
        bodies = self.select_bodies(names, selection_criteria)
        ## initialize plot
        x = self.ode_system_solver.t * time_scale
        fig, ax = plt.subplots(**kwargs)
        ## view per body / center-of-mass
        if quantity == 'energy': ## 'potential energy', 'kinetic energy', 'total energy'
            labels = np.char.title(self.scalar_quantities)
            if len(bodies) == 1:
                markers = ('.', '.', '.')
                facecolors = ('darkorange', 'steelblue', 'mediumviolet')
            else:
                markers = ('.', 'x', '_')
                facecolors = tuple([body['facecolor'] for body in bodies])
            for body in bodies:
                for qty, marker, label, facecolor in zip(self.scalar_quantities, markers, labels, facecolors):
                    y = self.get_scalar_coordinates(qty, body, ref_frame, quantity_scale)
                    ax.scatter(x, y,
                        marker=marker,
                        label=r'${%s}_{%s}$' % (label, body['name']),
                        color=facecolor,
                        s=2)
        else:
            for body in bodies:
                y = self.get_scalar_coordinates(quantity, body, ref_frame, quantity_scale)
                ax.scatter(x, y,
                    marker='.',
                    label=body['name'],
                    color=body['facecolor'],
                    s=2)
        ## axis labels
        xlabel = self.get_label('time', position_units, time_units, mass_units, energy_units, include_quantity=True)
        ylabel = self.get_label(quantity, position_units, time_units, mass_units, energy_units, include_quantity=True)
        ax.set_xlabel(xlabel, fontsize=self.labelsize)
        ax.set_ylabel(ylabel, fontsize=self.labelsize)
        ## axis title
        time_label = self.get_label('time', time_units=time_units)
        title = '${}$-Body Simulation (${}$ {})'.format(self.ode_system_solver.nbodies, duration, time_label)
        if ref_frame:
            title += '\nvia {} Frame'.format(ref_frame)
        ax.set_title(title, fontsize=self.titlesize)
        ## axis ticks
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax.tick_params(axis='both', labelsize=self.ticksize)
        ax.grid(color='k', linestyle=':', alpha=0.3)
        ## legend
        handles, labels = ax.get_legend_handles_labels()
        if subtitle is not None:
            if not isinstance(subtitle, str):
                raise ValueError("invalid type(subtitle): {}".format(type(subtitle)))
        self.subview_legend(fig, ax, handles, labels, title=subtitle)
        ## save or show figure
        if save:

            savename = '' if subtitle is None else '{}-'.format(subtitle.replace(' ', '_'))
            savename += '{}-'.format(self.ode_system_solver.method.replace('-', '_').replace(' ', '_'))
            savename += 'scalar_{}-'.format(quantity.replace(' ', '_'))
            savename += 'ref_{}'.format(ref_frame.replace('-', '_').replace(' ', '_'))
            savename = savename.replace(' ', '-')

            # savename = 'scalar_{}'.format(quantity.replace(' ', '_'))
            # if ref_frame:
            #     savename += '_ref-{}'.format(ref_frame.replace(' ', '_'))
            # if subtitle:
            #     savename += '_{}'.format(subtitle.replace(' ', '_'))
            # savename += '_{}'.format(self.ode_system_solver.method.replace(' ', '_'))

        else:
            savename = None
        self.display_image(fig, savename)

    def view_vector_quantity_animation(self, quantity, position_units=None, time_units=None, mass_units=None, energy_units=None, ref_frame='Center-of-Mass', names=None, selection_criteria='inclusive', subtitle=None, dims=None, elev=None, azim=None, save=False, extension=None, fps=None, **kwargs):
        ## verify quantity
        if quantity not in self.vector_quantities:
            raise ValueError("invalid quantity: {}".format(quantity))
        ## get time-scale
        if time_units is None:
            time_units = self.ode_system_solver.units['time']
        time_scale = self.get_scale('time', itime_units=self.ode_system_solver.units['time'], jtime_units=time_units)
        ## get duration
        duration = self.ode_system_solver.duration * time_scale
        if int(duration) == float(duration):
            duration = int(duration)
        ## get quantity_scale
        quantity_scale = self.get_scale(quantity,
            iposition_units=self.ode_system_solver.units['position'], jposition_units=position_units,
            itime_units=self.ode_system_solver.units['time'], jtime_units=time_units,
            imass_units=self.ode_system_solver.units['mass'], jmass_units=mass_units,
            ienergy_units=self.ode_system_solver.units['energy'], jenergy_units=energy_units)
        ## get bodies
        bodies = self.select_bodies(names, selection_criteria)
        ## get viewing dimensions
        dims = self.select_vector_dimensions(dims)
        ## initialize plot
        if dims.size == 3:
            fig = plt.figure(**kwargs)
            ax = fig.add_subplot(1, 1, 1, projection='3d')
        else:
            fig, ax = plt.subplots(**kwargs)
        ## axis labels
        axis_label = self.get_label(quantity, position_units, time_units, mass_units, energy_units, include_quantity=True)
        ax.set_xlabel(axis_label, fontsize=self.labelsize)
        ax.set_ylabel(axis_label, fontsize=self.labelsize)
        if dims.size == 3:
            ax.set_zlabel(axis_label, fontsize=self.labelsize)
        ## axis title
        time_label = self.get_label('time', time_units=time_units)
        title = '${}$-Body Simulation (${}$ {})'.format(self.ode_system_solver.nbodies, duration, time_label)
        if ref_frame:
            title += '\nvia {} Frame'.format(ref_frame)
        ax.set_title(title, fontsize=self.titlesize)
        ## axis limits
        lower_bounds, upper_bounds = [], []
        for body in bodies:
            vector_coordinates = self.get_vector_coordinates(quantity, body, ref_frame, quantity_scale)
            vmin, vmax = np.nanmin(vector_coordinates, axis=0), np.nanmax(vector_coordinates, axis=0)
            # vmin, vmax = np.nanmin(vector_coordinates, axis=1), np.nanmax(vector_coordinates, axis=1)
            # vmin, vmax = np.nanmin(vector_coordinates, axis=-1), np.nanmax(vector_coordinates, axis=-1)
            lower_bounds.append(vmin)
            upper_bounds.append(vmax)
        lower_bounds, upper_bounds = np.array(lower_bounds), np.array(upper_bounds)
        lower_bounds = np.nanmin(lower_bounds, axis=0) # * quantity_scale
        upper_bounds = np.nanmax(upper_bounds, axis=0) # * quantity_scale
        if dims.size > 1:
            lower_bounds, upper_bounds = lower_bounds[dims], upper_bounds[dims]
        else:
            lower_bounds, upper_bounds = np.array([lower_bounds[0], -1]), np.array([upper_bounds[0], 1])
        lower_bounds = [lb * 1.25 if lb < 0 else lb * 0.75 if lb > 0 else -1 for lb in lower_bounds]
        upper_bounds = [ub * 1.25 if ub > 0 else ub * 0.75 if ub < 0 else 1 for ub in upper_bounds]
        try:
            f_limits = (ax.set_xlim, ax.set_ylim, ax.set_zlim)
        except AttributeError:
            f_limits = (ax.set_xlim, ax.set_ylim)
        for f, lb, ub in zip(f_limits, lower_bounds, upper_bounds):
            f(lb, ub)
        ## axis ticks
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
        if dims.size == 3:
            ax.zaxis.set_minor_locator(ticker.AutoMinorLocator())
            for axis in ('x', 'y', 'z'):
                ax.tick_params(axis=axis, labelsize=self.ticksize)
        else:
            ax.tick_params(axis='both', labelsize=self.ticksize)
            ax.grid(color='k', linestyle=':', alpha=0.3)
        ## equalize axis-scaling
        try:
            ax.set_aspect('equal')
        except NotImplementedError:
            pass

        ## modify angular viewing position
        if elev is not None:
            ...
        else:
            ...
        if azim is not None:
            ...
        else:
            ...

        ## initialize zeroth frame
        coordinates, scats, lines = [], [], []
        for body in bodies:
            vector_coordinates = self.get_vector_coordinates(quantity, body, ref_frame, quantity_scale)
            coordinates.append(vector_coordinates)
            xdim = dims[0]
            xi = vector_coordinates[xdim, 0]
            if dims.size == 3:
                ydim, zdim = dims[1], dims[2]
                yi = vector_coordinates[ydim, 0]
                zi = vector_coordinates[zdim, 0]
                scat = ax.scatter(
                    xi,
                    yi,
                    zi,
                    marker='.',
                    color=body['facecolor'],
                    label=body['name'],
                    s=2)
                line, = ax.plot(
                    [],
                    [],
                    [],
                    color=body['facecolor'],
                    alpha=0.3)
            else:
                if dims.size == 1:
                    yi = 0
                else:
                    ydim = dims[1]
                    yi = vector_coordinates[ydim, 0]
                scat = ax.scatter(
                    xi,
                    yi,
                    marker='.',
                    color=body['facecolor'],
                    label=body['name'],
                    s=2)
                line, = ax.plot(
                    [],
                    [],
                    color=body['facecolor'],
                    alpha=0.3)
            scats.append(scat)
            lines.append(line)
        coordinates = np.array(coordinates)

        ## initialize legend
        if subtitle is not None:
            if not isinstance(subtitle, str):
                raise ValueError("invalid type(subtitle): {}".format(type(subtitle)))
        handles, labels = ax.get_legend_handles_labels()
        if subtitle is not None:
            if not isinstance(subtitle, str):
                raise ValueError("invalid type(subtitle): {}".format(type(subtitle)))
            # subtitle = 'Solution via {} Method'.format(self.ode_system_solver.method)
        leg = self.subview_legend(
            fig=fig,
            ax=ax,
            handles=handles,
            labels=labels,
            title=subtitle,
            bottom=0.325,
            textcolor='k',
            facecolor='silver',
            edgecolor='steelblue')

        ## initialize animation
        self.configure_frame_update(dims=dims)
        time_coordinates = time_scale * self.ode_system_solver.t
        args = (scats, lines, bodies, coordinates, dims, leg, subtitle, time_coordinates, time_units)
        anim = FuncAnimation(
            fig,
            self.update_frame,
            frames=self.ode_system_solver.t.size,
            fargs=args)
        ## save or show figure
        if save:

            savename = '' if subtitle is None else '{}-'.format(subtitle.replace(' ', '_'))
            savename += '{}-'.format(self.ode_system_solver.method.replace('-', '_').replace(' ', '_'))
            savename += 'vector_{}-dims_{}-'.format(quantity.replace(' ', '_'), '_'.join(dims.astype(str)))
            savename += 'ref_{}'.format(ref_frame.replace('-', '_').replace(' ', '_'))
            savename = savename.replace(' ', '-')

            # savename = 'vector_{}_dims-{}'.format(quantity.replace(' ', '_'), '_'.join(dims.astype(str)))
            # if ref_frame:
            #     savename += '_ref-{}'.format(ref_frame.replace(' ', '_'))
            # if subtitle:
            #     savename += '_{}'.format(subtitle.replace(' ', '_'))
            # savename += '_{}'.format(self.ode_system_solver.method.replace(' ', '_'))

        else:
            savename = None
        self.display_animation(
            fig=fig,
            anim=anim,
            fps=fps,
            savename=savename,
            extension=extension)

    def view_orbital_periods_table(self, true_periods, subtitle=None, save=False, **kwargs):
        ## initialize figure and axes
        fig, ax = plt.subplots(**kwargs)
        # fig.patch.set_visible(False)
        ax.axis('off')
        ax.axis('tight')
        ## initialize table parameters
        column_labels = ['Method', 'Orbital Period [s] (sim)', 'Orbital Period [s] (true)', 'Relative Error [%]']
        row_labels, cell_text, cell_colors = [], [], []
        for i, (name, p1) in enumerate(true_periods.items()):
            bodies = self.select_bodies(
                names=name,
                selection_criteria='inclusive')
            if len(bodies) != 1:
                raise ValueError("search criteria should only find one body per search name")
            body = bodies[0]
            p2 = body['orbital period']['value']
            method = body['orbital period']['method']
            rel_error = abs(p2 - p1) / p1 * 100
            row_labels.append(name)
            row_values = [
                method,
                '${:.3}$'.format(p2),
                '${:.3}$'.format(p1),
                '${:.3}$'.format(rel_error)]
            cell_text.append(row_values)
            facecolors = ['#fed9b8', '#b8d6fe', '#fed9b8', '#b8d6fe'] if i % 2 == 0 else ['#b8d6fe', '#fed9b8', '#b8d6fe', '#fed9b8']
            cell_colors.append(facecolors)
        ## initialize table
        ax.table(
            cellText=cell_text,
            colLabels=column_labels,
            rowLabels=row_labels,
            loc='center',
            rowColours=['#cfb8fe' for i in range(len(row_labels))],
            colColours=['#cfb8fe' for j in range(len(column_labels))],
            cellColours=cell_colors,
            rowLoc='center',
            colLoc='center',
            cellLoc='center')
        s = 'Table of Orbital Periods (via {})'.format(self.ode_system_solver.method)
        if subtitle is not None:
            s += '\n{}'.format(subtitle)
        ax.set_title(
            s,
            fontsize=self.titlesize)
        ## save or show figure
        if save:
            savename = '' if subtitle is None else '{}-'.format(subtitle.replace(' ', '_'))
            savename += '{}-'.format(self.ode_system_solver.method.replace('-', '_').replace(' ', '_'))
            savename += 'Table-Orbital_Period'
            savename = savename.replace(' ', '-')
        else:
            savename = None
        self.display_image(fig, savename)

class NBodySimulation(SimulationConfiguration):

    def __init__(self, bodies, gravitational_constant=6.67e-11, name='', savedir=None, epsilon=1e-6):
        super().__init__(bodies, gravitational_constant, name, savedir, epsilon)

    def run_simulation(self, duration, time_step, duration_unit='year', step_unit='day', method='velocity-verlet', period_method=None, f_window=None):
        if method not in list(self.ode_system_solver.keys()):
            raise ValueError("invalid method: {}".format(method))
        if len(list(self.bodies.keys())) > 0:
            raise ValueError("simulation results (via {} method) are already initialized".format(method))
        cls = self.ode_system_solver[method]
        cls.solve_system_odes(duration, time_step, duration_unit, step_unit)
        self.update_bodies(cls)
        self.update_system_center_of_mass(cls)
        self._ode_system_solver = cls
        if period_method is not None:
            self.update_orbital_periods(
                method=period_method,
                f_window=f_window)

class SimulationComparison(VisualConfiguration):

    def __init__(self, isim, jsim, savedir=None):
        super().__init__(savedir)
        ndim, units, inclusive_names, exclusive_names, ibodies, jbodies = self.verify_simulations(isim, jsim)
        self.isim = isim
        self.jsim = jsim
        self.ndim = ndim
        self.units = units
        self.inclusive_names = inclusive_names
        self.exclusive_names = exclusive_names
        self.ibodies = ibodies
        self.jbodies = jbodies
        self.t = isim.ode_system_solver.t
        self.duration = isim.ode_system_solver.duration

    @staticmethod
    def verify_simulations(isim, jsim):
        ## verify time-steps
        if not np.array_equal(isim.ode_system_solver.t, jsim.ode_system_solver.t):
            raise ValueError("the time-steps for each simulation should be the same")
        if isim.ode_system_solver.duration != jsim.ode_system_solver.duration:
            raise ValueError("the duration of each simulation should be the same")
        ## get and verify number of dimensions
        ndim = np.array([sim.ode_system_solver.ndim for sim in (isim, jsim)])
        if np.all(np.diff(ndim) != 0):
            raise ValueError("the number of dimensions should be the same for each simulation")
        ## get and verify units
        if isim.ode_system_solver.units != jsim.ode_system_solver.units:
            raise ValueError("units should be the same for both simulations")
        units = dict(isim.ode_system_solver.units)
        units['phase'] = 'radian'
        ## get inclusive / exclusive names
        inames = list(isim.bodies.keys())
        jnames = frozenset(jsim.bodies.keys())
        inclusive_names = [None] + [name for name in inames if name in jnames]
        exclusive_names = [name for name in inames if name not in jnames]
        exclusive_names += [name for name in jnames if name not in frozenset(inames)]
        ## get bodies
        ibodies = isim.select_bodies(inclusive_names, selection_criteria='inclusive')
        jbodies = jsim.select_bodies(inclusive_names, selection_criteria='inclusive')
        if len(ibodies) != len(jbodies):
            raise ValueError("number of inclusive bodies should be the same for both simulations")
        return ndim, units, tuple(inclusive_names), tuple(exclusive_names), ibodies, jbodies

    def view_offset_distance(self, quantity, position_units=None, time_units=None, mass_units=None, energy_units=None, ref_frame='Center-of-Mass', names=None, selection_criteria='inclusive', subtitle=None, save=False, **kwargs):
        ## get time-scale
        if time_units is None:
            time_units = self.units['time']
        time_scale = self.get_scale('time', itime_units=self.units['time'], jtime_units=time_units)
        ## get duration
        duration = self.duration * time_scale
        if int(duration) == float(duration):
            duration = int(duration)
        ## get quantity_scale
        quantity_scale = self.get_scale(quantity,
            iposition_units=self.units['position'], jposition_units=position_units,
            itime_units=self.units['time'], jtime_units=time_units,
            imass_units=self.units['mass'], jmass_units=mass_units,
            ienergy_units=self.units['energy'], jenergy_units=energy_units)
        ## get bodies
        ibodies = self.isim.select_bodies(names, selection_criteria)
        jbodies = self.jsim.select_bodies(names, selection_criteria)
        if len(ibodies) != len(jbodies):
            raise ValueError("number of inclusive bodies should be the same for both simulations")
        ## initialize plot
        fig, ax = plt.subplots(**kwargs)
        x = self.t * time_scale
        if quantity == 'energy':
            labels = np.char.title(self.isim.scalar_quantities)
            if len(bodies) == 1:
                markers = ('.', '.', '.')
                facecolors = ('darkorange', 'steelblue', 'mediumviolet')
            else:
                markers = ('.', 'x', '_')
                facecolors = tuple([body['facecolor'] for body in ibodies])
            for ibody, jbody in zip(ibodies, jbodies):
                for qty, marker, label, facecolor in zip(self.scalar_quantities, markers, labels, facecolors):
                    y = ...
                    ax.scatter(x, y,
                        marker=marker,
                        label=r'${%s}_{%s}$' % (label, body['name']),
                        color=facecolor,
                        s=2)
        else:
            if quantity in self.isim.vector_quantities:
                for ibody, jbody in zip(ibodies, jbodies):
                    ivector = self.isim.get_vector_coordinates(quantity, ibody, ref_frame, quantity_scale)
                    jvector = self.jsim.get_vector_coordinates(quantity, jbody, ref_frame, quantity_scale)
                    y = self.isim.ode_system_solver.get_scalar(ivector - jvector, axis=0)
                    ax.scatter(x, y,
                        marker='.',
                        label=ibody['name'],
                        color=ibody['facecolor'],
                        s=2)
            elif quantity in self.scalar_quantities:
                for ibody, jbody in zip(ibodies, jbodies):
                    iscalar = ...
                    jscalar = ...
                    y = np.abs(iscalar - jscalar)
                    ax.scatter(x, y,
                        marker='.',
                        label=ibody['name'],
                        color=ibody['facecolor'],
                        s=2)
            else:
                raise ValueError("invalid quantity: {}".format(quantity))
        ## axis labels
        xlabel = self.get_label('time', position_units, time_units, mass_units, energy_units, include_quantity=True)
        ylabel = 'Offset of {}'.format(self.get_label(quantity, position_units, time_units, mass_units, energy_units, include_quantity=True))
        ax.set_xlabel(xlabel, fontsize=self.labelsize)
        ax.set_ylabel(ylabel, fontsize=self.labelsize)
        ## axis title
        title = 'Comparison of ${}$ {} Simulations'.format(duration, self.get_label('time', time_units=time_units))
        if subtitle:
            if not isinstance(subtitle, str):
                raise ValueError("invalid type(subtitle): {}".format(type(subtitle)))
            title += '\n{}'.format(subtitle)
        if ref_frame:
            title += '\nvia {} Frame'.format(ref_frame)
        ax.set_title(title, fontsize=self.titlesize)
        ## axis ticks
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax.tick_params(axis='both', labelsize=self.ticksize)
        ax.grid(color='k', linestyle=':', alpha=0.3)
        ## legend
        handles, labels = ax.get_legend_handles_labels()
        self.subview_legend(fig, ax, handles, labels, title=subtitle)
        ## save or show figure
        if save:

            savename = '' if subtitle is None else '{}-'.format(subtitle.replace(' ', '_'))
            # savename += '{}-'.format(self.ode_system_solver.method.replace('-', '_').replace(' ', '_'))
            savename += '{}_AND_{}-'.format(self.isim.ode_system_solver.method.replace('-', '_').replace(' ', '_'), self.jsim.ode_system_solver.method.replace('-', '_').replace(' ', '_'))
            savename += 'offset_{}-'.format(quantity.replace(' ', '_'))
            savename += 'ref_{}'.format(ref_frame.replace('-', '_').replace(' ', '_'))
            # savename += "{}_AND_{}".format(
            #     self.isim.ode_system_solver.method,
            #     self.jsim.ode_system_solver.method)
            savename = savename.replace(' ', '-')

            # savename = 'offset_{}'.format(quantity.replace(' ', '_'))
            # if ref_frame:
            #     savename += '_ref-{}'.format(ref_frame.replace(' ', '_'))
            # if subtitle:
            #     savename += '_{}'.format(subtitle.replace(' ', '_'))
            # savename += "{}_AND_{}".format(
            #     self.isim.ode_system_solver.method,
            #     self.jsim.ode_system_solver.method)
            # savename = savename.replace(' ', '_')

        else:
            savename = None
        self.display_image(fig, savename)

    def view_offset_angle(self, quantity, position_units=None, time_units=None, mass_units=None, energy_units=None, phase_units=None, ref_frame='Center-of-Mass', names=None, selection_criteria='inclusive', subtitle=None, save=False, **kwargs):
        if quantity not in self.isim.vector_quantities:
            raise ValueError("invalid quantity: {}; offset phase angle can only be computed from vector quantities".format(quantity))
        ## get time-scale
        if time_units is None:
            time_units = self.units['time']
        time_scale = self.get_scale('time', itime_units=self.units['time'], jtime_units=time_units)
        ## get duration
        duration = self.duration * time_scale
        if int(duration) == float(duration):
            duration = int(duration)
        ## get quantity_scale
        quantity_scale = self.get_scale(quantity,
            iposition_units=self.units['position'], jposition_units=position_units,
            itime_units=self.units['time'], jtime_units=time_units,
            imass_units=self.units['mass'], jmass_units=mass_units,
            ienergy_units=self.units['energy'], jenergy_units=energy_units)
        ## autocorrect units of phase
        phase_scale = self.get_scale('phase',
            iphase_units=self.units['phase'], jphase_units=phase_units)
        ## get bodies
        ibodies = self.isim.select_bodies(names, selection_criteria)
        jbodies = self.jsim.select_bodies(names, selection_criteria)
        if len(ibodies) != len(jbodies):
            raise ValueError("number of inclusive bodies should be the same for both simulations")
        ## initialize plot
        fig, ax = plt.subplots(**kwargs)
        x = self.t * time_scale
        for ibody, jbody in zip(ibodies, jbodies):
            ivector = self.isim.get_vector_coordinates(quantity, ibody, ref_frame, quantity_scale)
            jvector = self.jsim.get_vector_coordinates(quantity, jbody, ref_frame, quantity_scale)
            imag = self.isim.ode_system_solver.get_scalar(ivector, axis=0)
            jmag = self.jsim.ode_system_solver.get_scalar(jvector, axis=0)
            dproduct = np.nansum(ivector * jvector, axis=0)
            inv_theta = np.clip(dproduct / (imag * jmag), -1.0, 1.0) ## clip bounds (float precision error; 1.0002 --> 1 useful inside cosine)
            loc = np.isnan(inv_theta) ## divide by zero-magnitude vector ==> no offset
            inv_theta[loc] = 1 ## arccos(1) is zero
            y = np.arccos(inv_theta) * phase_scale
            ## plot
            ax.scatter(x, y,
                marker='.',
                label=ibody['name'],
                color=ibody['facecolor'],
                s=2)
        ## axis labels
        xlabel = self.get_label('time', position_units, time_units, mass_units, energy_units, phase_units, include_quantity=True)
        ylabel = 'Phase Offset ({})'.format(self.get_label('phase', position_units, time_units, mass_units, energy_units, phase_units))
        ax.set_xlabel(xlabel, fontsize=self.labelsize)
        ax.set_ylabel(ylabel, fontsize=self.labelsize)
        ## axis title
        title = 'Comparison of ${}$ {} Simulations'.format(duration, self.get_label('time', time_units=time_units))
        if subtitle:
            if not isinstance(subtitle, str):
                raise ValueError("invalid type(subtitle): {}".format(type(subtitle)))
            title += '\n{}'.format(subtitle)
        if ref_frame:
            title += '\nvia {} Frame'.format(ref_frame)
        ax.set_title(title, fontsize=self.titlesize)
        ## axis ticks
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax.tick_params(axis='both', labelsize=self.ticksize)
        ax.grid(color='k', linestyle=':', alpha=0.3)
        ## legend
        handles, labels = ax.get_legend_handles_labels()
        self.subview_legend(fig, ax, handles, labels, title=subtitle)
        ## save or show figure
        if save:

            savename = '' if subtitle is None else '{}-'.format(subtitle.replace(' ', '_'))
            savename += '{}_AND_{}-'.format(self.isim.ode_system_solver.method.replace('-', '_').replace(' ', '_'), self.jsim.ode_system_solver.method.replace('-', '_').replace(' ', '_'))
            savename += 'offset_phase_{}-'.format(quantity.replace(' ', '_'))
            savename += 'ref_{}'.format(ref_frame.replace('-', '_').replace(' ', '_'))
            # savename += "{}_AND_{}".format(
            #     self.isim.ode_system_solver.method,
            #     self.jsim.ode_system_solver.method)
            savename = savename.replace(' ', '-')

            # savename = 'offset_phase_{}'.format(quantity.replace(' ', '_'))
            # if ref_frame:
            #     savename += '_ref-{}'.format(ref_frame.replace(' ', '_'))
            # if subtitle:
            #     savename += '_{}'.format(subtitle.replace(' ', '_'))
            # savename += "{}_{}".format(
            #     self.isim.ode_system_solver.method,
            #     self.jsim.ode_system_solver.method)
            # savename = savename.replace(' ', '_')

        else:
            savename = None
        self.display_image(fig, savename)






##
