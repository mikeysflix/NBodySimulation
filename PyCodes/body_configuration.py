import copy
from units_configuration import *

class Body():

    def __init__(self, name, mass, position, velocity=None, facecolor='k'):
        super().__init__()
        self.name = name
        self.mass = mass
        self.facecolor = facecolor
        self._position = None
        self._velocity = None
        self.initialize_position_and_velocity(position, velocity)
        self._extra_info = dict()

    @property
    def position(self):
        return self._position

    @property
    def velocity(self):
        return self._velocity

    @property
    def extra_info(self):
        return self._extra_info

    @staticmethod
    def at_least_3d(args):
        if not isinstance(args, (tuple, list, np.ndarray)):
            raise ValueError("invalid type(args): {}".format(type(args)))
        size = len(args)
        if size < 3:
            if size == 0:
                raise ValueError("parameter cannot be empty")
            elif size == 1:
                args = (args[0], 0, 0)
            else:
                args = (args[0], args[1], 0)
        return np.array(args)

    def initialize_position_and_velocity(self, position, velocity=None):
        position = self.at_least_3d(position)
        if velocity is None:
            velocity = np.zeros(position.size)
        else:
            velocity = self.at_least_3d(velocity)
        if np.all(np.diff(np.array([position.size, velocity.size])) != 0):
            raise ValueError("position and velocity should have the same number of dimensions")
        self._position = position
        self._velocity = velocity

class MultiBodySystem(Quantity):

    def __init__(self, bodies, gravitational_constant, name):
        super().__init__()
        self._bodies = bodies
        self.gravitational_constant = gravitational_constant
        self.name = name
        self.nbodies = len(bodies)
        self.ndim = bodies[0].position.size
        standard_mu = []
        self.total_mass = 0
        self.reduced_mass = 1
        for body in bodies:
            self.total_mass += body.mass
            standard_mu.append(body.mass * gravitational_constant)
            self.reduced_mass *= body.mass
        self.reduced_mass /= self.total_mass
        self.standard_mu = np.array(standard_mu)
        self.units = {
            'position' : 'meter',
            'time' : 'second',
            'mass' : 'kilogram',
            'energy' : 'joule'}

    @property
    def bodies(self):
        return self._bodies




















##
