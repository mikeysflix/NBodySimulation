from astropy.coordinates import solar_system_ephemeris, get_body_barycentric_posvel
solar_system_ephemeris.set('jpl')
from astropy.time import Time
import datetime
from body_configuration import *

# class SolarSystemBodies(Quantity):
#
#     def __init__(self):
#         super().__init__()
#         self.ndim = 3
#         self.pos_scale = self.get_scale('position', iposition_units='Astronomical Unit', jposition_units='meter')
#         self.mass_scale = self.get_scale('mass', imass_units='earth mass', jmass_units='kilogram')
#         self.vel_scale = self.get_scale('velocity', iposition_units='mile', jposition_units='meter', itime_units='hour', jtime_units='second')
#
#     @property
#     def sun(self):
#         return Body(
#             name='Sun',
#             mass=self.get_scale('mass', imass_units='solar mass', jmass_units='kilogram'),
#             position=np.zeros(self.ndim),
#             facecolor='yellow')
#
#     @property
#     def mercury(self):
#         return Body(
#             name='Mercury',
#             mass=0.05227 * self.mass_scale,
#             position=np.array([0.4392, 0, 0]) * self.pos_scale,
#             velocity=np.array([0, 106000, 0]) * self.vel_scale,
#             facecolor='mediumvioletred',)
#
#     @property
#     def venus(self):
#         return Body(
#             name='Venus',
#             mass=0.815 * self.mass_scale,
#             position=np.array([0, 0.723, 0]) * self.pos_scale,
#             velocity=np.array([78300, 0, 0]) * self.vel_scale,
#             facecolor='darkorange',)
#
#     @property
#     def earth(self):
#         return Body(
#             name='Earth',
#             mass=self.mass_scale,
#             position=np.array([1, 0, 0]) * self.pos_scale,
#             velocity=np.array([0, 66600, 0]) * self.vel_scale, # np.array([0, 29.78, 0]),
#             facecolor='blue')
#
#     @property
#     def moon(self):
#         return Body(
#             name='Moon',
#             mass=self.get_scale('mass', imass_units='lunar mass', jmass_units='kilogram'),
#             position=np.array([1, 0.002682, 0]) * self.pos_scale,
#             velocity=np.array([2280, 0, 0]) * self.vel_scale,
#             facecolor='gray')
#
#     @property
#     def mars(self):
#         return Body(
#             name='Mars',
#             mass=0.1704 * self.mass_scale,
#             position=np.array([0, 1.653, 0]) * self.pos_scale,
#             velocity=np.array([53900, 0, 0]) * self.vel_scale,
#             facecolor='darkred')
#
#     @property
#     def jupiter(self):
#         return Body(
#             name='Jupiter',
#             mass=317.8 * self.mass_scale,
#             position=np.array([5.255, 0, 0]) * self.pos_scale,
#             velocity=np.array([0, 29200, 0]) * self.vel_scale,
#             facecolor='red')
#
#     @property
#     def saturn(self):
#         return Body(
#             name='Saturn',
#             mass=95.16 * self.mass_scale,
#             position=np.array([0, 10.04, 0]) * self.pos_scale,
#             velocity=np.array([21600, 0, 0]) * self.vel_scale,
#             facecolor='gold')
#
#     @property
#     def uranus(self):
#         return Body(
#             name='Uranus',
#             mass=14.54 * self.mass_scale,
#             position=np.array([19.83, 0, 0]) * self.pos_scale,
#             velocity=np.array([0, 15200, 0]) * self.vel_scale,
#             facecolor='limegreen')
#
#     @property
#     def neptune(self):
#         return Body(
#             name='Neptune',
#             mass=17.15 * self.mass_scale,
#             position=np.array([0, 29.93, 0]) * self.pos_scale,
#             velocity=np.array([12200, 0, 0]) * self.vel_scale,
#             facecolor='darkblue')
#
#     @property
#     def pluto(self):
#         return Body(
#             name='Pluto',
#             mass=0.002192 * self.mass_scale,
#             position=...,
#             velocity=...,
#             facecolor='k')
#
#     def get_bodies(self, names):
#         bodies = []
#         for name in names:
#             body = getattr(self, name.lower())
#             bodies.append(body)
#         return bodies

class SolarSystemBodies(Quantity):

    def __init__(self):
        """

        """
        super().__init__()
        self.mass_scale = self.get_scale('mass', imass_units='earth mass', jmass_units='kilogram')
        self.pos_scale = self.get_scale('position', iposition_units='kilometer', jposition_units='meter')
        self.vel_scale = self.get_scale('velocity', iposition_units='kilometer', jposition_units='meter', itime_units='day', jtime_units='second')
        self.body_mapping = {
            'Sun' : {
                'facecolor' : 'yellow',
                'mass' : self.get_scale('mass', imass_units='solar mass', jmass_units='kilogram')},
            'Mercury' : {
                'facecolor' : 'mediumvioletred',
                'mass' : 0.05227 * self.mass_scale},
            'Venus' : {
                'facecolor' : 'darkorange',
                'mass' : 0.815 * self.mass_scale},
            'Earth' : {
                'facecolor' : 'blue',
                'mass' : self.mass_scale},
            'Mars' : {
                'facecolor' : 'darkred',
                'mass' : 0.1704 * self.mass_scale},
            'Jupiter' : {
                'facecolor' : 'red',
                'mass' : 317.8 * self.mass_scale},
            'Saturn' : {
                'facecolor' : 'gold',
                'mass' : 95.16 * self.mass_scale},
            'Uranus' : {
                'facecolor' : 'limegreen',
                'mass' : 14.54 * self.mass_scale},
            'Neptune' : {
                'facecolor' : 'darkblue',
                'mass' : 17.15 * self.mass_scale}}

    def get_bodies(self, names, date_and_time=None):
        ## autocorrect date/time
        if date_and_time is None:
            date_and_time = datetime.datetime.now()
        dt = Time(
            val=date_and_time,
            format='datetime')
        ## verify names
        if isinstance(names, str):
            names = [names]
        elif not isinstance(names, (tuple, list, np.ndarray)):
            raise ValueError("invalid type(names): {}".format(names))
        available_names = list(self.body_mapping.keys())
        ## add bodies
        bodies = []
        for name in names:
            if name in available_names:
                pos_and_vel = get_body_barycentric_posvel(
                    body=name,
                    time=dt)
                body_attrs = self.body_mapping[name]
                position = pos_and_vel[0].xyz.value * self.pos_scale
                velocity = pos_and_vel[1].xyz.value * self.vel_scale
                body = Body(
                    name=name,
                    position=position,
                    velocity=velocity,
                    **body_attrs)
                bodies.append(body)
        if len(bodies) == 0:
            raise ValueError("zero bodies were selected")
        return bodies

class CircularTwoBodySystem(Quantity):

    def __init__(self):
        """

        """
        ## https://math.stackexchange.com/questions/1478440/defining-the-initial-conditions-for-a-planetary-motion-to-have-a-circular-orbit
        super().__init__()
        self.ndim = 3
        self.gravitational_constant = 6.67e-11

    def get_bodies(self):
        """

        """
        pos_primary = np.zeros(self.ndim)
        mass_primary = self.get_scale('mass', imass_units='solar mass', jmass_units='kilogram')
        _x = self.get_scale('position', iposition_units='Astronomical Unit', jposition_units='meter')
        pos_satellite = np.array([_x, 0, 0])
        _v = np.sqrt(self.gravitational_constant * mass_primary / _x)
        vel_satellite = np.array([0, _v, 0])
        if np.dot(pos_satellite, vel_satellite) != 0:
            raise ValueError("position vector and velocity vector (via initial condition of satellite) should be perpendicular")
        primary = Body(
            name='Primary',
            mass=mass_primary,
            position=np.zeros(self.ndim),
            facecolor='darkorange')
        satellite = Body(
            name='Secondary',
            mass=self.get_scale('mass', imass_units='earth mass', jmass_units='kilogram'),
            position=pos_satellite,
            velocity=vel_satellite,
            facecolor='steelblue')
        return [primary, satellite]

class BinaryStarSystemOrbit(Quantity):

    def __init__(self):
        """

        """
        super().__init__()

    def get_bodies(self):
        """

        """
        primary = ...
        secondary = ...
        return [primary, secondary]

class FigureEightOrbit(Quantity):

    def __init__(self):
        """

        """
        super().__init__()

    def get_bodies(self):
        """

        """
        first = ...
        second = ...
        third = ...
        return [first, second, third]




##
