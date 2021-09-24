import numpy as np

class BaseQuantity():

    def __init__(self):
        super().__init__()
        self.position = {
            'unit' : ('meter', 'kilometer', 'mile', 'Astronomical Unit', 'Light Year', 'parsec'),
            'scale' : (1, 1000, 1609.34, 1.496e11, 9.461e15, 3.086e16),
            'label' : ('m', 'km', 'mi', 'AU', 'LY', 'pc')}
        self.time = {
            'unit' : ('second', 'minute', 'hour', 'day', 'year'),
            'scale' : (1, 60, 3600, 3600*24, 3600*24*365),
            'label' : ('s', 'min', 'hr', 'day', 'yr')}
        self.mass = {
            'unit' : ('kilogram', 'gram', 'lunar mass', 'earth mass', 'solar mass'),
            'scale' : (1, 1/1000, 7.3459e22, 5.9722e24, 1.988e30),
            'label' : ('kg', 'g', r'$M_{☽︎}$', r'$M_{⨁}$', r'$M_{☉}$')}
        self.energy = {
            'unit' : ('joule', 'kilojoule', 'megajoule', 'gigajoule', 'electron-volt', 'kiloelectron-volt', 'megaelectron-volt', 'gigaelectron-volt', 'erg'),
            'scale' : (1, 1000, 1000**2, 1000**3, 1.60217653e-19, 1.60217653e-16, 1.60217653e-13, 1.60217653e-10, 1e-7),
            'label' : ('J', 'kJ', 'MJ', 'GJ', 'eV', 'keV', 'MeV', 'GeV', 'erg')}
        self.phase = {
            'unit' : ('radian', 'degree', 'arcsecond'),
            'scale' : (1, np.pi/180, 1/206265),
            'label' : ('rad', '°', 'arcsec')}
        self.base_units = {
            'position' : self.position,
            'time' : self.time,
            'mass' : self.mass,
            'energy' : self.energy,
            'phase' : self.phase}

    def get_base_label(self, quantity, units):
        loc = self.base_units[quantity]['unit'].index(units)
        return self.base_units[quantity]['label'][loc]

    def get_base_conversion(self, quantity, original_units=None, prime_units=None):
        if (original_units is None) and (prime_units is None):
            return 1
        elif (original_units is None) or (prime_units is None):
            raise ValueError("cannot convert to/from units of type <None>")
        else:
            i = self.base_units[quantity]['unit'].index(original_units)
            j = self.base_units[quantity]['unit'].index(prime_units)
            scales = self.base_units[quantity]['scale']
            return scales[j] / scales[i]

class CompositeQuantity(BaseQuantity):

    def __init__(self):
        super().__init__()

    def get_velocity_label(self, position_units, time_units):
        labels = []
        for quantity, units in zip(('position', 'time'), (position_units, time_units)):
            i = self.base_units[quantity]['unit'].index(units)
            label = self.base_units[quantity]['label'][i]
            labels.append(label)
        return r'$\frac{%s}{%s}$' % (labels[0], labels[1])

    def get_acceleration_label(self, position_units, time_units):
        labels = []
        for quantity, units in zip(('position', 'time'), (position_units, time_units)):
            i = self.base_units[quantity]['unit'].index(units)
            label = self.base_units[quantity]['label'][i]
            labels.append(label)
        return r'$\frac{%s}{%s^2}$' % (labels[0], labels[1])

    def get_momentum_label(self, position_units, time_units, mass_units):
        velocity_label = self.get_velocity_label(position_units, time_units)
        mass_label = self.get_base_label('mass', mass_units)
        return '{} ⋅ {}'.format(mass_label, velocity_label)

    def get_angular_momentum_label(self, position_units, time_units, mass_units):
        mass_label = self.get_base_label('mass', mass_units)
        position_label = self.get_base_label('position', position_units)
        time_label = self.get_base_label('time', time_units)
        frac_label = r'$\frac{%s}{%s^2}$' % (position_label, time_label)
        return '{} ⋅ {}'.format(mass_label, frac_label)

    def get_force_label(self, position_units, time_units, mass_units):
        if (position_units == 'meter') and (time_units == 'second') and (mass_units == 'kilogram'):
            return 'N'
        else:
            mass_label = self.get_base_label('mass', mass_units)
            acc_label = self.get_acceleration_label(position_units, time_units)
            return '{} ⋅ {}'.format(mass_label, acc_label)

    def get_velocity_scale(self, iposition_units=None, jposition_units=None, itime_units=None, jtime_units=None):
        scale = 1
        for quantity, original_units, prime_units in zip(('position', 'time'), (iposition_units, itime_units), (jposition_units, jtime_units)):
            ratio = self.get_base_conversion(quantity, original_units, prime_units)
            if quantity == 'position':
                scale *= ratio
            else: # 'time'
                scale /= ratio
        return scale

    def get_acceleration_scale(self, iposition_units=None, jposition_units=None, itime_units=None, jtime_units=None):
        scale = 1
        for quantity, original_units, prime_units in zip(('position', 'time'), (iposition_units, itime_units), (jposition_units, jtime_units)):
            ratio = self.get_base_conversion(quantity, original_units, prime_units)
            if quantity == 'position':
                scale *= ratio
            else: # 'time'
                scale /= np.square(ratio)
        return scale

    def get_momentum_scale(self, iposition_units=None, jposition_units=None, itime_units=None, jtime_units=None, imass_units=None, jmass_units=None):
        scale = self.get_velocity_scale(iposition_units, jposition_units, itime_units, jtime_units)
        scale *= self.get_base_conversion('mass', imass_units, jmass_units)
        return scale

    def get_angular_momentum_scale(self, iposition_units=None, jposition_units=None, itime_units=None, jtime_units=None, imass_units=None, jmass_units=None):
        dm = self.get_base_conversion('mass', imass_units, jmass_units)
        dx = self.get_base_conversion('position', iposition_units, jposition_units)
        dt = self.get_base_conversion('time', itime_units, jtime_units)
        return dm * np.square(dx) / dt

    def get_force_scale(self, iposition_units=None, jposition_units=None, itime_units=None, jtime_units=None, imass_units=None, jmass_units=None):
        dm = self.get_base_conversion('mass', imass_units, jmass_units)
        da = self.get_acceleration_scale(iposition_units, jposition_units, itime_units, jtime_units)
        return dm * da

class Quantity(CompositeQuantity):

    def __init__(self):
        super().__init__()
        self.quantities = {
            'base' : ('position', 'time', 'mass', 'energy', 'phase'),
            'composite' : ('velocity', 'acceleration', 'momentum', 'angular velocity', 'angular acceleration', 'angular momentum')}

    def get_label(self, quantity, position_units=None, time_units=None, mass_units=None, energy_units=None, phase_units=None, include_quantity=False):
        if quantity == 'eccentricity':
            return 'Eccentricity'
        else:
            units_map = {
                'position' : position_units,
                'time' : time_units,
                'mass' : mass_units,
                'energy' : energy_units,
                'phase' : phase_units}
            if quantity in ('velocity', 'acceleration', 'angular velocity', 'angular acceleration'):
                s = 'get_{}_label'.format(quantity.replace(' ', '_'))
                f = getattr(self, s)
                label = f(position_units, time_units)
            elif quantity in ('momentum', 'angular momentum', 'force'):
                s = 'get_{}_label'.format(quantity.replace(' ', '_'))
                f = getattr(self, s)
                label = f(position_units, time_units, mass_units)
            elif ('energy' in quantity) or (quantity == 'lagrangian'): ## 'potential energy', 'kinetic energy', 'total energy'
                units = units_map['energy']
                label = self.get_base_label('energy', units)
            else:
                units = units_map[quantity]
                label = self.get_base_label(quantity, units)
            if include_quantity:
                return '{} ({})'.format(quantity.title(), label)
            else:
                return label

    def get_scale(self, quantity, iposition_units=None, jposition_units=None, itime_units=None, jtime_units=None, imass_units=None, jmass_units=None, ienergy_units=None, jenergy_units=None, iphase_units=None, jphase_units=None):
        units_map = {
            'position' : (iposition_units, jposition_units),
            'time' : (itime_units, jtime_units),
            'mass' : (imass_units, jmass_units),
            'energy' : (ienergy_units, jenergy_units),
            'phase' : (iphase_units, jphase_units)}
        if quantity in ('velocity', 'acceleration', 'angular velocity', 'angular acceleration'):
            f = getattr(self, 'get_{}_scale'.format(quantity.replace(' ', '_')))
            scale = f(iposition_units, jposition_units, itime_units, jtime_units)
        elif quantity in ('momentum', 'angular momentum', 'force'):
            f = getattr(self, 'get_{}_scale'.format(quantity.replace(' ', '_')))
            scale = f(iposition_units, jposition_units, itime_units, jtime_units, imass_units, jmass_units)
        elif ('energy' in quantity) or (quantity == 'lagrangian'): ## 'potential energy', 'kinetic energy', 'total energy'
            original_units, prime_units = units_map['energy']
            scale = self.get_base_conversion('energy', original_units, prime_units)
        elif quantity == 'eccentricity':
            scale = 1
        else:
            original_units, prime_units = units_map[quantity]
            scale = self.get_base_conversion(quantity, original_units, prime_units)
        return 1 / scale







#
