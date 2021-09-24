from body_configuration import *

class SolverConfiguration(MultiBodySystem):

    def __init__(self, bodies, gravitational_constant, name):
        super().__init__(bodies, gravitational_constant, name)
        self._t = None
        self._dt = None
        self._duration = None
        self._initial_condition = dict()
        self._solution = dict()

    @property
    def t(self):
        return self._t

    @property
    def dt(self):
        return self._dt

    @property
    def duration(self):
        return self._duration

    @property
    def initial_condition(self):
        return self._initial_condition

    @property
    def solution(self):
        return self._solution

    @staticmethod
    def get_scalar(vector, axis): ## CARTESIAN
        return np.sqrt(np.nansum(np.square(vector), axis=axis))

    def compute_accelerations(self, positions):
        ## FIX ME !!
        # ... vectorize calculation of acceleration
        ## FIX ME !!
        acceleration = []
        for i in range(self.nbodies):
            dacc = 0
            for j in range(self.nbodies):
                if i != j:
                    dpos = positions[i, :] - positions[j, :]
                    r = np.sum(np.square(dpos))
                    dacc -= self.standard_mu[j] * dpos / np.sqrt(r**3)
            acceleration.append(dacc)
        return np.array(acceleration)

    def initialize_equalized_time_steps(self, duration, time_step, duration_unit='year', step_unit='day'):
        step_scale = self.get_scale('time', itime_units=step_unit, jtime_units=self.units['time'])
        dur_scale = self.get_scale('time', itime_units=duration_unit, jtime_units=self.units['time'])
        self._dt = time_step * step_scale
        self._duration = dur_scale * duration
        nsteps = int(np.floor(self.duration / self.dt)) # +1
        self._t = np.arange(nsteps) * self.dt

    def initialize_initial_condition(self):
        pos, vel = [], []
        for body in self.bodies:
            pos.append(body.position)
            vel.append(body.velocity)
        pos, vel = np.array(pos), np.array(vel)
        acc = self.compute_accelerations(pos)
        initial_condition = {
            'position' : np.array(pos),
            'velocity' : np.array(vel),
            'acceleration' : np.array(acc)}
        self._initial_condition.update(initial_condition)

    def update_initial_condition_by_halfstep(self, parameter):
        if parameter == 'position':
            self._initial_condition['position'] += self.initial_condition['velocity'] * self.dt / 2
        elif parameter == 'velocity':
            self._initial_condition['velocity'] += self.initial_condition['acceleration'] * self.dt / 2
        else:
            raise ValueError("invalid parameter: {}".format(parameter))

class ModifiedEulerIntegrator(SolverConfiguration):

    def __init__(self, bodies, gravitational_constant, name):
        super().__init__(bodies, gravitational_constant, name)
        self.method = 'Euler-Cromer'
        self.is_symplectic = True
        self.order = 1

    def solve_system_odes(self, duration, time_step, duration_unit, step_unit):
        ## initialize time-steps with unit-conversions
        self.initialize_equalized_time_steps(duration, time_step, duration_unit, step_unit)
        ## initialize i=0 step
        self.initialize_initial_condition()
        ri = np.copy(self.initial_condition['position'])
        vi = np.copy(self.initial_condition['velocity'])
        ai = np.copy(self.initial_condition['acceleration'])
        ## initialize solution-space
        pos = np.zeros((self.nbodies, self.ndim, self.t.size))
        vel = np.zeros((self.nbodies, self.ndim, self.t.size))
        acc = np.zeros((self.nbodies, self.ndim, self.t.size))
        pos[:, :, 0] = ri
        vel[:, :, 0] = vi
        acc[:, :, 0] = ai
        ## iterate time-steps
        for i in range(1, self.t.size):
            ak = self.compute_accelerations(ri)
            vk = vi + ak * self.dt
            rk = ri + vk * self.dt
            pos[:, :, i] = rk
            vel[:, :, i] = vk
            acc[:, :, i] = ak
            ri, vi, ai = rk, vk, ak
        ## update solution
        self._solution.update({
            'position' : pos,
            'velocity' : vel,
            'acceleration' : acc})

class VerletLeapFrogIntegrator(SolverConfiguration):

    def __init__(self, bodies, gravitational_constant, name):
        super().__init__(bodies, gravitational_constant, name)
        self.method = 'Velocity-Verlet (Leap-Frog)'
        self.is_symplectic = True
        self.order = 2

    #
    # i --> i, current-step
    # j --> i + 1/2, half-step
    # k --> i + 1, full-step
    #

    def kick(self, vi, ai):
        vj = vi + ai * self.dt / 2
        return vj

    def drift(self, ri, vj):
        rk = ri + vj * self.dt
        ak = self.compute_accelerations(rk)
        vk = vj + ak * self.dt / 2
        return rk, vk, ak

    def solve_system_odes(self, duration, time_step, duration_unit, step_unit):
        ## initialize time-steps with unit-conversions
        self.initialize_equalized_time_steps(duration, time_step, duration_unit, step_unit)
        ## initialize i=0 step
        self.initialize_initial_condition()
        ri = np.copy(self.initial_condition['position'])
        vi = np.copy(self.initial_condition['velocity'])
        ai = np.copy(self.initial_condition['acceleration'])
        ## initialize i=1/2 step
        self.update_initial_condition_by_halfstep('position')
        vj = np.copy(self.initial_condition['velocity'])
        ## initialize solution-space
        pos = np.zeros((self.nbodies, self.ndim, self.t.size))
        vel = np.zeros((self.nbodies, self.ndim, self.t.size))
        acc = np.zeros((self.nbodies, self.ndim, self.t.size))
        pos[:, :, 0] = ri
        vel[:, :, 0] = vi
        acc[:, :, 0] = ai
        ## iterate time-steps (kick + drift)
        for i in range(1, self.t.size):
            vj = self.kick(vi, ai)
            rk, vk, ak = self.drift(ri, vj)
            pos[:, :, i] = rk
            vel[:, :, i] = vk
            acc[:, :, i] = ak
            ri, vi, ai = rk, vk, ak
        ## update solution
        self._solution.update({
            'position' : pos,
            'velocity' : vel,
            'acceleration' : acc})

class YoshidaLeapFrogIntegrator(SolverConfiguration):

    def __init__(self, bodies, gravitational_constant, name):
        super().__init__(bodies, gravitational_constant, name)
        self.method = 'Yoshida (Leap-Frog)'
        self.is_symplectic = True
        self.order = 4

    def solve_system_odes(self, duration, time_step, duration_unit, step_unit):
        ## initialize time-steps with unit-conversions
        self.initialize_equalized_time_steps(duration, time_step, duration_unit, step_unit)
        ## initialize i=0 step
        self.initialize_initial_condition()
        ri = np.copy(self.initial_condition['position'])
        vi = np.copy(self.initial_condition['velocity'])
        ai = np.copy(self.initial_condition['acceleration'])
        ## initialize solution-space
        pos = np.zeros((self.nbodies, self.ndim, self.t.size))
        vel = np.zeros((self.nbodies, self.ndim, self.t.size))
        acc = np.zeros((self.nbodies, self.ndim, self.t.size))
        pos[:, :, 0] = ri
        vel[:, :, 0] = vi
        acc[:, :, 0] = ai
        ## weight interval between steps [i], [k=i+1]
        w1 = 1 / (2 -2**(1/3))
        w0 = 2**(1/3) * w1
        c1 = c4 = w1/2
        c2 = c3 = (w0 + w1) / 2
        d1 = d3 = w1
        d2 = w0
        ## iterate time-steps (kick + drift)
        for i in range(1, self.t.size):
            ## use weighted sub-intervals
            r1 = ri + c1 * vi * self.dt
            a1 = self.compute_accelerations(r1)
            v1 = vi + d1 * a1 * self.dt
            r2 = r1 + c2 * v1 * self.dt
            a2 = self.compute_accelerations(r2)
            v2 = v1 + d2 * a2 * self.dt
            r3 = r2 + c3 * v2 * self.dt
            a3 = self.compute_accelerations(r3)
            v3 = v2 + d3 * a3 * self.dt
            ## aggregate
            rk = r3 + c4 * v3 * self.dt
            vk = v3
            ak = self.compute_accelerations(rk)
            ## store and re-loop
            pos[:, :, i] = rk
            vel[:, :, i] = vk
            acc[:, :, i] = ak
            ri, vi, ai = rk, vk, ak
        ## update solution
        self._solution.update({
            'position' : pos,
            'velocity' : vel,
            'acceleration' : acc})

class RungeKuttaFourIntegrator(SolverConfiguration):

    def __init__(self, bodies, gravitational_constant, name):
        super().__init__(bodies, gravitational_constant, name)
        self.method = 'RK-4'
        self.is_symplectic = False
        self.order = 4

        ## https://stackoverflow.com/questions/52334558/runge-kutta-4th-order-method-to-solve-second-order-odes
        ## https://math.stackexchange.com/questions/721076/help-with-using-the-runge-kutta-4th-order-method-on-a-system-of-2-first-order-od


    def solve_system_odes(self, duration, time_step, duration_unit, step_unit):
        ## initialize time-steps with unit-conversions
        self.initialize_equalized_time_steps(duration, time_step, duration_unit, step_unit)
        ## initialize i=0 step
        self.initialize_initial_condition()
        ri = np.copy(self.initial_condition['position'])
        vi = np.copy(self.initial_condition['velocity'])
        ai = np.copy(self.initial_condition['acceleration'])
        ## initialize solution-space
        pos = np.zeros((self.nbodies, self.ndim, self.t.size))
        vel = np.zeros((self.nbodies, self.ndim, self.t.size))
        acc = np.zeros((self.nbodies, self.ndim, self.t.size))
        pos[:, :, 0] = ri
        vel[:, :, 0] = vi
        acc[:, :, 0] = ai
        ## iterate time-steps
        for i in range(1, self.t.size):
            ka1 = self.compute_accelerations(ri) * self.dt
            kv1 = self.dt * ka1
            ka2 = ...
            kv2 = ...
            ka3 = ...
            kv3 = ...
            ka4 = ...
            kv4 = ...
            ak = self.dt * (k1 + 2*k2 + 2*k3 + k4) / 6
            vk = vi + ak * self.dt
            rk = ri + vk * self.dt
            pos[:, :, i] = rk
            vel[:, :, i] = vk
            acc[:, :, i] = ak
            ri, vi, ai = rk, vk, ak
        ## update solution
        self._solution.update({
            'position' : pos,
            'velocity' : vel,
            'acceleration' : acc})

class AdamsBashforthMoultonIntegrator(SolverConfiguration):

    def __init__(self, bodies, gravitational_constant, name):
        super().__init__(bodies, gravitational_constant, name)
        self.method = 'Adams-Bashforth-Moulton'
        self.is_symplectic = False
        self.order = ...

    def solve_system_odes(self, duration, time_step, duration_unit, step_unit):
        ...























##
