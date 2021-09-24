from backend_methods import *

## specify directory to save figures
savedir = '...' # None

## initialize viewing units
units = {
    'position_units' : 'Astronomical Unit',
    'time_units' : 'year',
    'mass_units' : 'earth mass',
    'energy_units' : 'joule'}

## select initial date/time
date_and_time = datetime.datetime(
    year=1900,
    month=1,
    day=1)
# date_and_time = None


## select simulations to run
run_circular_orbit = True # False
run_solar_system = True # False

if __name__ == "__main__":

    if run_circular_orbit:
        run_two_body_circular_orbits(
            duration=10,
            duration_unit='year',
            time_step=0.2,
            step_unit='day',
            methods=('euler-cromer', 'velocity-verlet'), #, 'yoshida'),
            savedir=savedir,
            **units)

    if run_solar_system:
        run_solar_system_orbits(
            duration=100,
            duration_unit='year',
            time_step=0.2,
            step_unit='day',
            methods=('euler-cromer', 'velocity-verlet'), #, 'yoshida'),
            date_and_time=date_and_time,
            savedir=savedir,
            **units)














##
