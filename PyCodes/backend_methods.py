from simulation_configuration import *
from sample_data import *

def run_two_body_circular_orbits(duration, duration_unit, time_step, step_unit, methods, savedir=None, **units):

    ## identify simulation results by subtitle
    subtitle = 'Circular Orbit'

    ## verify integration methods
    if isinstance(methods, str):
        methods = [methods]
    elif not isinstance(methods, (tuple, list, np.ndarray)):
        raise ValueError("invalid type(methods): {}".format(type(methods)))

    ## run simulation
    bodies = CircularTwoBodySystem().get_bodies()

    ## select method to solve ode system
    for method in methods:
        ## initialize simulation
        sim = NBodySimulation(
            bodies,
            name='Circular Orbit',
            savedir=savedir)
        ## run simulation
        sim.run_simulation(
            duration,
            time_step,
            duration_unit,
            step_unit,
            method)

        ## view vector quantities
        for quantity in sim.vector_quantities:
            ## 2-D (xy) and ## 3-D
            for dims in [(0,1), (0,1,2)]:
                for ref_frame in ('Center-of-Mass', 'Primary', 'Secondary'):
                    sim.view_vector_quantity(
                        quantity=quantity,
                        ref_frame=ref_frame,
                        subtitle=subtitle,
                        dims=dims,
                        figsize=(12, 7),
                        save=True,
                        **units)

        ## view animation of vector position
        for extension in ('.mp4',): # '.gif'):
            for dims in [(0,1), (0,1,2)]:
                for ref_frame in ('Center-of-Mass', 'Primary', 'Secondary'):
                    sim.view_vector_quantity_animation(
                        quantity='position',
                        ref_frame=ref_frame,
                        subtitle=subtitle,
                        dims=dims,
                        # elev=None,
                        # azim=None,
                        save=True,
                        extension=extension,
                        fps=60,
                        **units)

        ## view scalar quantities
        for quantity in list(sim.scalar_quantities) + list(sim.vector_quantities):
            for ref_frame in ('Center-of-Mass', 'Primary', 'Secondary'):
                sim.view_scalar_quantity(
                    quantity=quantity,
                    ref_frame=ref_frame,
                    subtitle=subtitle,
                    figsize=(12,7),
                    save=True,
                    **units)

def run_solar_system_orbits(duration, duration_unit, time_step, step_unit, methods, date_and_time=None, savedir=None, **units):

    ## initialize true orbital periods
    solar_system_periods = {
        'Sun' : 7.45e15,
        'Mercury' : 7.600544e6,
        'Venus' : 1.9414149e7,
        'Earth' : 3.1558149e7,
        'Mars' : 5.9355036e7,
        'Jupiter' : 3.7435566e8,
        'Saturn' : 9.2929236e8,
        'Uranus' : 2.65137e9,
        'Neptune' : 5.2004186e9}

    ## identify simulation results by subtitle
    base_subtitle = 'Solar System'
    offset_subtitle = 'Solar System (with and without Neptune)'

    ## verify integration methods
    if isinstance(methods, str):
        methods = [methods]
    elif not isinstance(methods, (tuple, list, np.ndarray)):
        raise ValueError("invalid type(methods): {}".format(type(methods)))

    ## initialize simulations
    for method in methods:

        ## run simulation (with Neptune)
        bodies_with_neptune = SolarSystemBodies().get_bodies(
            names=(
                'Sun',
                'Mercury',
                'Venus',
                'Earth',
                'Mars',
                'Jupiter',
                'Saturn',
                'Uranus',
                'Neptune'),
            date_and_time=date_and_time)
        nbody_with_neptune = NBodySimulation(
            bodies_with_neptune,
            name='Solar System',
            savedir=savedir)
        nbody_with_neptune.run_simulation(
            duration,
            time_step,
            duration_unit,
            step_unit,
            method,
            period_method='fourier',
            f_window=np.hamming)

        ## view table of orbital periods
        nbody_with_neptune.view_orbital_periods_table(
            true_periods=solar_system_periods,
            subtitle=base_subtitle,
            save=True,
            figsize=(7, 7))

        ## run simulation (without Neptune)
        bodies_without_neptune = SolarSystemBodies().get_bodies(
            names=(
                'Sun',
                'Mercury',
                'Venus',
                'Earth',
                'Mars',
                'Jupiter',
                'Saturn',
                'Uranus'),
            date_and_time=date_and_time)
        nbody_without_neptune = NBodySimulation(
            bodies_without_neptune,
            name='Solar System',
            savedir=savedir)
        nbody_without_neptune.run_simulation(
            duration,
            time_step,
            duration_unit,
            step_unit,
            method)

        ## initialize comparison of both simulations
        nbody_comparison = SimulationComparison(
            nbody_with_neptune,
            nbody_without_neptune,
            savedir=savedir)

        ## view offset phase angle via vector quantities
        for quantity in nbody_comparison.isim.vector_quantities:
            for ref_frame in ('Center-of-Mass', 'Sun', 'Earth', 'Uranus'):
                nbody_comparison.view_offset_angle(
                    quantity,
                    phase_units='arcsecond',
                    ref_frame=ref_frame,
                    names='Uranus',
                    subtitle=offset_subtitle,
                    figsize=(12,7),
                    save=True,
                    **units)

        ## view offset magnitude of vector/scalar quantities
        for quantity in list(nbody_comparison.isim.vector_quantities): # + ['energy']:
            nbody_comparison.view_offset_distance(
                quantity,
                ref_frame=ref_frame,
                names='Uranus',
                subtitle=offset_subtitle,
                figsize=(12,7),
                save=True,
                **units)

        ## view vector quantities of simulation (with neptune)
        for quantity in nbody_comparison.isim.vector_quantities:
            for dims in [(0,1), (0,1,2)]:
                for ref_frame in ('Center-of-Mass', 'Sun', 'Earth', 'Uranus'):
                    nbody_with_neptune.view_vector_quantity(
                        quantity,
                        ref_frame=ref_frame,
                        subtitle=base_subtitle,
                        dims=dims,
                        figsize=(12,7),
                        save=True,
                        **units)

        ## view animation of vector position
        for extension in ('.mp4',): # '.gif'):
            for dims in [(0,1), (0,1,2)]:
                for ref_frame in ('Center-of-Mass', 'Sun', 'Earth', 'Uranus'):
                    nbody_with_neptune.view_vector_quantity_animation(
                        quantity='position',
                        ref_frame=ref_frame,
                        subtitle=base_subtitle,
                        dims=dims,
                        # elev=None,
                        # azim=None,
                        save=True,
                        extension=extension,
                        fps=60,
                        **units)

        ## view scalar quantities of simulation (with neptune)
        for quantity in list(nbody_with_neptune.scalar_quantities) + list(nbody_with_neptune.vector_quantities):
            for ref_frame in ('Center-of-Mass', 'Sun', 'Earth', 'Uranus'):
                nbody_with_neptune.view_scalar_quantity(
                    quantity,
                    ref_frame=ref_frame,
                    subtitle=base_subtitle,
                    figsize=(12,7),
                    save=True,
                    **units)





##
