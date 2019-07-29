m"""
Driver for Monte Carlo simulations of particles in fluid.
"""

def _parse_arguments(**kwargs):
    """Checks whether the user has provided any arguments directly to
    the main function through the kwargs keyword, and if not, attempt
    to find them from the command line arguments using argparse.
    Notes
    -----
    This method requires the following arguments:

    build_method, file_name, num_particles, simulation_cutoff, reduced_temperature,
    reduced_density, n_steps, freq, energy_function, max_displacement

    to be either all present in the kwargs dictionary, or given on the command line.
    Parameters
    ----------
    kwargs: dict of str, values, optional
        A dictionary of passed arguments (which may be empty), where the key is
        the name of the argument, and the value is the argument value.
    Returns
    -------
    argparse.Namespace
        A namespace containing the arguments which either were directly passed
        in as kwargs, or this which were taken as input from the command line.
        e.g. Namespace(build_method = 'random', file_name = None, num_particles = 100, simulation_cutoff = 3.0,
        reduced_temperature = 0.9, reduced_density = 0.9, n_steps = 1000000, freq = 1000, energy_function = 'LJ', max_displacement = 0.1)
    """
    arguments = None

    import argparse

    if not kwargs:

        # parsing command-line arguments
        parser = argparse.ArgumentParser(
            prog='mcfluid',
            description='Monte Carlo simulation of particles in fluid',
        )
        parser.add_argument(
            '--build_method',
            required=True,
            action='store',
            choices=['random', 'file'],
            nargs='?',
            help='the method to get the LJ fluid coordinates, either randomly generated or from a user-provided file. The default is "random".'
        )
        parser.add_argument(
            '--filename',
            required=True,
            action='store',
            nargs='?',
            type=str,
            help='the filename of the coordinates provided by users.')
        parser.add_argument(
            '--num_particles',
            required=True,
            action='store',
            nargs='?',
            type=str,
            help='the number of particles in the system. The default is 20.')
        parser.add_argument(
            '--simulation_cutoff',
            required=True,
            action='store',
            nargs='?',
            type=str,
            help='the number of particles in the system. The default is 20.')
        parser.add_argument(
            '--reduced_temperature',
            required=True,
            action='store',
            nargs='?',
            type=str,
            help='temperature in reduced units.')
        parser.add_argument(
            '--reduced_density',
            required=True,
            action='store',
            nargs='?',
            type=str,
            help='density in reduced units.')
        parser.add_argument(
            '--n_steps',
            required=True,
            action='store',
            nargs='?',
            type=str,
            help='total number of Monte Carlo simulation steps.')
        parser.add_argument(
            '--freq',
            required=True,
            action='store',
            nargs='?',
            type=str,
            help='output frequency.')
        parser.add_argument(
            '--energy_function',
            required=True,
            action='store',
            choices=['LJ', 'Buckingham', 'UnitlessLJ'],
            nargs='?',
            help='energy function used to calculate interactions among particles in fluid.'
        )
        parser.add_argument(
            '--max_displacement',
            required=True,
            action='store',
            nargs='?',
            type=str,
            help='initial size of displacement.')

        arguments = parser.parse_args()

    else:

        # Check to make sure that all of the required args are present
        # in the kwargs dictionary.

        # You might want to add your own validation logic here.
        assert 'build_method' in kwargs
        assert 'file_name' in kwargs
        assert 'num_particles' in kwargs
        assert 'simulation_cutoff' in kwargs
        assert 'reduced_temperature' in kwargs
        assert 'reduced_density' in kwargs
        assert 'n_steps' in kwargs
        assert 'freq' in kwargs
        assert 'energy_function' in kwargs
        assert 'max_displacement' in kwargs

        arguments = argparse.Namespace(**kwargs)

    return arguments


# main function
def mcfluid(**kwargs):
    """
    Perform Monte Carlo simulations on particles in fluid.
    Parameters
    ----------
    args: a Namespace object from argparse
        Information from parsing the command line
        e.g. Namespace(build_method = 'random', file_name = None, num_particles = 100, simulation_cutoff = 3.0,
        reduced_temperature = 0.9, reduced_density = 0.9, n_steps = 1000000, freq = 1000, energy_function = 'LJ', max_displacement = 0.1)
    Returns
    -------

    """

    # absolute import (with kinomodel installed)
    from mm_2019_sss_2 import system, energy, monte_carlo

    args = _parse_arguments(**kwargs)

    coordinates = generate_initial_state(method=args.build_method, num_particles=args.num_particles, box_length=args.box_length)

    if args.energy_function == "LJ":
        pass

    elif args.energy_function == "Buckinghamt":
        pass

    elif args.energy_function == "UnitlessLJ":
        pass
    else:
        raise Exception("Unknown energy function '{}'".format(args.feature))
        return