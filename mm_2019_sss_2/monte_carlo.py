class MonteCarlo:
    # TODO: which path?
    #  (1) to include other classes (System, EnergyFunction) into this class---how to specify the class name?
    #  (2) to import other modules in the same directory (System, EnergyFunction)
    def __init__(self, new_system : object = None, energy : object = None):
        # get parameters from the System object
        self.num_particles = new_system.n_particles
        self.coordinates = new_system.coordinates
        self.box_length = new_system.box_length
        # get functions from the EnergyFunction object
        self.calculate_initial_energy = energy.calculate_initial_energy
        self.calculate_pair_energy = energy.calculate_pair_energy
        self.calculate_tail_correction = energy.calculate_tail_correction

    def _parse_arguments(self, **kwargs):
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

        self.arguments = None

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
                defalte='random',
                help='the method to get the LJ fluid coordinates, either randomly generated or from a user-provided file. The default is "random".'
            )
            parser.add_argument(
                '--filename',
                required=True,
                action='store',
                nargs='?',
                type=str,
                defalte=None,
                help='the filename of the coordinates provided by users.')
            parser.add_argument(
                '--num_particles',
                required=True,
                action='store',
                nargs='?',
                type=int,
                defalte=20,
                help='the number of particles in the system. The default is 20.')
            parser.add_argument(
                '--simulation_cutoff',
                required=True,
                action='store',
                nargs='?',
                type=float,
                defalte=3.0,
                help='the number of particles in the system. The default is 20.')
            parser.add_argument(
                '--reduced_temperature',
                required=True,
                action='store',
                nargs='?',
                type=float,
                defalte=0.9,
                help='temperature in reduced units.')
            parser.add_argument(
                '--reduced_density',
                required=True,
                action='store',
                nargs='?',
                type=float,
                defalte=0.9,
                help='density in reduced units.')
            parser.add_argument(
                '--n_steps',
                required=True,
                action='store',
                nargs='?',
                type=int,
                defalte=1000000,
                help='total number of Monte Carlo simulation steps.')
            parser.add_argument(
                '--freq',
                required=True,
                action='store',
                nargs='?',
                type=int,
                defalte=1000,
                help='output frequency.')
            parser.add_argument(
                '--energy_function',
                required=True,
                action='store',
                choices=['LJ', 'Buckingham', 'UnitlessLJ'],
                nargs='?',
                defalte='LJ',
                help='energy function used to calculate interactions among particles in fluid.'
            )
            parser.add_argument(
                '--max_displacement',
                required=True,
                action='store',
                nargs='?',
                type=float,
                defalte=0.1,
                help='initial size of displacement.')

            parser.add_argument(
                'tune_displacement',
                required=True,
                action='store',
                choices=['LJ', 'Buckingham', 'UnitlessLJ'],
                nargs='?',
                type=bool,
                defalte=True,
                help='whether to adjust the size of displacement based on accept rate.'
            )

            self.arguments = parser.parse_args()

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
            assert 'tune_displacement' in kwargs

            self.arguments = argparse.Namespace(**kwargs)

        return self.arguments

    def accept_or_reject(self, delta_e, beta):
        """Accept or reject a move based on the energy difference and system \
             temperature.

        This function uses a random numbers to adjust the acceptance criteria.

        Parameters
        ----------
        delta_e : float
            The difference between the proposed and current energies.
        beta : float
            The inverse value of the reduced temperature.

        Returns
        -------
        accept : booleen
            Either a "True" or "False" to determine whether to reject the trial.
        """
        import numpy as np

        if delta_e < 0.0:
            accept = True
        else:
            random_number = np.random.rand(1)
            p_acc = np.exp(-beta * delta_e)

            if random_number < p_acc:
                accept = True
            else:
                accept = False

        return accept


    def adjust_displacement(self, max_displacement, n_accept, n_trials):
        """Change the acceptance criteria to get the desired rate.

        When the acceptance rate is too high, the maximum displacement is adjusted \
             to be higher.
        When the acceptance rate is too low, the maximum displacement is \
             adjusted lower.

        Parameters
        ----------
        n_trials : integer
            The number of trials that have been performed when the function is \
                 initiated.
        n_accept : integer
            The current number of accepted trials when the function is initiated.
        max_displacement : float
            The specified maximum value for the displacement of the trial.

        Returns
        -------
        max_displacement : float
            The adjusted displacement based on the acceptance rate.
        n_trials : integer, 0
            The new number of trials.
        n_accept : integer, 0
            The new number of trials.
        """

        acc_rate = float(n_accept) / float(n_trials)
        if (acc_rate < 0.38):
            max_displacement *= 0.8

        elif (acc_rate > 0.42):
            max_displacement *= 1.2

        n_trials = 0
        n_accept = 0

        return max_displacement, n_trials, n_accept

    ## THE MAIN FUNCTION
    def run_simulation(self):
        """This is the main function of the class that runs Monte Carlo simulation on particles in fluid."""

        import numpy as np

        args = self._parse_arguments(**kwargs)

        # set the initial total pair energy between particles in the system
        total_pair_energy = self.calculate_initial_energy(self.coordinates, self.box_length, args.simulation_cutoff)
        # TODO: need to add a function that explicitly calculates the tail correction
        tail_correction = self.calculate_tail_correction(self.coordinates, self.box_length, args.simulation_cutoff)

        # set up an array to store energy values
        energy_array = np.zeros(args.n_steps)

        # start the Monte Carlo iterations
        n_trials = 0
        for i_step in range(args.n_steps):
            n_trials += 1
            i_particle = np.random.randint(self.num_particles)
            random_displacement = (2.0 * np.random.rand(3) - 1.0) * max_displacement
            current_energy = self.calculate_pair_energy(self.coordinates, self.box_length, args.simulation_cutoff, i_particle)
            proposed_coordinates = self.coordinates.copy()
            proposed_coordinates[i_particle] += random_displacement
            proposed_coordinates -= self.box_length * np.round(proposed_coordinates / self.box_length)
            proposed_energy = self.calculate_pair_energy(self.coordinates, self.box_length, args.simulation_cutoff, i_particle)
            delta_e = proposed_energy - current_energy
            beta = 1.0 / args.reduced_temperature
            accept = self.accept_or_reject(delta_e, beta)
            if accept:
                total_pair_energy += delta_e
                n_accept += 1
            self.coordinates[i_particle] += random_displacement
            total_energy = (total_pair_energy + tail_correction) / self.num_particles
            energy_array[i_step] = total_energy
            if np.mod(i_step + 1, freq) == 0:
                print(i_step + 1, energy_array[i_step])
            if args.tune_displacement:
                max_displacement, n_trials, n_accept = self.adjust_displacement(n_trials, n_accept, max_displacement)
        return
