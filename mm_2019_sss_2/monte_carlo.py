class MonteCarlo:
    def __init__(self, system_setup: object = None,
                 energy_functions : object = None):
        self.num_particles = system_setup.n_particles
        self.coordinates = system_setup.coordinates
        self.box_length = system_setup.box_length

    def _accept_or_reject_(self, delta_e: float, beta: float):
        """Accept or reject a move based on the energy difference and system
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
        accept : boolean
            Either a "True" or "False" to determine whether to reject the trial.
        """
        import numpy as np

        self.delta_e = delta_e
        self.beta = beta

        if self.delta_e < 0.0:
            self.accept = True
        else:
            random_number = np.random.rand(1)
            p_acc = np.exp(-self.beta * self.delta_e)

            if random_number < p_acc:
                self.accept = True
            else:
                self.accept = False


    def adjust_displacement(self, max_displacement, n_accept, n_trials):
        """Change the acceptance criteria to get the desired rate.

        When the acceptance rate is too high, the maximum displacement is
        adjusted to be higher. When the acceptance rate is too low, the
        maximum displacement is adjusted lower.

        Parameters
        ----------
        n_trials : integer
            The number of trials that have been performed when the function is
            initiated.
        n_accept : integer
            The current number of accepted trials when the function is
            initiated.
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
        self.max_displacement = max_displacement
        self.n_accept = n_accept
        self.n_trials = n_trials

        acc_rate = float(self.n_accept) / float(self.n_trials)
        if (acc_rate < 0.38):
            self.max_displacement *= 0.8

        elif (acc_rate > 0.42):
            self.max_displacement *= 1.2

        self.n_trials = 0
        self.n_accept = 0

        return self.max_displacement, self.n_trials, self.n_accept

    ## EQUIVALENT TO THE MAIN FUNCTION
    def run_simulation(self):
        self.total_pair_energy = self.energy_functions.calculate_initial_energy(self.coordinates, self.box_length, simulation_cutoff)

        n_trials = 0
        for i_step in range(n_steps):
            n_trials += 1
            i_particle = np.random.randint(num_particles)
            random_displacement = (2.0 * np.random.rand(3) - 1.0) * max_displacement
            current_energy = energy_functions.calculate_pair_energy(self.coordinates, self.box_length, self.simulation_cutoff, self.i_particle)
            proposed_coordinates = coordinates.copy()
            proposed_coordinates[i_particle] += random_displacement
            proposed_coordinates -= box_length * np.round(proposed_coordinates / box_length)
            proposed_energy = energy_functions.calculate_pair_energy(self.coordinates, self.box_length, self.simulation_cutoff, self.i_particle)
            delta_e = proposed_energy - current_energy
            accept = accept_or_reject(delta_e, beta)
            if accept:
                total_pair_energy += delta_e
                n_accept += 1
            coordinates[i_particle] += random_displacement
            total_energy = (total_pair_energy + tail_correction) / num_particles
            energy_array[i_step] = total_energy
            if np.mod(i_step + 1, freq) == 0:
                print(i_step + 1, energy_array[i_step])
            if tune_displacement:
                max_displacement, n_trials, n_accept = adjust_displacement(n_trials, n_accept, max_displacement)
    return
