class MonteCarlo:
    def __init__(self, accept, n_trials, n_accept):
        self.accept = accept
        self.n_trials = n_trials
        self.n_accept = n_accept

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

        return self.accept


    def adjust_displacement(self, max_displacement):
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
        self.max_displacement = max_displacement

        acc_rate = float(self.n_accept) / float(self.n_trials)
        if (acc_rate < 0.38):
            self.max_displacement *= 0.8




        elif (acc_rate > 0.42):
            self.max_displacement *= 1.2

        self.n_trials = 0
        self.n_accept = 0

        return self.max_displacement, self.n_trials, self.n_accept